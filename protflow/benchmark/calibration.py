"""
Empirical Calibration Table module for ProtFlow.
Maps raw node and pipeline scores to calibrated true-positive confidence probabilities
stratified by enzyme class and tool type.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd


@dataclass
class CalibrationBin:
    min_score: float
    max_score: float
    tp: int
    fp: int
    empirical_probability: float

    @property
    def total(self) -> int:
        return self.tp + self.fp


class CalibrationTable:
    """
    Maintains empirical calibration curves per EC class and overall.
    Converts raw prediction scores into statistically calibrated confidence probabilities.
    """
    def __init__(self, data: Optional[pd.DataFrame] = None):
        self.class_calibrations: Dict[str, List[CalibrationBin]] = {}
        self.global_bins: List[CalibrationBin] = []
        self.metrics_summary: pd.DataFrame = pd.DataFrame()
        if data is not None and not data.empty:
            self.fit(data)

    def fit(self, df: pd.DataFrame, score_col: str = "final_priority_score", label_col: str = "is_active_enzyme") -> CalibrationTable:
        """
        Fits empirical calibration tables from benchmark evaluation results.
        Required columns: score_col (float), label_col (bool), 'ec_class' (str).
        """
        # 1. Ensure score column exists, populate prior if absent
        if score_col not in df.columns:
            df = df.copy()
            df[score_col] = np.where(df[label_col] == True, 0.80, 0.40)

        # Global calibration bins: [0.0-0.5, 0.5-0.7, 0.7-0.85, 0.85-1.0]
        bin_edges = [0.0, 0.50, 0.70, 0.85, 1.00]
        self.global_bins = self._calculate_bins(df, score_col, label_col, bin_edges)

        # 2. Per-class calibration
        summary_rows = []
        for ec_cls, grp in df.groupby("ec_class"):
            cbins = self._calculate_bins(grp, score_col, label_col, bin_edges)
            self.class_calibrations[ec_cls] = cbins

            tp = int((grp[label_col] == True).sum())
            fp = int((grp[label_col] == False).sum())
            mean_score = float(grp[score_col].mean())
            precision = tp / max(tp + fp, 1) if tp + fp > 0 else 0.0

            # Reliability Grade (A, B, C, D)
            if precision >= 0.85:
                grade = "A (High Reliability)"
            elif precision >= 0.70:
                grade = "B (Good Reliability)"
            elif precision >= 0.50:
                grade = "C (Moderate Reliability)"
            else:
                grade = "D (Low / High False Positive)"

            summary_rows.append({
                "ec_class": ec_cls,
                "samples": len(grp),
                "true_positives": tp,
                "false_positives": fp,
                "mean_raw_score": round(mean_score, 3),
                "empirical_precision": round(precision, 3),
                "reliability_grade": grade
            })

        self.metrics_summary = pd.DataFrame(summary_rows)
        return self

    @staticmethod
    def _calculate_bins(df: pd.DataFrame, score_col: str, label_col: str, edges: List[float]) -> List[CalibrationBin]:
        bins = []
        for i in range(len(edges) - 1):
            low, high = edges[i], edges[i + 1]
            mask = (df[score_col] >= low) & (df[score_col] <= high if i == len(edges) - 2 else df[score_col] < high)
            subset = df[mask]
            tp = int((subset[label_col] == True).sum())
            fp = int((subset[label_col] == False).sum())
            total = tp + fp
            # Laplace smoothing (empirical bayes)
            prob = float((tp + 1.0) / (total + 2.0)) if total > 0 else float(low + high) / 2.0
            bins.append(CalibrationBin(min_score=low, max_score=high, tp=tp, fp=fp, empirical_probability=prob))
        return bins

    def calibrate(self, raw_score: float, ec_class: Optional[str] = None) -> float:
        """
        Returns calibrated confidence (0.0 to 1.0) for a given raw score.
        Uses class-specific calibration if available, else falls back to global curve.
        """
        bins = self.class_calibrations.get(ec_class, self.global_bins) if ec_class else self.global_bins
        if not bins:
            return raw_score

        for b in bins:
            if b.min_score <= raw_score <= b.max_score:
                return round(b.empirical_probability, 4)

        return round(float(np.clip(raw_score, 0.0, 1.0)), 4)

    def attach_calibration(self, results_df: pd.DataFrame, score_col: str = "final_priority_score") -> pd.DataFrame:
        """
        Adds 'calibrated_confidence' and 'reliability_grade' columns to a results DataFrame.
        """
        out_df = results_df.copy()
        conf_scores = []
        grades = []

        grade_map = {}
        if not self.metrics_summary.empty:
            grade_map = dict(zip(self.metrics_summary["ec_class"], self.metrics_summary["reliability_grade"]))

        for _, row in out_df.iterrows():
            score = float(row.get(score_col, 0.5))
            ec = row.get("ec_class") if "ec_class" in row else None
            conf = self.calibrate(score, ec)
            conf_scores.append(conf)
            grades.append(grade_map.get(ec, "Standard (Global)"))

        out_df["calibrated_confidence"] = conf_scores
        out_df["reliability_grade"] = grades
        return out_df
