"""
Pipeline Optimization & Grid Sweeper module for ProtFlow.
Evaluates multiple pipeline configurations against labeled enzyme benchmarks
to discover the Pareto-optimal parameter set.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import logging
import time
import pandas as pd
import numpy as np

logger = logging.getLogger("protflow.benchmark.optimizer")

from protflow.pipeline import Pipeline
from protflow.runner import Runner
from protflow.types import Sequence, SequenceCollection, Table
from protflow.benchmark.dataset import load_benchmark_dataframe, SWISS_PROT_BENCHMARK
from protflow.benchmark.calibration import CalibrationTable


@dataclass
class VariantConfig:
    name: str
    pocket_tool: str = "p2rank"  # p2rank or fpocket
    min_plddt: float = 70.0
    weight_quality: float = 0.30
    weight_pocket: float = 0.30
    weight_solubility: float = 0.20
    weight_novelty: float = 0.20


@dataclass
class VariantResult:
    variant: VariantConfig
    accuracy: float
    precision: float
    recall: float
    f1_score: float
    avg_latency_ms: float
    predictions_df: pd.DataFrame


class PipelineOptimizer:
    """
    Evaluates and compares pipeline variants against Swiss-Prot benchmark ground-truth.
    """
    def __init__(self, base_pipeline_yaml: Path, runner: Optional[Runner] = None):
        self.base_pipeline_yaml = Path(base_pipeline_yaml)
        self.runner = runner or Runner(enable_cache=True)

    def evaluate_variant(self, variant: VariantConfig, benchmark_df: Optional[pd.DataFrame] = None) -> VariantResult:
        df = benchmark_df if benchmark_df is not None else load_benchmark_dataframe()

        # Load pipeline and apply variant config
        pipeline = Pipeline.from_yaml(self.base_pipeline_yaml, runner=self.runner)

        # 1. Update pocket tool
        pocket_step = next((s for s in pipeline.steps if s.name == "pocket_finding"), None)
        if pocket_step:
            manifest_name = f"manifests/05_pocket_{variant.pocket_tool}.yaml"
            pocket_step.manifest_path = manifest_name

        # 2. Update plddt threshold
        plddt_step = next((s for s in pipeline.steps if s.name == "quality_assessment"), None)
        if plddt_step:
            plddt_step.parameters["min_plddt"] = variant.min_plddt

        # Run pipeline over all benchmark proteins
        pred_rows = []
        latencies = []

        for _, row in df.iterrows():
            seq_obj = Sequence(protein_id=row["protein_id"], sequence=row["sequence"], description=row["description"])
            t0 = time.perf_counter()

            # Execute pipeline steps
            try:
                coll = SequenceCollection(protein_id=seq_obj.protein_id, sequences={seq_obj.protein_id: seq_obj})
                state = pipeline.run({"protein": seq_obj, "proteins": coll})
                p_table = state["priority_table"].df
                raw_score = float(p_table["final_priority_score"].iloc[0])

                # Recalculate with variant specific weights if custom
                q = float(p_table["structure_quality"].iloc[0])
                pkt = float(p_table["pocket_score"].iloc[0])
                sol = float(p_table["solubility_est"].iloc[0])
                nov = float(p_table["novelty_score"].iloc[0])
                adj_score = (
                    variant.weight_quality * q +
                    variant.weight_pocket * pkt +
                    variant.weight_solubility * sol +
                    variant.weight_novelty * nov
                )

                # Classify candidate: threshold at 0.60
                is_predicted_active = adj_score >= 0.60
                latencies.append((time.perf_counter() - t0) * 1000.0)

            except Exception as e:
                logger.warning("Pipeline variant evaluation failed for protein '%s': %s", row["protein_id"], e)
                adj_score = 0.0
                is_predicted_active = False
                latencies.append(10.0)

            pred_rows.append({
                "protein_id": row["protein_id"],
                "ec_class": row["ec_class"],
                "is_active_enzyme": bool(row["is_active_enzyme"]),
                "final_priority_score": round(adj_score, 4),
                "predicted_active": is_predicted_active
            })

        preds_df = pd.DataFrame(pred_rows)

        # Calculate classification metrics
        y_true = preds_df["is_active_enzyme"].values
        y_pred = preds_df["predicted_active"].values

        tp = int(np.sum((y_true == True) & (y_pred == True)))
        fp = int(np.sum((y_true == False) & (y_pred == True)))
        fn = int(np.sum((y_true == True) & (y_pred == False)))
        tn = int(np.sum((y_true == False) & (y_pred == False)))

        accuracy = float((tp + tn) / max(len(y_true), 1))
        precision = float(tp / max(tp + fp, 1))
        recall = float(tp / max(tp + fn, 1))
        f1 = float(2 * (precision * recall) / max(precision + recall, 1e-6))
        avg_lat = float(np.mean(latencies)) if latencies else 0.0

        return VariantResult(
            variant=variant,
            accuracy=round(accuracy, 3),
            precision=round(precision, 3),
            recall=round(recall, 3),
            f1_score=round(f1, 3),
            avg_latency_ms=round(avg_lat, 1),
            predictions_df=preds_df
        )

    def sweep(self, variants: List[VariantConfig]) -> Tuple[VariantResult, pd.DataFrame, CalibrationTable]:
        """
        Sweeps across candidate variants, ranks them, and fits calibration table on the best.
        """
        results: List[VariantResult] = [self.evaluate_variant(v) for v in variants]
        # Rank by F1 descending, then precision
        results.sort(key=lambda r: (r.f1_score, r.precision), reverse=True)

        summary_rows = []
        for r in results:
            summary_rows.append({
                "variant": r.variant.name,
                "pocket_tool": r.variant.pocket_tool,
                "min_plddt": r.variant.min_plddt,
                "f1_score": r.f1_score,
                "precision": r.precision,
                "recall": r.recall,
                "accuracy": r.accuracy,
                "avg_latency_ms": r.avg_latency_ms
            })

        summary_df = pd.DataFrame(summary_rows)
        best = results[0]

        # Fit calibration table using best variant predictions
        cal_table = CalibrationTable(best.predictions_df)

        return best, summary_df, cal_table
