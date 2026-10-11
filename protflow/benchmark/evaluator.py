"""
Automatic Node Evaluator for ProtFlow.
When a new tool manifest is added, compares it against the baseline node
on the standard benchmark dataset and produces an automated quality comparison report.
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict, List, Optional
import time
import pandas as pd
import numpy as np

from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.types import BaseBioType, Structure, Pocket


@dataclass
class ToolBenchmarkReport:
    candidate_name: str
    baseline_name: str
    samples_tested: int
    candidate_avg_ms: float
    baseline_avg_ms: float
    speedup: float
    mean_candidate_score: float
    mean_baseline_score: float
    recommendation: str
    comparison_df: pd.DataFrame


class NodeEvaluator:
    """Automates benchmarking when a new tool manifest is added to the system."""
    def __init__(self, runner: Optional[Runner] = None):
        self.runner = runner or Runner(enable_cache=False)

    def compare_nodes(
        self,
        candidate_manifest: Manifest,
        baseline_manifest: Manifest,
        test_inputs_list: List[Dict[str, BaseBioType]]
    ) -> ToolBenchmarkReport:
        rows = []
        cand_times = []
        base_times = []

        for item in test_inputs_list:
            first_in = next(iter(item.values()))
            pid = first_in.protein_id

            # Run baseline
            t0 = time.perf_counter()
            out_base = self.runner.run(baseline_manifest, item)
            t_base_ms = (time.perf_counter() - t0) * 1000.0
            base_times.append(t_base_ms)

            # Run candidate
            t0 = time.perf_counter()
            out_cand = self.runner.run(candidate_manifest, item)
            t_cand_ms = (time.perf_counter() - t0) * 1000.0
            cand_times.append(t_cand_ms)

            # Extract primary score from first output object
            base_obj = next(iter(out_base.values()))
            cand_obj = next(iter(out_cand.values()))

            base_val = float(base_obj.df["score"].iloc[0]) if hasattr(base_obj, "df") and "score" in base_obj.df.columns else 0.0
            cand_val = float(cand_obj.df["score"].iloc[0]) if hasattr(cand_obj, "df") and "score" in cand_obj.df.columns else 0.0

            rows.append({
                "protein_id": pid,
                f"{baseline_manifest.name}_score": round(base_val, 3),
                f"{candidate_manifest.name}_score": round(cand_val, 3),
                "score_diff": round(cand_val - base_val, 3),
                f"{baseline_manifest.name}_ms": round(t_base_ms, 1),
                f"{candidate_manifest.name}_ms": round(t_cand_ms, 1),
            })

        comp_df = pd.DataFrame(rows)
        avg_base_ms = float(np.mean(base_times)) if base_times else 1.0
        avg_cand_ms = float(np.mean(cand_times)) if cand_times else 1.0
        speedup = round(avg_base_ms / max(avg_cand_ms, 0.001), 4)

        mean_base = float(comp_df[f"{baseline_manifest.name}_score"].mean())
        mean_cand = float(comp_df[f"{candidate_manifest.name}_score"].mean())

        # Formulate plain English recommendation
        if mean_cand >= mean_base and speedup >= 1.0:
            rec = f"RECOMMENDED: '{candidate_manifest.name}' matches or exceeds baseline score (+{round(mean_cand-mean_base, 3)}) with {speedup}x speed."
        elif mean_cand >= mean_base:
            rec = f"ACCEPTABLE: '{candidate_manifest.name}' improves score (+{round(mean_cand-mean_base, 3)}) but is slightly slower ({speedup}x)."
        else:
            rec = f"INFERIOR: Baseline '{baseline_manifest.name}' demonstrates higher confidence on benchmark proteins."

        return ToolBenchmarkReport(
            candidate_name=candidate_manifest.name,
            baseline_name=baseline_manifest.name,
            samples_tested=len(test_inputs_list),
            candidate_avg_ms=round(avg_cand_ms, 1),
            baseline_avg_ms=round(avg_base_ms, 1),
            speedup=speedup,
            mean_candidate_score=round(mean_cand, 3),
            mean_baseline_score=round(mean_base, 3),
            recommendation=rec,
            comparison_df=comp_df
        )
