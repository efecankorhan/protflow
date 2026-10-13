"""Prioritization combiner."""

from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from protflow.types import Table


# --- 8. Prioritization Combiner ---

def parse_prioritization(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Combines all upstream metrics into a single ranked priority score:
    - Structure quality (mean pLDDT)
    - Pocket confidence & volume
    - Heterologous expression / solubility
    - Novelty component
    """
    inputs = context["inputs"]
    q_table = inputs.get("quality_table")
    pid = q_table.protein_id if q_table else "SAMPLE"

    q_df = q_table.df if q_table and hasattr(q_table, "df") else pd.DataFrame()
    mean_plddt = float(q_df["mean_plddt"].iloc[0]) if not q_df.empty and "mean_plddt" in q_df.columns else 0.0

    pockets_obj = inputs.get("pockets")
    p_df = pockets_obj.df if pockets_obj and hasattr(pockets_obj, "df") else pd.DataFrame()
    pocket_score = float(p_df["score"].iloc[0]) if not p_df.empty and "score" in p_df.columns else 0.0
    pocket_vol = float(p_df["volume"].iloc[0]) if not p_df.empty and "volume" in p_df.columns else 0.0

    sasa_obj = inputs.get("sasa_table")
    s_df = sasa_obj.df if sasa_obj and hasattr(sasa_obj, "df") else pd.DataFrame()
    sasa_hydro = float(s_df["hydrophobic_ratio"].iloc[0]) if not s_df.empty and "hydrophobic_ratio" in s_df.columns else 0.40

    foldseek_obj = inputs.get("foldseek_table")
    f_df = foldseek_obj.df if foldseek_obj and hasattr(foldseek_obj, "df") else pd.DataFrame()
    novelty = float(f_df["novelty_score"].iloc[0]) if not f_df.empty and "novelty_score" in f_df.columns else 0.0

    # Expression & solubility estimation: balanced hydrophobicity (~0.35 - 0.45 ideal)
    solubility_est = max(0.0, 1.0 - abs(sasa_hydro - 0.40) * 2.0)
    quality_norm = mean_plddt / 100.0

    # Check if upstream filter passed
    passed_plddt = bool(q_df["passed_plddt_filter"].iloc[0]) if not q_df.empty and "passed_plddt_filter" in q_df.columns else True

    if not passed_plddt:
        final_score = 0.0
        status = "DROPPED_FILTER"
    else:
        # Composite priority score (0.0 to 1.0)
        final_score = (
            0.30 * quality_norm +
            0.30 * pocket_score +
            0.20 * solubility_est +
            0.20 * novelty
        )
        status = "CANDIDATE" if final_score >= 0.70 else "DEPURGED"

    summary_df = pd.DataFrame([{
        "protein_id": pid,
        "final_priority_score": round(final_score, 4),
        "structure_quality": round(quality_norm, 3),
        "pocket_score": round(pocket_score, 3),
        "pocket_volume": pocket_vol,
        "solubility_est": round(solubility_est, 3),
        "novelty_score": round(novelty, 3),
        "status": status
    }])

    return {"priority_table": Table(protein_id=pid, df=summary_df)}
