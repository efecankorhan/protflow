"""Prioritization combiner."""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import json
import re
import numpy as np
import pandas as pd
import requests
import gemmi

from protflow.types import (
    Sequence,
    Structure,
    ResidueScore,
    Pocket,
    Alignment,
    Annotation,
    Table,
    Provenance
)
from protflow.adapters.residue import AA3_TO_1, AA1_TO_3


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
    pid = inputs["quality_table"].protein_id

    mean_plddt = float(inputs["quality_table"].df["mean_plddt"].iloc[0])
    pocket_score = float(inputs["pockets"].df["score"].iloc[0])
    pocket_vol = float(inputs["pockets"].df["volume"].iloc[0])
    sasa_hydro = float(inputs["sasa_table"].df["hydrophobic_ratio"].iloc[0])
    novelty = float(inputs["foldseek_table"].df["novelty_score"].iloc[0])

    # Expression & solubility estimation: balanced hydrophobicity (~0.35 - 0.45 ideal)
    solubility_est = max(0.0, 1.0 - abs(sasa_hydro - 0.40) * 2.0)
    quality_norm = mean_plddt / 100.0

    # Check if upstream filter passed
    q_df = inputs["quality_table"].df
    passed_plddt = bool(q_df["passed_plddt_filter"].iloc[0]) if "passed_plddt_filter" in q_df.columns else True

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
