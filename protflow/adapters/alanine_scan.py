"""
Wet Lab Alanine Scanning & Experimental Feedback Loop for ProtFlow.
Generates alanine mutation panels for top candidate enzyme pockets
and integrates experimental wet lab validation back into benchmark calibration.
"""

from __future__ import annotations
from typing import Any, Dict, List, Tuple
import pandas as pd
from protflow.benchmark.calibration import CalibrationTable


def generate_alanine_mutants(
    priority_df: pd.DataFrame,
    pockets_df: pd.DataFrame,
    top_n: int = 10,
    max_mutants_per_protein: int = 3
) -> pd.DataFrame:
    """
    Selects top candidates from the priority table and designs alanine scanning mutants
    targeting key residues in the predicted catalytic pocket.
    """
    if priority_df.empty or "final_priority_score" not in priority_df.columns or pockets_df.empty:
        return pd.DataFrame()

    candidates = priority_df.sort_values(by="final_priority_score", ascending=False).head(top_n)
    rows = []

    for _, cand in candidates.iterrows():
        pid = cand["protein_id"]
        score = cand["final_priority_score"]

        # Find matching pocket residues
        cand_pockets = pockets_df[pockets_df["protein_id"] == pid]
        if cand_pockets.empty or "residues" not in cand_pockets.columns:
            continue
        res_val = cand_pockets["residues"].iloc[0]
        if not res_val:
            continue
        res_str = str(res_val)

        residue_items = [r.strip() for r in res_str.split(",") if r.strip()][:max_mutants_per_protein]

        for item in residue_items:
            # Format: 'A:24' or 'A:24_GLU'
            chain = item.split(":")[0] if ":" in item else "A"
            pos_part = item.split(":")[1] if ":" in item else item
            pos = pos_part.split("_")[0] if "_" in pos_part else pos_part
            wt_res = pos_part.split("_")[1] if "_" in pos_part else "X"

            mutation_code = f"{chain}:{pos} {wt_res}->ALA"
            rows.append({
                "protein_id": pid,
                "candidate_priority_score": score,
                "chain_id": chain,
                "position": pos,
                "wildtype": wt_res,
                "mutant": "ALA",
                "mutation_code": mutation_code,
                "experimental_assay": "Alanine-Scanning Spectrophotometric Activity Assay",
                "validation_rationale": "Disrupts binding pocket cavity to test catalytic necessity"
            })

    return pd.DataFrame(rows)


def record_wet_lab_feedback(
    experimental_records: List[Dict[str, Any]],
    current_benchmark_df: pd.DataFrame,
    calibration_table: CalibrationTable
) -> Tuple[pd.DataFrame, CalibrationTable]:
    """
    Appends wet lab alanine scan assay results into the benchmark ledger
    and recalibrates the empirical confidence table.
    """
    new_rows = []
    for exp in experimental_records:
        new_rows.append({
            "protein_id": exp["protein_id"],
            "uniprot_acc": exp.get("uniprot_acc", "WET_LAB_EXP"),
            "name": exp.get("name", exp["protein_id"]),
            "ec_class": exp.get("ec_class", "EC 3 (Hydrolase)"),
            "is_active_enzyme": bool(exp.get("activity_validated", True)),
            "sequence": exp.get("sequence", "MAGLVAILEKFGDELRWAP"),
            "known_pocket_residues": exp.get("pocket_residues", ""),
            "description": f"Wet lab validated candidate: {exp.get('notes', '')}",
            "final_priority_score": float(exp.get("measured_score", 0.85))
        })

    updated_df = pd.concat([current_benchmark_df, pd.DataFrame(new_rows)], ignore_index=True)
    # Refit calibration table with newly gained empirical lab ground-truth
    calibration_table.fit(updated_df, score_col="final_priority_score", label_col="is_active_enzyme")
    return updated_df, calibration_table
