"""pLDDT quality filter."""

from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from protflow.types import Structure, Table, Provenance


# --- 4. pLDDT Quality Filter ---

def parse_plddt_filter(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Filters structures with mean pLDDT below threshold (default: 70.0).
    Passes through filtered Structure / StructureCollection and outputs quality_table.
    """
    scores_obj = context["inputs"]["residue_scores"]
    struct_input = context["inputs"].get("structure")
    pid = scores_obj.protein_id
    min_thresh = float(context["parameters"].get("min_plddt", 70.0))
    df = scores_obj.df

    if struct_input is not None and hasattr(struct_input, "structures"):
        from protflow.types import StructureCollection
        passed_structures = {}
        rows = []
        for s_pid, s_obj in struct_input.structures.items():
            s_df = df[df["protein_id"] == s_pid] if "protein_id" in df.columns else df
            s_mean = float(s_df["score"].mean()) if not s_df.empty else 0.0
            s_passed = s_mean >= min_thresh
            if s_passed:
                passed_structures[s_pid] = s_obj
            rows.append({
                "protein_id": s_pid,
                "mean_plddt": s_mean,
                "min_threshold": min_thresh,
                "passed_plddt_filter": s_passed,
                "status": "PASSED" if s_passed else "DROPPED",
                "drop_reason": "" if s_passed else f"Low pLDDT ({s_mean:.1f} < {min_thresh})"
            })
        res_df = pd.DataFrame(rows)
        table_obj = Table(protein_id=struct_input.protein_id, df=res_df)
        dropped_ct = len(struct_input.structures) - len(passed_structures)
        table_obj.provenance = Provenance(
            tool_name="plddt_filter",
            parameters={
                "input_count": len(struct_input.structures),
                "passed_count": len(passed_structures),
                "dropped_count": dropped_ct,
                "reasons": {f"pLDDT < {min_thresh}": dropped_ct} if dropped_ct else {}
            }
        )
        return {
            "quality_table": table_obj,
            "structure": StructureCollection(protein_id=struct_input.protein_id, structures=passed_structures)
        }

    mean_score = float(df["score"].mean()) if not df.empty else 0.0
    passed = mean_score >= min_thresh

    res_df = pd.DataFrame([{
        "protein_id": pid,
        "mean_plddt": mean_score,
        "min_threshold": min_thresh,
        "passed_plddt_filter": passed,
        "status": "PASSED" if passed else "DROPPED",
        "drop_reason": "" if passed else f"Low pLDDT ({mean_score:.1f} < {min_thresh})"
    }])
    table_obj = Table(protein_id=pid, df=res_df)
    table_obj.provenance = Provenance(
        tool_name="plddt_filter",
        parameters={
            "input_count": 1,
            "passed_count": 1 if passed else 0,
            "dropped_count": 0 if passed else 1,
            "reasons": {f"pLDDT < {min_thresh}": 1} if not passed else {}
        }
    )

    result = {"quality_table": table_obj}
    if struct_input is not None:
        result["structure"] = struct_input if passed else Structure(protein_id=pid, file_path="")
    return result



