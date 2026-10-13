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
        if hasattr(struct_input, "structures"):
            from protflow.types import StructureCollection
            out_structures = StructureCollection(
                protein_id=struct_input.protein_id,
                structures={k: v for k, v in struct_input.structures.items() if passed}
            )
            result["structure"] = out_structures
        else:
            result["structure"] = struct_input if passed else Structure(protein_id=pid, file_path="")
    return result



