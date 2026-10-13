"""Subcellular localization filter."""

from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from protflow.types import Sequence, Table, Provenance


# --- 2. Subcellular Localization / Filtering (DeepLoc / SignalP) ---

def parse_subcellular_filter(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Adapter for Subcellular Localization filter.
    Scores N-terminal hydrophobic segment to detect signal peptides / transmembrane domains.
    Supports single Sequence and SequenceCollection, passing through intracellular proteins.
    """
    seq_input = context["inputs"].get("proteins") or context["inputs"]["protein"]

    if isinstance(seq_input, Sequence):
        items = [seq_input]
        is_collection = False
    elif hasattr(seq_input, "sequences"):
        items = list(seq_input.sequences.values())
        is_collection = True
    else:
        items = [seq_input]
        is_collection = False

    rows = []
    passed_items = []
    dropped_count = 0
    hydrophobic = set("LIVFWMAG")

    for seq_obj in items:
        seq = seq_obj.sequence.upper()
        pid = seq_obj.protein_id
        n_term = seq[:25]
        n_hydro_count = sum(1 for aa in n_term if aa in hydrophobic)
        hydro_fraction = n_hydro_count / max(len(n_term), 1)

        is_secreted = hydro_fraction > 0.60
        loc = "Extracellular" if is_secreted else "Cytoplasm"
        is_intracellular = not is_secreted
        conf = 0.92 if is_intracellular else 0.85
        status = "PASSED" if is_intracellular else "DROPPED"
        reason = "" if is_intracellular else "Signal peptide (Extracellular)"

        if is_intracellular:
            passed_items.append(seq_obj)
        else:
            dropped_count += 1

        rows.append({
            "protein_id": pid,
            "localization": loc,
            "is_intracellular": is_intracellular,
            "confidence": conf,
            "status": status,
            "drop_reason": reason
        })

    df = pd.DataFrame(rows)
    table_obj = Table(protein_id=items[0].protein_id if items else "NONE", df=df)
    table_obj.provenance = Provenance(
        tool_name="subcellular_filter",
        parameters={
            "input_count": len(items),
            "passed_count": len(passed_items),
            "dropped_count": dropped_count,
            "reasons": {"Extracellular": dropped_count} if dropped_count else {}
        }
    )

    if is_collection:
        from protflow.types import SequenceCollection
        out_collection = SequenceCollection(
            protein_id=seq_input.protein_id,
            sequences={s.protein_id: s for s in passed_items}
        )
        return {"proteins": out_collection, "filter_table": table_obj}
    else:
        out_prot = passed_items[0] if passed_items else Sequence(protein_id=seq_input.protein_id, sequence="", description="DROPPED_BY_FILTER")
        return {"protein": out_prot, "filter_table": table_obj}



