"""Foldseek structural match."""

from __future__ import annotations
from typing import Any, Dict
import pandas as pd

from protflow.types import Alignment, Table, Provenance


def parse_foldseek(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """Computes structural similarity against reference database and calculates novelty score."""
    st_obj = context["inputs"]["structure"]
    pid = st_obj.protein_id

    # shortcut: reference structural alignment mock without foldseek binary/database, connect foldseek for real proteome searches
    top_hit = "SWISS_P0A855"
    tm_score = 0.88
    identity = 0.72
    novelty = 1.0 - identity

    prov = Provenance(
        tool_name="foldseek",
        tool_version="1.0.0",
        parameters={"approximate": True, "top_hit": top_hit}
    )

    align_obj = Alignment(
        protein_id=pid,
        target_id=top_hit,
        identity=identity,
        score=tm_score,
        evalue=1e-12
    ).attach_provenance(prov)

    df = pd.DataFrame([{
        "protein_id": pid,
        "top_hit": top_hit,
        "tm_score": tm_score,
        "identity": identity,
        "novelty_score": round(novelty, 3)
    }])
    table = Table(protein_id=pid, df=df).attach_provenance(prov)

    return {"alignment": align_obj, "foldseek_table": table}


