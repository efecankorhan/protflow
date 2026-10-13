"""
Standard adapter parsers for converting raw tool outputs to ProtFlow types.
"""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Union
import json
import pandas as pd
from protflow.types import Pocket


def parse_dummy_pocket_finder(raw_data: Union[str, Dict[str, Any]], context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Adapter that parses dummy/reference pocket finder output into a Pocket object.
    Raw data can be a json string, dict, or path to output file.
    """
    protein_id = context.get("protein_id", "unknown")
    if isinstance(raw_data, (str, Path)) and Path(raw_data).exists():
        raw_text = Path(raw_data).read_text(encoding="utf-8")
        data = json.loads(raw_text)
    elif isinstance(raw_data, str):
        data = json.loads(raw_data)
    else:
        data = raw_data

    # Expected format: list of pockets or dict with "pockets"
    pockets = data.get("pockets", []) if isinstance(data, dict) else data
    rows = []
    for p in pockets:
        rows.append({
            "protein_id": protein_id,
            "pocket_id": p.get("pocket_id", 1),
            "score": float(p.get("score", 0.0)),
            "center_x": float(p.get("center_x", 0.0)),
            "center_y": float(p.get("center_y", 0.0)),
            "center_z": float(p.get("center_z", 0.0)),
            "volume": float(p.get("volume", 0.0)),
            "residues": p.get("residues", "")
        })

    df = pd.DataFrame(rows)
    return {"pockets": Pocket(protein_id=protein_id, df=df)}
