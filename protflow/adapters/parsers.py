"""
Standard adapter parsers for converting raw tool outputs to ProtFlow types.
"""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, Union
import json
import pandas as pd
from protflow.types import Pocket, ResidueScore, Structure, Sequence, Table


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


def parse_alphafold_db_http_response(raw_data: Union[str, Dict[str, Any]], context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Adapter that parses AlphaFold DB API response to extract structure URL and pLDDT residue scores.
    """
    protein_id = context.get("protein_id", "unknown")
    if isinstance(raw_data, str):
        data = json.loads(raw_data)
    else:
        data = raw_data

    # AFDB response is a list of entries
    entry = data[0] if isinstance(data, list) and data else data
    cif_url = entry.get("cifUrl", "")
    mean_plddt = float(entry.get("globalMetricValue", 0.0))

    # Return structure record and metadata table
    meta_df = pd.DataFrame([{
        "protein_id": protein_id,
        "cif_url": cif_url,
        "global_plddt": mean_plddt,
        "gene": entry.get("gene", ""),
        "organism": entry.get("organismScientificName", "")
    }])

    return {
        "structure_meta": Table(protein_id=protein_id, df=meta_df)
    }
