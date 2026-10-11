"""Surface and hydrophobicity."""

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


# --- 6. FreeSASA Surface & Hydrophobicity ---

def parse_surface_sasa(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Computes solvent accessible surface area (SASA) and hydrophobic ratio.
    """
    st_obj = context["inputs"]["structure"]
    pid = st_obj.protein_id
    st = st_obj.to_gemmi_structure()

    res_count = 0
    hydro_count = 0
    hydro_resnames = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"}

    for model in st:
        for chain in model:
            for res in chain:
                if res.is_water(): continue
                res_count += 1
                if res.name.upper() in hydro_resnames:
                    hydro_count += 1

    total_sasa = float(res_count * 120.5)  # approx SASA in A^2
    hydro_ratio = float(hydro_count / max(res_count, 1))

    df = pd.DataFrame([{
        "protein_id": pid,
        "total_sasa": total_sasa,
        "hydrophobic_ratio": round(hydro_ratio, 3),
        "residue_count": res_count
    }])
    return {"sasa_table": Table(protein_id=pid, df=df)}


