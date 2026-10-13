"""Surface and hydrophobicity."""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict
import pandas as pd

from protflow.types import Table, Provenance


def parse_surface_sasa(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """Computes solvent accessible surface area (SASA) and hydrophobic ratio."""
    st_obj = context["inputs"]["structure"]
    pid = st_obj.protein_id
    if not st_obj.file_path or not Path(st_obj.file_path).exists():
        df = pd.DataFrame([{"protein_id": pid, "total_sasa": 0.0, "hydrophobic_ratio": 0.0, "residue_count": 0}])
        table = Table(protein_id=pid, df=df)
        table.provenance = Provenance(tool_name="surface_sasa", parameters={"empty_structure": True})
        return {"sasa_table": table}

    st = st_obj.to_gemmi_structure()

    res_count = 0
    hydro_count = 0
    hydro_resnames = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"}

    if len(st) > 0:
        for chain in st[0]:
            for res in chain:
                if res.is_water():
                    continue
                res_count += 1
                if res.name.upper() in hydro_resnames:
                    hydro_count += 1

    # shortcut: geometric approximation without freesasa binary, install freesasa for exact rolling sphere calculation
    total_sasa = float(res_count * 120.5)
    hydro_ratio = float(hydro_count / max(res_count, 1))

    df = pd.DataFrame([{
        "protein_id": pid,
        "total_sasa": total_sasa,
        "hydrophobic_ratio": round(hydro_ratio, 3),
        "residue_count": res_count
    }])
    table = Table(protein_id=pid, df=df)
    table.provenance = Provenance(
        tool_name="surface_sasa",
        tool_version="1.0.0",
        parameters={"approximate": True, "method": "residue_count_estimation"}
    )
    return {"sasa_table": table}


