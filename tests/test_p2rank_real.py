"""Tests for real P2Rank pocket prediction and PyMOL PML export."""

from pathlib import Path
import pytest
from protflow.adapters.pockets import parse_p2rank, _find_p2rank_and_jre
from protflow.types import Structure, Pocket

p2_dir, jre_dir = _find_p2rank_and_jre()
HAVE_P2RANK = (p2_dir is not None and jre_dir is not None)


@pytest.mark.skipif(not HAVE_P2RANK, reason="P2Rank or JRE not found")
def test_p2rank_real_prediction_and_pymol_export(tmp_path):
    root = Path(__file__).resolve().parents[1]
    pdb_file = root / "tools" / "p2rank_2.5.1" / "test_data" / "1fbl.pdb"
    st = Structure(protein_id="1FBL", file_path=str(pdb_file), format="pdb")

    out = parse_p2rank(None, {
        "inputs": {"structure": st},
        "output_dir": str(tmp_path)
    })

    assert "pockets" in out
    pockets = out["pockets"]
    assert isinstance(pockets, Pocket)
    assert len(pockets.df) >= 3
    assert "center_x" in pockets.df.columns
    assert "residues" in pockets.df.columns
    assert pockets.df["score"].iloc[0] > 5.0

    # Check PyMOL .pml visualization script in provenance
    prov = getattr(pockets, "provenance", None)
    assert prov is not None
    pml = prov.parameters.get("pymol_pml", "")
    assert "from pymol import cmd" in pml
    assert "load" in pml
