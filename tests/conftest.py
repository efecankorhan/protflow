"""
Pytest fixtures and test helpers for ProtFlow.
"""

import os
from pathlib import Path
import gemmi
import pytest
from protflow.types import Sequence, Structure

# Tests run offline: opt in to synthetic structures (off by default in production).
os.environ.setdefault("PROTFLOW_ALLOW_SYNTHETIC", "1")


@pytest.fixture
def fixtures_dir(tmp_path) -> Path:
    fdir = tmp_path / "fixtures"
    fdir.mkdir(parents=True, exist_ok=True)
    return fdir


@pytest.fixture
def sample_data(fixtures_dir: Path):
    # 1. FASTA
    fasta_path = fixtures_dir / "test_prot.fasta"
    seq_content = ">P01234 Test Protein\nMAG\n"
    fasta_path.write_text(seq_content, encoding="utf-8")
    seq_obj = Sequence(protein_id="P01234", sequence="MAG", description="Test Protein")

    # 2. CIF
    cif_path = fixtures_dir / "test_prot.cif"
    st = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    for idx, (resname, el, pos) in enumerate([
        ("MET", "C", (0.0, 0.0, 0.0)),
        ("ALA", "C", (1.5, 0.0, 0.0)),
        ("GLY", "C", (3.0, 0.0, 0.0)),
    ], start=1):
        res = gemmi.Residue()
        res.name = resname
        res.seqid = gemmi.SeqId(str(idx))
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element(el)
        atom.pos = gemmi.Position(*pos)
        res.add_atom(atom)
        chain.add_residue(res)

    model.add_chain(chain)
    st.add_model(model)
    doc = st.make_mmcif_document()
    doc.write_file(str(cif_path))

    struct_obj = Structure(protein_id="P01234", file_path=str(cif_path), format="cif")

    return {
        "fasta_path": fasta_path,
        "sequence": seq_obj,
        "cif_path": cif_path,
        "structure": struct_obj,
    }
