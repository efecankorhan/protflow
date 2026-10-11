"""MMseqs2 homology node: real search against a small target FASTA."""

import pytest

from protflow.adapters.homology import find_mmseqs, parse_mmseqs_homology
from protflow.types import Sequence, SequenceCollection

try:
    find_mmseqs()
    HAVE_MMSEQS = True
except FileNotFoundError:
    HAVE_MMSEQS = False

KINASE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"
UNRELATED = "WWHHCCPPYYMMWWHHCCPPYYMMWWHHCCPPYYMMWWHHCCPPYYMM"


@pytest.mark.skipif(not HAVE_MMSEQS, reason="mmseqs not installed")
def test_mmseqs_homology_hits_and_no_hits(tmp_path):
    target = tmp_path / "t.fasta"
    target.write_text(f">sp|T1|KINASE_TEST Test kinase\n{KINASE}\n")
    coll = SequenceCollection(protein_id="G", sequences={
        "hit": Sequence(protein_id="hit", sequence=KINASE),
        "orphan": Sequence(protein_id="orphan", sequence=UNRELATED),
    })
    out = parse_mmseqs_homology(None, {"inputs": {"proteins": coll}, "parameters": {"target_fasta": str(target)}})
    df = out["homology"].df
    assert list(df["query"]) == ["hit"] and df["pident"].iloc[0] == 100.0
    p = out["homology"].provenance.parameters
    assert (p["input_count"], p["queries_with_hits"], p["no_hit"]) == (2, 1, 1)
