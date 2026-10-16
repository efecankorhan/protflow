"""Unit tests for MAFFT / Multiple Sequence Alignment adapter."""

import pytest
import pandas as pd
from protflow.types import Sequence, SequenceCollection, Alignment, Table
from protflow.adapters.msa import (
    _blosum62_score,
    _needleman_wunsch,
    _run_progressive_msa,
    compute_conservation,
    parse_mafft
)


def test_blosum62_and_needleman_wunsch():
    # Identity match should be positive
    assert _blosum62_score('W', 'W') == 11
    assert _blosum62_score('A', 'A') == 4
    # Mismatch
    assert _blosum62_score('W', 'D') == -4

    # Pairwise alignment
    s1 = "HEAGAWGHEE"
    s2 = "PRAWHEAE"
    a1, a2, score = _needleman_wunsch(s1, s2)
    assert len(a1) == len(a2)
    assert a1.replace('-', '') == s1
    assert a2.replace('-', '') == s2
    assert score > -50


def test_progressive_msa_and_conservation():
    seqs = {
        "protA": "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSH",
        "protB": "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKTH",
        "protC": "MVLSEGEWQLVLNVWAKVEADVAGHGQDILIRLFKSH",
    }
    aligned = _run_progressive_msa(seqs)
    assert len(aligned) == 3
    lens = [len(v) for v in aligned.values()]
    assert len(set(lens)) == 1  # All aligned sequences must have equal length

    cons_df = compute_conservation(aligned)
    assert not cons_df.empty
    assert "conservation_score" in cons_df.columns
    assert "entropy" in cons_df.columns
    assert "consensus_aa" in cons_df.columns

    # First position 'M' is 100% conserved
    pos1 = cons_df.iloc[0]
    assert pos1["consensus_aa"] == 'M'
    assert pos1["conservation_score"] == 1.0
    assert pos1["entropy"] == 0.0


def test_parse_mafft_end_to_end():
    collection = SequenceCollection(protein_id="test_family", sequences={
        "seq1": Sequence(protein_id="seq1", sequence="ACDEFGHIKLMNPQRSTVWY"),
        "seq2": Sequence(protein_id="seq2", sequence="ACDEFGHIKLMNPQRSTVWY"),
        "seq3": Sequence(protein_id="seq3", sequence="ACDEFGHIKAMNPQRSTVWY"),
    })

    out = parse_mafft(None, {"inputs": {"sequences": collection}})

    assert "alignment" in out
    assert "aligned_sequences" in out
    assert "conservation" in out

    assert isinstance(out["alignment"], Alignment)
    assert isinstance(out["aligned_sequences"], SequenceCollection)
    assert isinstance(out["conservation"], Table)

    assert len(out["aligned_sequences"].sequences) == 3
    assert out["alignment"].identity > 0.9
    assert len(out["conservation"].df) == 20


def test_parse_famsa_alias_and_manifest():
    from protflow.adapters.msa import parse_famsa
    from protflow.manifest import Manifest

    m = Manifest.from_yaml("manifests/12_msa_famsa.yaml")
    assert m.name == "msa_famsa"
    assert m.inputs["sequences"] == "SequenceCollection"

    collection = SequenceCollection(protein_id="test_family", sequences={
        "seq1": Sequence(protein_id="seq1", sequence="ACDEFGHIKLMNPQRSTVWY"),
        "seq2": Sequence(protein_id="seq2", sequence="ACDEFGHIKAMNPQRSTVWY"),
    })
    out = parse_famsa(None, {"inputs": {"sequences": collection}})
    assert "alignment" in out
    assert "conservation" in out


def test_parse_mafft_single_sequence():
    from protflow.adapters.msa import parse_mafft
    collection = SequenceCollection(protein_id="single", sequences={
        "seq1": Sequence(protein_id="seq1", sequence="ACDEFGHIKLMNPQRSTVWY"),
    })
    out = parse_mafft(None, {"inputs": {"sequences": collection}})
    assert out["aligned_sequences"].sequences["seq1"].sequence == "ACDEFGHIKLMNPQRSTVWY"
    assert out["alignment"].identity == 1.0
    assert out["aligned_sequences"].provenance.tool_name == "passthrough"


def test_progressive_msa_multi_sequence_insertions():
    seqs = {
        "s1": "MKTLLVV",
        "s2": "MKAAATLLVV",
        "s3": "MKTLLGGGVV",
        "s4": "MKTLLVV",
    }
    aligned = _run_progressive_msa(seqs)
    assert len(set(len(v) for v in aligned.values())) == 1
    for k, orig in seqs.items():
        assert aligned[k].replace("-", "") == orig
    assert aligned["s1"] == aligned["s4"]


