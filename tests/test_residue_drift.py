"""
Unit tests specifically targeting residue drift, chain ID mismatch,
and numbering offsets between tool outputs and 3D structures.
"""

import pandas as pd
import pytest
from protflow.adapters.residue import (
    validate_residue_table,
    extract_structure_residues,
    ResidueDriftError
)
from protflow.types import Sequence


def test_valid_residue_table(sample_data):
    st = sample_data["structure"]
    seq = sample_data["sequence"]

    df = pd.DataFrame([
        {"chain_id": "A", "auth_seq_id": "1", "resname": "MET", "plddt": 95.0},
        {"chain_id": "A", "auth_seq_id": "2", "resname": "ALA", "plddt": 91.0},
        {"chain_id": "A", "auth_seq_id": "3", "resname": "GLY", "plddt": 88.0},
    ])

    validated = validate_residue_table(df, st, sequence=seq)
    assert len(validated) == 3
    assert list(validated["one_letter"]) == ["M", "A", "G"]


def test_detects_shifted_numbering_drift(sample_data):
    st = sample_data["structure"]
    seq = sample_data["sequence"]

    # Silent index drift: Tool output starts from 10 instead of 1
    drifted_df = pd.DataFrame([
        {"chain_id": "A", "auth_seq_id": "10", "resname": "MET"},
        {"chain_id": "A", "auth_seq_id": "11", "resname": "ALA"},
    ])

    with pytest.raises(ResidueDriftError, match="does not exist in structure"):
        validate_residue_table(drifted_df, st, sequence=seq)


def test_detects_residue_mutation_mismatch(sample_data):
    st = sample_data["structure"]
    seq = sample_data["sequence"]

    # Position 2 is ALA in structure, but tool erroneously says VAL (or 'V')
    mismatched_df = pd.DataFrame([
        {"chain_id": "A", "auth_seq_id": "1", "resname": "MET"},
        {"chain_id": "A", "auth_seq_id": "2", "resname": "VAL"},
    ])

    with pytest.raises(ResidueDriftError, match="Residue mismatch at A:2"):
        validate_residue_table(mismatched_df, st, sequence=seq)


def test_detects_chain_mismatch(sample_data):
    st = sample_data["structure"]

    # Structure has Chain A, tool says Chain B
    wrong_chain_df = pd.DataFrame([
        {"chain_id": "B", "auth_seq_id": "1", "resname": "MET"},
    ])

    with pytest.raises(ResidueDriftError, match="does not exist in structure"):
        validate_residue_table(wrong_chain_df, st)
