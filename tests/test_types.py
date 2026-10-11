"""
Tests for ProtFlow 7 core data types and provenance metadata.
"""

from pathlib import Path
import pandas as pd
import pytest
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


def test_sequence_fasta_io(tmp_path):
    seq = Sequence(protein_id="PROT_1", sequence="MKFLV", description="kinase domain")
    prov = Provenance(tool_name="prodigal", tool_version="2.6.3")
    seq.attach_provenance(prov)

    out_file = tmp_path / "out.fasta"
    seq.save_fasta(out_file)

    loaded = Sequence.from_fasta(out_file)
    assert loaded.protein_id == "PROT_1"
    assert loaded.sequence == "MKFLV"
    assert "kinase domain" in loaded.description


def test_structure_cif(sample_data):
    st = sample_data["structure"]
    chains = st.get_chains()
    assert chains == ["A"]


def test_residue_score_parquet_io(tmp_path):
    df = pd.DataFrame([
        {"protein_id": "P1", "chain_id": "A", "auth_seq_id": "1", "resname": "MET", "score_type": "plddt", "score": 92.5},
        {"protein_id": "P1", "chain_id": "A", "auth_seq_id": "2", "resname": "ALA", "score_type": "plddt", "score": 88.0},
    ])
    rs = ResidueScore(protein_id="P1", df=df)
    out_file = tmp_path / "scores.parquet"
    rs.save_parquet(out_file)

    loaded = ResidueScore.from_parquet(out_file)
    assert loaded.protein_id == "P1"
    assert len(loaded.df) == 2
    assert loaded.df["score"].iloc[0] == 92.5


def test_pocket_parquet_io(tmp_path):
    df = pd.DataFrame([
        {
            "protein_id": "P1",
            "pocket_id": 1,
            "score": 0.85,
            "center_x": 10.0,
            "center_y": 15.2,
            "center_z": -3.4,
            "volume": 450.0,
            "residues": "A:1,A:2"
        }
    ])
    pkt = Pocket(protein_id="P1", df=df)
    out_file = tmp_path / "pockets.parquet"
    pkt.save_parquet(out_file)

    loaded = Pocket.from_parquet(out_file)
    assert loaded.protein_id == "P1"
    assert loaded.df["score"].iloc[0] == 0.85
    assert loaded.df["volume"].iloc[0] == 450.0


def test_alignment_record():
    aln = Alignment(
        protein_id="P1",
        target_id="SWISS_P999",
        query_seq="MKFLV",
        target_seq="MK-LV",
        identity=0.8,
        score=150.0,
        evalue=1e-5
    )
    d = aln.to_dict()
    assert d["protein_id"] == "P1"
    assert d["identity"] == 0.8


def test_annotation_gff3_io(tmp_path):
    records = [{
        "seqid": "P1",
        "source": "SignalP",
        "type": "signal_peptide",
        "start": 1,
        "end": 20,
        "score": "0.98",
        "strand": "+",
        "phase": ".",
        "attributes": {"ID": "SP_1", "Name": "Sec/SPI"}
    }]
    annot = Annotation(protein_id="P1", records=records)
    out_gff = tmp_path / "annot.gff3"
    annot.save_gff3(out_gff)

    loaded = Annotation.from_gff3(out_gff)
    assert loaded.protein_id == "P1"
    assert len(loaded.records) == 1
    assert loaded.records[0]["type"] == "signal_peptide"
    assert loaded.records[0]["attributes"]["ID"] == "SP_1"


def test_table_parquet_io(tmp_path):
    df = pd.DataFrame([
        {"protein_id": "P1", "solubility": 0.76, "codon_cai": 0.82}
    ])
    tbl = Table(protein_id="P1", df=df)
    out_file = tmp_path / "table.parquet"
    tbl.save_parquet(out_file)

    loaded = Table.from_parquet(out_file)
    assert loaded.protein_id == "P1"
    assert loaded.df["solubility"].iloc[0] == 0.76
