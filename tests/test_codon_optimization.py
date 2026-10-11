"""
Tests for Codon Optimization and Host Profile Adaptation:
1. Translates protein into host-optimized codons (E. coli, Human/CHO, Yeast).
2. Avoids forbidden motifs (cryptic splice sites, SD pauses, BsaI/EcoRI restriction sites).
3. CAI calculation and %GC tracking.
4. Manifest runner execution through ProtFlow engine.
5. Multi-format export package (.fasta, .gbk, .csv, .parquet).
"""

from pathlib import Path
import pytest
import pandas as pd

from protflow.types import Sequence, Table
from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.adapters.codon import (
    HOST_PROFILES,
    compute_cai,
    optimize_protein_sequence,
    parse_codon_optimization
)
from protflow.export import export_sequence_to_genbank, export_construct_package


def test_codon_optimization_ecoli():
    # Protein sequence with amino acids
    prot = Sequence(protein_id="ENZ_01", sequence="MKFLVALA", description="Test enzyme")
    
    dna, df, stats = optimize_protein_sequence(prot.sequence, host="ecoli")
    
    assert len(dna) == (len(prot.sequence) + 1) * 3  # including stop codon
    assert stats["cai_score"] >= 0.85
    assert 45.0 <= stats["gc_percent"] <= 60.0
    assert len(df) == len(prot.sequence) + 1
    assert "codon" in df.columns
    assert "gc_percent" in df.columns


def test_motif_avoidance_breaks_forbidden_sites():
    # Sequence of Leucines which naturally produces repetitive optimal CTG codons
    prot_seq = "LLLL"

    dna, df, stats = optimize_protein_sequence(
        protein_seq=prot_seq,
        host="ecoli",
        extra_forbidden=["CTGCTG"]
    )

    assert "CTGCTG" not in dna
    assert stats["motifs_avoided"] >= 1


def test_codon_manifest_runner_execution(tmp_path):
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    manifest = Manifest.from_yaml("manifests/10_codon_optimization.yaml")
    
    prot = Sequence(protein_id="TEST_PROT", sequence="MKFLVALAEG", description="Enzyme construct")
    outputs = runner.run(manifest, {"protein": prot}, output_dir=tmp_path / "out")
    
    assert "optimized_cds" in outputs
    assert "codon_table" in outputs
    assert isinstance(outputs["optimized_cds"], Sequence)
    assert isinstance(outputs["codon_table"], Table)
    
    cds = outputs["optimized_cds"]
    assert cds.sequence.startswith("ATG")
    assert "CAI=" in cds.description
    
    table_df = outputs["codon_table"].df
    assert "cai_score" in table_df.columns
    assert "total_gc" in table_df.columns


def test_export_construct_package(tmp_path):
    cds = Sequence(protein_id="CONSTRUCT_01", sequence="ATGAAGTTTCTGGTTGCACTGGCAGAAGGTTAA", description="Optimized gene")
    df = pd.DataFrame([{"pos": 1, "aa": "M", "codon": "ATG", "cai": 1.0}])
    table = Table(protein_id="CONSTRUCT_01", df=df)
    
    files = export_construct_package(cds, table, tmp_path / "exports")
    
    assert files["fasta"].exists()
    assert files["genbank"].exists()
    assert files["csv"].exists()
    assert files["parquet"].exists()
    
    gbk_content = files["genbank"].read_text(encoding="utf-8")
    assert "LOCUS" in gbk_content
    assert "ORIGIN" in gbk_content
    assert "CONSTRUCT_01" in gbk_content
