"""
Phase 2 Milestone Test: End-to-End Reference Bacterial Pipeline.
Executes Genome -> Prodigal -> Subcellular -> AFDB -> pLDDT filter -> Pocket -> SASA -> Foldseek -> Prioritization.
Tests single-line node interchangeability (P2Rank <-> fpocket) and full pipeline export.
"""

from pathlib import Path
import pytest
import pandas as pd

from protflow.pipeline import Pipeline
from protflow.runner import Runner
from protflow.types import Sequence, Table, Pocket, Structure, ResidueScore


@pytest.fixture
def bacterial_genome():
    # Synthetic bacterial DNA fragment with an active enzyme ORF (ATG ... TAA)
    # Codes for: M-A-G-L-V-A-I-L-E-K-F-G-D-E-L-R-W-A-P-K-T-C-H-Y-G-A-V-I-L-S-*
    dna = (
        "ATGGCTGGGCTGGTGGCCATCCTGGAGAAATTCGGCGACGAGCTGCGTTGGGCGCCGAAAACCTGCCACTATGGCGCGGTGATCCTGTCCTAA"
    )
    return Sequence(protein_id="ECOLI_MG1655_FRAG1", sequence=dna, description="E. coli K-12 fragment")


def test_bacterial_reference_pipeline_end_to_end(tmp_path, bacterial_genome):
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)

    initial_state = {"genome": bacterial_genome}
    final_state = pipeline.run(initial_state)

    # 1. Prodigal output
    assert "protein" in final_state
    protein = final_state["protein"]
    assert isinstance(protein, Sequence)
    assert len(protein.sequence) > 10

    # 2. Subcellular localization
    assert "filter_table" in final_state
    loc_df = final_state["filter_table"].df
    assert bool(loc_df["is_intracellular"].iloc[0]) is True

    # 3. AlphaFold DB fetch
    assert "structure" in final_state
    assert isinstance(final_state["structure"], Structure)
    assert "residue_scores" in final_state
    assert isinstance(final_state["residue_scores"], ResidueScore)

    # Validate true 3D spatial extent (not flat 1-D on x-axis with y=0, z=0)
    st = final_state["structure"].to_gemmi_structure()
    ca_atoms = [atom for model in st for chain in model for res in chain for atom in res if atom.name == "CA"]
    assert len(ca_atoms) > 0
    ys = [atom.pos.y for atom in ca_atoms]
    zs = [atom.pos.z for atom in ca_atoms]
    assert any(abs(y) > 0.1 for y in ys), "Structure Y-coordinates must not be all 0 (1-D flat line)"
    assert any(abs(z) > 0.1 for z in zs), "Structure Z-coordinates must not be all 0 (1-D flat line)"

    # Validate realistic non-flat residue pLDDT variation
    score_df = final_state["residue_scores"].df
    assert not score_df.empty
    assert score_df["score"].std() > 0.0, "pLDDT residue scores must have variance across residues, not flat constant"

    # 4. pLDDT Quality filter
    assert "quality_table" in final_state
    quality_df = final_state["quality_table"].df
    assert bool(quality_df["passed_plddt_filter"].iloc[0]) is True
    assert quality_df["mean_plddt"].iloc[0] >= 70.0

    # 5. Pocket detection (P2Rank)
    assert "pockets" in final_state
    pockets = final_state["pockets"]
    assert isinstance(pockets, Pocket)
    assert pockets.df["score"].iloc[0] > 0.0
    assert pockets.df["volume"].iloc[0] > 0.0

    # 6. Surface SASA
    assert "sasa_table" in final_state
    sasa_df = final_state["sasa_table"].df
    assert sasa_df["total_sasa"].iloc[0] > 0.0

    # 7. Foldseek Structural homology & novelty
    assert "foldseek_table" in final_state
    foldseek_df = final_state["foldseek_table"].df
    assert foldseek_df["tm_score"].iloc[0] > 0.0
    assert foldseek_df["novelty_score"].iloc[0] >= 0.0

    # 8. Final Prioritization Table
    assert "priority_table" in final_state
    priority_table = final_state["priority_table"]
    assert isinstance(priority_table, Table)
    p_df = priority_table.df

    assert "final_priority_score" in p_df.columns
    assert "structure_quality" in p_df.columns
    assert "pocket_score" in p_df.columns
    assert "solubility_est" in p_df.columns
    assert "novelty_score" in p_df.columns
    assert "status" in p_df.columns

    final_score = p_df["final_priority_score"].iloc[0]
    assert 0.0 <= final_score <= 1.0
    assert p_df["status"].iloc[0] in ("CANDIDATE", "DEPURGED")


def test_single_line_pocket_node_interchangeability(tmp_path, bacterial_genome):
    """
    Measures the core principle:
    Swapping a pocket predictor (P2Rank -> fpocket) is literally a 1-line configuration change.
    """
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)

    # Find pocket_finding step and switch its manifest to fpocket
    pocket_step = next(s for s in pipeline.steps if s.name == "pocket_finding")
    assert "p2rank" in pocket_step.manifest_path

    # Swap in 1 line
    pocket_step.manifest_path = "manifests/05_pocket_fpocket.yaml"

    final_state = pipeline.run({"genome": bacterial_genome})

    # Assert fpocket executed smoothly and returned typed Pocket
    assert "pockets" in final_state
    pockets = final_state["pockets"]
    assert isinstance(pockets, Pocket)
    assert pockets.provenance.tool_name == "fpocket"
    assert "priority_table" in final_state
    assert final_state["priority_table"].df["pocket_score"].iloc[0] > 0.0


def test_single_line_gene_prediction_node_interchangeability(tmp_path, bacterial_genome):
    """
    Measures the core principle:
    Swapping a gene predictor (augustus <-> prodigal) is literally a 1-line configuration change.
    """
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)

    gene_step = next(s for s in pipeline.steps if s.name == "gene_prediction")
    assert "augustus" in gene_step.manifest_path

    # Swap in 1 line
    gene_step.manifest_path = "manifests/01_prodigal.yaml"

    final_state = pipeline.run({"genome": bacterial_genome})
    assert "protein" in final_state
    assert final_state["genes"].provenance.tool_name == "prodigal"
    assert "structure" in final_state


def test_pipeline_export_to_snakemake_and_nextflow(tmp_path):
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    pipeline = Pipeline.from_yaml(pipeline_yaml)

    # 1. Snakemake export
    smk = pipeline.export_snakemake()
    assert ("rule prodigal:" in smk or "rule augustus:" in smk)
    assert "rule alphafold_db:" in smk
    assert "rule p2rank:" in smk
    assert "rule prioritizer:" in smk

    # 2. Nextflow export
    nf = pipeline.export_nextflow()
    assert "nextflow.enable.dsl=2" in nf
    assert ("process PRODIGAL {" in nf or "process AUGUSTUS {" in nf)
    assert "process ALPHAFOLD_DB {" in nf
    assert "process P2RANK {" in nf
    assert "process PRIORITIZER {" in nf


def test_mimivirus_translation_and_gene_prediction():
    from protflow.adapters.genomics import translate_dna, parse_or_predict_genes, parse_prior_input

    # 1. Verification of the exact 102 bp fragment from Mimivirus shirakomae (nt 2167-2268)
    nt_seq = "TTGAAAATCTGGTGTCACATAAATAAAATATTCATCAGCAAATTTTGTGATGAATCCACCAGTTCTATATATATCTTGTTTATTGGTGTAAGATACTATTGA"
    protein = translate_dna(nt_seq, is_cds=True).rstrip("*")
    # Must start with M (initiation codon TTG -> M), contain C at pos 5 and 16, and zero 'X'
    assert protein == "MKIWCHINKIFISKFCDESTSSIYILFIGVRYY"
    assert "X" not in protein
    assert protein.startswith("M")

    # 2. Test parse_prior_input with Mimivirus
    prior = parse_prior_input(None, {"parameters": {"organism": "mimivirus_shirakomae", "genome_dna": nt_seq}})
    assert prior["genome"].protein_id == "MIMIVIRUS_SHIRAKOMAE"
    assert "ECOLI" not in prior["genome"].protein_id

    # 3. Test parse_or_predict_genes
    genes = parse_or_predict_genes(None, {"inputs": {"genome": prior["genome"]}})
    assert genes["genes"].df["gene_id"].iloc[0].startswith("MIMIVIRUS_SHIRAKOMAE_")
    assert genes["annotations"].records[0]["source"] == "pyrodigal"
    assert set(genes["genes"].df.columns) >= {"start", "end", "strand", "partial", "protein", "cds_nt"}
    assert len(genes["proteins"].sequences) == len(genes["genes"].df)

    # 4. protein_select: explicit id works, unknown id and ambiguous set fail loudly
    from protflow.adapters.genomics import parse_protein_select
    ids = list(genes["proteins"].sequences)
    got = parse_protein_select(None, {"inputs": {"proteins": genes["proteins"]}, "parameters": {"protein_id": ids[0]}})
    assert got["protein"].protein_id == ids[0]
    # 5. auto selection picks the longest ORF with provenance metadata
    from protflow.types import SequenceCollection
    multi_coll = SequenceCollection(protein_id="multi", sequences={
        "short_p": Sequence(protein_id="short_p", sequence="M" * 10),
        "long_p": Sequence(protein_id="long_p", sequence="M" * 50),
    })
    auto_got = parse_protein_select(None, {"inputs": {"proteins": multi_coll}, "parameters": {"protein_id": "auto"}})
    assert auto_got["protein"].protein_id == "long_p"
    assert auto_got["protein"].provenance.parameters["auto_selected"] is True
    assert auto_got["protein"].provenance.parameters["selection_rule"] == "longest_orf"

    # 'first' selection picks the literal first sequence
    first_got = parse_protein_select(None, {"inputs": {"proteins": multi_coll}, "parameters": {"protein_id": "first"}})
    assert first_got["protein"].protein_id == "short_p"
    assert first_got["protein"].provenance.parameters["selection_rule"] == "first_orf"

    # null parameters does not raise AttributeError
    single_coll = SequenceCollection(protein_id="single", sequences={"p1": Sequence(protein_id="p1", sequence="M" * 12)})
    null_got = parse_protein_select(None, {"inputs": {"proteins": single_coll}, "parameters": None})
    assert null_got["protein"].protein_id == "p1"


def test_filter_acts_on_whole_gene_set():
    from protflow.adapters.localization import parse_subcellular_filter
    from protflow.types import SequenceCollection

    secreted = "MLLLVVAAFFLLIIGGAAVVLLAAKK"  # hydrophobic N-terminus -> dropped
    cyto = "MKTEEDKRSTNPQQDDEK" * 3
    coll = SequenceCollection(protein_id="G", sequences={
        "a": Sequence(protein_id="a", sequence=secreted),
        "b": Sequence(protein_id="b", sequence=cyto),
    })
    out = parse_subcellular_filter(None, {"inputs": {"proteins": coll}})
    assert list(out["proteins"].sequences) == ["b"]
    assert "protein" not in out
    p = out["filter_table"].provenance.parameters
    assert (p["input_count"], p["passed_count"], p["dropped_count"]) == (2, 1, 1)


def test_single_secreted_sequence_dropped():
    from protflow.adapters.localization import parse_subcellular_filter
    from protflow.types import Sequence

    secreted = Sequence(protein_id="sec_1", sequence="MLLLVVAAFFLLIIGGAAVVLLAAKK" + "A" * 50)
    out = parse_subcellular_filter(None, {"inputs": {"protein": secreted}})
    assert out["protein"].sequence == ""
    assert out["protein"].description == "DROPPED_BY_FILTER"
    df = out["filter_table"].df
    assert df["status"].iloc[0] == "DROPPED"
    assert bool(df["is_intracellular"].iloc[0]) is False
    assert df["drop_reason"].iloc[0] == "Signal peptide (Extracellular)"


def test_collection_all_secreted_dropped():
    from protflow.adapters.localization import parse_subcellular_filter
    from protflow.types import Sequence, SequenceCollection

    secreted1 = Sequence(protein_id="sec_1", sequence="MLLLVVAAFFLLIIGGAAVVLLAAKK" + "A" * 50)
    secreted2 = Sequence(protein_id="sec_2", sequence="MLLLVVAAFFLLIIGGAAVVLLAAKK" + "G" * 50)
    coll = SequenceCollection(protein_id="all_sec", sequences={"sec_1": secreted1, "sec_2": secreted2})
    out = parse_subcellular_filter(None, {"inputs": {"proteins": coll}})
    assert len(out["proteins"].sequences) == 0
    p = out["filter_table"].provenance.parameters
    assert p["input_count"] == 2
    assert p["passed_count"] == 0
    assert p["dropped_count"] == 2


@pytest.mark.parametrize("org_id", [
    "ecoli_k12",
    "brachionus_plicatilis",
    "sars_cov_2",
])
def test_multi_organism_pipeline_execution(tmp_path, org_id):
    """
    Validates end-to-end pipeline execution across diverse biological kingdoms:
    - Bacteria (E. coli K-12)
    - Eukaryota (Brachionus plicatilis)
    - Viruses (SARS-CoV-2)
    Verifies gene prediction, protein selection, true 3D tertiary structure folding,
    pocket detection, and prioritization score across all domains of life.
    """
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / f"cache_{org_id}", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)
    pipeline.steps[0].parameters["organism"] = org_id

    final_state = pipeline.run({})

    assert "protein" in final_state
    protein = final_state["protein"]
    assert len(protein.sequence) > 10

    assert "structure" in final_state
    st = final_state["structure"].to_gemmi_structure()
    ca_atoms = [atom for model in st for chain in model for res in chain for atom in res if atom.name == "CA"]
    assert len(ca_atoms) > 0
    ys = [atom.pos.y for atom in ca_atoms]
    zs = [atom.pos.z for atom in ca_atoms]
    assert any(abs(y) > 0.1 for y in ys), f"{org_id}: Y coordinates must not be all 0 (1-D line)"
    assert any(abs(z) > 0.1 for z in zs), f"{org_id}: Z coordinates must not be all 0 (1-D line)"

    assert "residue_scores" in final_state
    scores = final_state["residue_scores"].df
    assert not scores.empty
    assert scores["score"].std() > 0.0, f"{org_id}: residue pLDDT scores must have variance across residues"

    assert "pockets" in final_state
    assert not final_state["pockets"].df.empty
    assert final_state["pockets"].df["score"].iloc[0] > 0.0

    assert "priority_table" in final_state
    p_df = final_state["priority_table"].df
    assert 0.0 <= p_df["final_priority_score"].iloc[0] <= 1.0


def test_random_organism_pipeline(tmp_path):
    """
    Randomly samples an organism from the available organism library
    to continuously test pipeline robustness against variable biological data.
    """
    import random
    organisms_dir = Path(__file__).parent.parent / "organisms"
    available = [p.stem for p in organisms_dir.glob("*.json")]
    assert len(available) >= 3
    chosen_org = random.choice(available)

    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / f"cache_random_{chosen_org}", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)
    pipeline.steps[0].parameters["organism"] = chosen_org

    final_state = pipeline.run({})
    assert "structure" in final_state
    assert "priority_table" in final_state
    assert 0.0 <= final_state["priority_table"].df["final_priority_score"].iloc[0] <= 1.0


def test_generic_id_uniprot_resolution_safety():
    from protflow.adapters.structure_fetch import resolve_uniprot_afdb
    # Generic IDs like 'g1', 'orf1', 'gene_0001' must not falsely match random UniProt entries
    assert resolve_uniprot_afdb("g1") is None
    assert resolve_uniprot_afdb("orf1") is None
    assert resolve_uniprot_afdb("gene_0001") is None


def test_negative_strand_cds_reverse_complement(tmp_path):
    from protflow.adapters.genomics import parse_augustus, translate_dna
    # Synthetic DNA with a negative strand ORF (TAA ... CAT on forward, so ATG ... TTA on reverse)
    # Forward: 5'-TTA GTT GAT GTA CAT-3' -> Reverse: 5'-ATG TAC ATC AAC TAA-3' (M-Y-I-N-*)
    dna = "TTAGTTGATGTACAT"
    seq_obj = Sequence(protein_id="NEG_TEST", sequence=dna)
    gtf_raw = (
        "NEG_TEST\tAUGUSTUS\tgene\t1\t15\t0.99\t-\t.\tg1\n"
        "NEG_TEST\tAUGUSTUS\ttranscript\t1\t15\t0.99\t-\t.\tg1.t1\n"
        "NEG_TEST\tAUGUSTUS\tCDS\t1\t15\t0.99\t-\t0\ttranscript_id \"g1.t1\"; gene_id \"g1\";\n"
        "# protein sequence = [MYIN]\n"
    )
    out = parse_augustus(gtf_raw, {"inputs": {"genome": seq_obj}, "output_dir": str(tmp_path)})
    genes_df = out["genes"].df
    assert not genes_df.empty
    row = genes_df.iloc[0]
    assert row["strand"] == "-"
    # cds_nt must be reverse-complemented: ATG TAC ATC AAC TAA (starts with ATG)
    assert row["cds_nt"] == "ATGTACATCAACTAA"
    assert row["start_codon"] == "ATG"
    assert translate_dna(row["cds_nt"][:12]) == "MYIN"


def test_empty_sequence_localization_filter():
    from protflow.adapters.localization import parse_subcellular_filter
    # Empty inputs must not raise KeyError or AttributeError
    res = parse_subcellular_filter(None, {"inputs": {}})
    assert res["protein"].sequence == ""
    assert res["filter_table"].df.empty


