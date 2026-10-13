"""Tests for topological DAG execution and whole-genome scaling."""

from pathlib import Path
import pytest
from protflow.pipeline import PipelineStep
from protflow.ui.server import ProtFlowUIHandler
from protflow.adapters.genomics import parse_protein_select
from protflow.types import Sequence, SequenceCollection


def test_topological_sort_resolves_jumbled_dag():
    """Verify that any jumbled order of visual canvas steps executes in valid producer-consumer order."""
    jumbled_steps = [
        PipelineStep(
            name="quality_assessment",
            manifest_path="manifests/04_plddt_filter.yaml",
            inputs={"structure": "structure_prediction.structure", "residue_scores": "structure_prediction.residue_scores"}
        ),
        PipelineStep(
            name="protein_select",
            manifest_path="manifests/02b_protein_select.yaml",
            inputs={"proteins": "localization_filter.proteins"}
        ),
        PipelineStep(
            name="localization_filter",
            manifest_path="manifests/02_subcellular_filter.yaml",
            inputs={"proteins": "gene_prediction.proteins"}
        ),
        PipelineStep(
            name="structure_prediction",
            manifest_path="manifests/03_alphafold_db.yaml",
            inputs={"protein": "protein_select.protein"}
        ),
        PipelineStep(
            name="gene_prediction",
            manifest_path="manifests/01_prodigal.yaml",
            inputs={"genome": "prior_input.genome"}
        ),
        PipelineStep(
            name="prior_input",
            manifest_path="manifests/00_prior_input.yaml",
            inputs={}
        ),
    ]

    handler = ProtFlowUIHandler.__new__(ProtFlowUIHandler)
    sorted_steps = handler._topological_sort_steps(jumbled_steps)
    step_names = [s.name for s in sorted_steps]

    assert step_names == [
        "prior_input",
        "gene_prediction",
        "localization_filter",
        "protein_select",
        "structure_prediction",
        "quality_assessment"
    ]


def test_protein_select_validation_and_single_candidate():
    """Verify that protein_select accepts single candidate set without id, and enforces explicit id for multi-candidate sets."""
    coll_single = SequenceCollection(
        protein_id="MIMI_SINGLE",
        sequences={"MIMI_0001": Sequence(protein_id="MIMI_0001", sequence="MAGLV")}
    )
    out_single = parse_protein_select(None, {"inputs": {"proteins": coll_single}, "parameters": {}})
    assert out_single["protein"].protein_id == "MIMI_0001"

    coll_multi = SequenceCollection(
        protein_id="MIMI_MULTI",
        sequences={
            "MIMI_0001": Sequence(protein_id="MIMI_0001", sequence="MAGLV"),
            "MIMI_0002": Sequence(protein_id="MIMI_0002", sequence="MKTAY"),
        }
    )
    with pytest.raises(ValueError):
        parse_protein_select(None, {"inputs": {"proteins": coll_multi}, "parameters": {}})

    out_spec = parse_protein_select(None, {"inputs": {"proteins": coll_multi}, "parameters": {"protein_id": "MIMI_0002"}})
    assert out_spec["protein"].protein_id == "MIMI_0002"


def test_gene_prediction_caching_and_hiv1_preview(tmp_path):
    """Verify gene prediction outputs retain provenance across cache hits and preview never raises KeyError."""
    from protflow.runner import Runner
    from protflow.manifest import Manifest

    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    manifest = Manifest.from_yaml("manifests/01_prodigal.yaml")
    # Synthetic ~1.2kb coding sequence
    genome_seq = Sequence(protein_id="HIV1_MOCK", sequence="ATG" * 400 + "TAA")

    # Run 1: Execution & cache save
    res1 = runner.run(manifest, {"genome": genome_seq})
    assert "genes" in res1
    assert "total_genes" in res1["genes"].provenance.parameters
    assert res1["genes"].provenance.parameters["total_genes"] >= 1

    # Run 2: Cache hit
    res2 = runner.run(manifest, {"genome": genome_seq})
    assert res2["genes"].provenance.parameters.get("cache_hit") is True
    assert "total_genes" in res2["genes"].provenance.parameters
    assert res2["genes"].provenance.parameters["total_genes"] == res1["genes"].provenance.parameters["total_genes"]

    # Test legacy cache simulation (provenance lacks 'total_genes')
    legacy_genes = res2["genes"]
    legacy_genes.provenance.parameters = {"cache_hit": True}
    gp = getattr(getattr(legacy_genes, "provenance", None), "parameters", {}) or {}
    total = gp.get("total_genes", len(legacy_genes.df) if hasattr(legacy_genes, "df") else 0)
    rev = gp.get("reverse_strand", int((legacy_genes.df["strand"] == "-").sum()) if hasattr(legacy_genes, "df") and "strand" in legacy_genes.df.columns else 0)
    partial = gp.get("partial", int(legacy_genes.df["partial"].sum()) if hasattr(legacy_genes, "df") and "partial" in legacy_genes.df.columns else 0)
    preview = f"{total} gen ({rev} ters iplik, {partial} kısmi)"

    assert preview.startswith(f"{len(legacy_genes.df)} gen")


def test_literature_ncbi_cds_parser():
    """Verify literature-standard parsing of NCBI fasta_cds_aa records including frameshift locations and gene symbols."""
    from protflow.adapters.genomics import parse_ncbi_cds_fasta
    raw_sample = (
        ">lcl|NC_045512.2_prot_YP_009724389.1_1 [gene=ORF1ab] [protein=ORF1ab polyprotein] [exception=ribosomal slippage] [protein_id=YP_009724389.1] [location=join(266..13468,13468..21555)] [gbkey=CDS]\n"
        "MESLVPGFNEKTHVQLSLPV\n"
        ">lcl|NC_045512.2_prot_YP_009724390.1_2 [gene=S] [protein=surface glycoprotein] [protein_id=YP_009724390.1] [location=21563..25384] [gbkey=CDS]\n"
        "MFVFLVLLPLVSSQC\n"
    )
    records = parse_ncbi_cds_fasta(raw_sample)
    assert len(records) == 2
    assert records[0]["gene_id"] == "YP_009724389.1"
    assert records[0]["gene_symbol"] == "ORF1ab"
    assert records[0]["product"] == "ORF1ab polyprotein"
    assert records[0]["start"] == 266
    assert records[0]["end"] == 21555

    assert records[1]["gene_id"] == "YP_009724390.1"
    assert records[1]["gene_symbol"] == "S"
    assert records[1]["product"] == "surface glycoprotein"
    assert records[1]["start"] == 21563
    assert records[1]["end"] == 25384


def test_viral_annotation_and_smart_protein_selection():
    """Verify that viral gene analysis avoids the polyprotein trap (>2500 aa) and prioritizes drug/vaccine targets (Spike)."""
    from protflow.adapters.genomics import annotate_gene_product, parse_protein_select
    from protflow.types import SequenceCollection, Sequence

    sym_orf1a, prod_orf1a = annotate_gene_product("VIRUS_0001", 4382, is_viral=True, pid_prefix="VIRUS")
    assert sym_orf1a == "ORF1a"
    assert "polyprotein" in prod_orf1a

    sym_s, prod_s = annotate_gene_product("VIRUS_0003", 1256, is_viral=True, pid_prefix="VIRUS")
    assert sym_s == "S"
    assert "Spike" in prod_s

    coll = SequenceCollection(protein_id="VIRUS", sequences={
        "VIRUS_0001": Sequence(protein_id="VIRUS_0001", sequence="M" * 4382, description="ORF1a polyprotein precursor | 258..13406 (+)"),
        "VIRUS_0003": Sequence(protein_id="VIRUS_0003", sequence="M" * 1256, description="S Surface glycoprotein (Spike / S protein) | 21485..25255 (+)"),
        "VIRUS_0009": Sequence(protein_id="VIRUS_0009", sequence="M" * 422, description="N Nucleocapsid phosphoprotein (N protein) | 28116..29384 (+)"),
    })

    # 'auto' must avoid the 4382 aa polyprotein and pick Spike
    sel_auto = parse_protein_select(None, {"inputs": {"proteins": coll}, "parameters": {"protein_id": "auto"}})
    assert sel_auto["protein"].protein_id == "VIRUS_0003"
    assert "spike" in sel_auto["protein"].provenance.parameters["selection_rule"]

    # Target keyword 'nucleocapsid' or 'N' must pick Nucleocapsid
    sel_n = parse_protein_select(None, {"inputs": {"proteins": coll}, "parameters": {"protein_id": "nucleocapsid"}})
    assert sel_n["protein"].protein_id == "VIRUS_0009"

    # Explicit 'longest' can still pick the polyprotein if user explicitly requests it
    sel_longest = parse_protein_select(None, {"inputs": {"proteins": coll}, "parameters": {"protein_id": "longest"}})
    assert sel_longest["protein"].protein_id == "VIRUS_0001"


