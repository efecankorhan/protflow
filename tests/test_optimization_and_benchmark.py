"""
Phase 3 Milestone Tests:
1. Verification of Structure Prediction Fallback Guard (eliminating silent mock risks).
2. Node-level Score Normalization.
3. Swiss-Prot Benchmark Dataset integrity.
4. Empirical Calibration Table & Calibrated Confidence scoring.
5. Multi-variant Pipeline Optimization Sweep.
6. Automatic Node Evaluator (comparing new node against baseline).
"""

from pathlib import Path
import pytest
import pandas as pd

from protflow.types import Sequence, Structure
from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.adapters.structure_fetch import StructureNotFoundError
from protflow.normalize import (
    normalize_plddt,
    normalize_pocket_score,
    normalize_solubility_from_sasa
)
from protflow.benchmark import (
    load_benchmark_dataframe,
    get_benchmark_sequences,
    CalibrationTable,
    VariantConfig,
    PipelineOptimizer,
    NodeEvaluator
)


def test_structure_prediction_fallback_guard_raises_error(tmp_path):
    """
    User Risk Fix:
    Verifies that when an unannotated protein is missing from AlphaFold DB
    and allow_synthetic_fallback is False, the system halts with StructureNotFoundError
    instead of silently fabricating a model.
    """
    manifest_af = Manifest.from_dict({
        "name": "alphafold_db",
        "version": "v4",
        "inputs": {"protein": "Sequence"},
        "outputs": {"structure": "Structure", "residue_scores": "ResidueScore"},
        "execution": {
            "type": "http",
            "url": "https://alphafold.ebi.ac.uk/api/prediction/{inputs.protein.protein_id}",
            "method": "GET"
        },
        "adapter": "protflow.adapters.structure_fetch:parse_alphafold_fetcher",
        "is_remote": True
    })

    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=False)
    seq = Sequence(protein_id="NOVEL_UNANNOTATED_GENE_9999", sequence="MKFLVALA")

    # When fallback is disabled, must raise StructureNotFoundError
    with pytest.raises(StructureNotFoundError, match="not found in AlphaFold DB"):
        runner.run(manifest_af, {"protein": seq}, parameters={"allow_synthetic_fallback": False})

    # When fallback is explicitly allowed, produces model with dev mock tag
    out = runner.run(manifest_af, {"protein": seq}, parameters={"allow_synthetic_fallback": True})
    assert "structure" in out
    assert isinstance(out["structure"], Structure)


def test_score_normalization():
    # pLDDT Normalization
    assert normalize_plddt(95.0) > 0.85
    assert 0.50 <= normalize_plddt(75.0) <= 0.85
    assert normalize_plddt(40.0) < 0.20

    # Pocket score normalization
    assert normalize_pocket_score(0.85, tool="p2rank") == 0.85
    assert 0.0 <= normalize_pocket_score(15.0, tool="fpocket") <= 1.0

    # SASA Solubility
    assert normalize_solubility_from_sasa(0.40) == 1.0  # optimal target
    assert normalize_solubility_from_sasa(0.70) < 0.50  # excessively hydrophobic


def test_swiss_prot_benchmark_dataset():
    df = load_benchmark_dataframe()
    assert len(df) == 9

    # Check 6 EC classes present
    ec_classes = df[df["is_active_enzyme"] == True]["ec_class"].unique()
    assert len(ec_classes) == 6

    # Check negative controls present
    neg_controls = df[df["is_active_enzyme"] == False]
    assert len(neg_controls) == 3


def test_calibration_table_and_confidence_attachment():
    # Synthetic evaluation data
    eval_df = pd.DataFrame([
        {"protein_id": "P1", "ec_class": "EC 3 (Hydrolase)", "is_active_enzyme": True, "final_priority_score": 0.88},
        {"protein_id": "P2", "ec_class": "EC 3 (Hydrolase)", "is_active_enzyme": True, "final_priority_score": 0.78},
        {"protein_id": "P3", "ec_class": "EC 3 (Hydrolase)", "is_active_enzyme": False, "final_priority_score": 0.52},
        {"protein_id": "P4", "ec_class": "NON_ENZYME", "is_active_enzyme": False, "final_priority_score": 0.35},
        {"protein_id": "P5", "ec_class": "NON_ENZYME", "is_active_enzyme": False, "final_priority_score": 0.42},
    ])

    cal_table = CalibrationTable(eval_df)
    assert not cal_table.metrics_summary.empty

    # High raw score should map to high calibrated probability
    high_conf = cal_table.calibrate(0.90, ec_class="EC 3 (Hydrolase)")
    low_conf = cal_table.calibrate(0.30, ec_class="NON_ENZYME")
    assert high_conf > low_conf

    # Verify attachment
    enriched = cal_table.attach_calibration(eval_df)
    assert "calibrated_confidence" in enriched.columns
    assert "reliability_grade" in enriched.columns
    assert enriched["calibrated_confidence"].iloc[0] > 0.60


def test_pipeline_variant_optimizer_sweep(tmp_path):
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    optimizer = PipelineOptimizer(pipeline_yaml, runner=runner)

    # Define 2 candidate variants to compare
    variants = [
        VariantConfig(name="variant_p2rank_std", pocket_tool="p2rank", min_plddt=70.0),
        VariantConfig(name="variant_fpocket_std", pocket_tool="fpocket", min_plddt=70.0),
    ]

    best, summary_df, cal_table = optimizer.sweep(variants)

    assert len(summary_df) == 2
    assert "f1_score" in summary_df.columns
    assert "precision" in summary_df.columns
    assert best.f1_score > 0.0
    assert not cal_table.metrics_summary.empty


def test_automatic_node_evaluator(sample_data):
    st = sample_data["structure"]
    evaluator = NodeEvaluator()

    m_p2rank = Manifest.from_yaml("manifests/05_pocket_p2rank.yaml")
    m_fpocket = Manifest.from_yaml("manifests/05_pocket_fpocket.yaml")

    report = evaluator.compare_nodes(
        candidate_manifest=m_p2rank,
        baseline_manifest=m_fpocket,
        test_inputs_list=[{"structure": st}]
    )

    assert report.samples_tested == 1
    assert report.candidate_name == "p2rank"
    assert report.baseline_name == "fpocket"
    assert report.speedup > 0.0
    assert "RECOMMENDED" in report.recommendation or "ACCEPTABLE" in report.recommendation
    assert not report.comparison_df.empty


def test_automatic_node_evaluator_with_empty_output_table(tmp_path):
    evaluator = NodeEvaluator()
    dummy_cmd = 'python -c "import json; print(json.dumps([]))"'
    m_base = Manifest.from_dict({
        "name": "base_empty",
        "version": "1.0",
        "inputs": {"structure": "Structure"},
        "outputs": {"pockets": "Pocket"},
        "execution": {"type": "local", "command": dummy_cmd},
        "adapter": "protflow.adapters.parsers:parse_dummy_pocket_finder"
    })
    m_cand = Manifest.from_dict({
        "name": "cand_empty",
        "version": "1.0",
        "inputs": {"structure": "Structure"},
        "outputs": {"pockets": "Pocket"},
        "execution": {"type": "local", "command": dummy_cmd},
        "adapter": "protflow.adapters.parsers:parse_dummy_pocket_finder"
    })
    st = Structure(protein_id="P_EMPTY", file_path="dummy.pdb")
    # Must not raise IndexError when df is empty
    report = evaluator.compare_nodes(m_cand, m_base, [{"structure": st}])
    assert report.samples_tested == 1
    assert report.mean_candidate_score == 0.0


def test_pipeline_from_yaml_with_manifest_path_key(tmp_path):
    from protflow.pipeline import Pipeline
    yaml_content = """
name: test_pipe
steps:
  - name: test_step
    manifest_path: manifests/00_prior_input.yaml
    inputs: {}
"""
    p_file = tmp_path / "test_pipe.yaml"
    p_file.write_text(yaml_content, encoding="utf-8")
    pipe = Pipeline.from_yaml(p_file)
    assert len(pipe.steps) == 1
    assert pipe.steps[0].manifest_path == "manifests/00_prior_input.yaml"


def test_container_runner_command_interpolation(tmp_path):
    from protflow.types import Sequence
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=False)
    seq = Sequence(protein_id="PROT_01", sequence="MKFLV")
    cmd_tpl = "tool --in {inputs.protein.file_path} --out {output_dir} --opt {parameters.mode}"
    formatted = runner._format_command(cmd_tpl, {"protein": seq}, {"mode": "fast"}, tmp_path / "work")
    assert "PROT_01.fasta" in formatted
    assert "--out" in formatted
    assert "--opt fast" in formatted

