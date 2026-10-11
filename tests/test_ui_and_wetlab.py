"""
Tests for Phase 4:
- Wet Lab Alanine Scanning mutation generator and feedback loop.
- UI backend API handlers and manifest live testing.
"""

from http.server import HTTPServer
import threading
import urllib.request
import json
import pytest
import pandas as pd

from protflow.adapters.alanine_scan import (
    generate_alanine_mutants,
    record_wet_lab_feedback
)
from protflow.benchmark import CalibrationTable, load_benchmark_dataframe
from protflow.ui.server import ProtFlowUIHandler


def test_generate_alanine_mutants():
    priority_df = pd.DataFrame([
        {"protein_id": "CAND_01", "final_priority_score": 0.88},
        {"protein_id": "CAND_02", "final_priority_score": 0.76},
    ])
    pockets_df = pd.DataFrame([
        {"protein_id": "CAND_01", "residues": "A:24_TRP,A:35_GLU,A:52_ASP"},
        {"protein_id": "CAND_02", "residues": "A:12_HIS,A:45_LYS"},
    ])

    mutants_df = generate_alanine_mutants(priority_df, pockets_df, top_n=2)
    assert len(mutants_df) == 5  # 3 from CAND_01, 2 from CAND_02
    assert "mutation_code" in mutants_df.columns
    assert "A:24 TRP->ALA" in mutants_df["mutation_code"].values
    assert "A:35 GLU->ALA" in mutants_df["mutation_code"].values


def test_wet_lab_feedback_loop():
    base_bench = load_benchmark_dataframe()
    initial_count = len(base_bench)
    cal_table = CalibrationTable(base_bench)

    # Simulate wet lab experimental measurement
    lab_results = [
        {
            "protein_id": "WETLAB_HIT_01",
            "activity_validated": True,
            "measured_score": 0.92,
            "notes": "Alanine mutation at Glu35 completely abolished enzymatic activity (true catalytic residue)"
        }
    ]

    updated_df, updated_cal = record_wet_lab_feedback(lab_results, base_bench, cal_table)
    assert len(updated_df) == initial_count + 1
    assert "WETLAB_HIT_01" in updated_df["protein_id"].values
    assert not updated_cal.metrics_summary.empty


@pytest.fixture(scope="module")
def ui_test_server():
    server = HTTPServer(("127.0.0.1", 0), ProtFlowUIHandler)
    port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{port}"
    server.shutdown()


def test_ui_api_endpoints(ui_test_server):
    # 1. GET /
    req = urllib.request.urlopen(f"{ui_test_server}/")
    assert req.status == 200
    html = req.read().decode("utf-8")
    assert "ProtFlow" in html
    assert "Pipeline Studio" in html

    # 2. GET /api/organisms
    req = urllib.request.urlopen(f"{ui_test_server}/api/organisms")
    org_data = json.loads(req.read().decode("utf-8"))
    assert "organisms" in org_data
    assert len(org_data["organisms"]) >= 4

    # 3. GET /api/manifests
    req = urllib.request.urlopen(f"{ui_test_server}/api/manifests")
    data = json.loads(req.read().decode("utf-8"))
    assert "manifests" in data
    assert len(data["manifests"]) >= 7

    # 3. GET /api/pipeline
    req = urllib.request.urlopen(f"{ui_test_server}/api/pipeline")
    pdata = json.loads(req.read().decode("utf-8"))
    assert pdata["name"] == "bacterial_enzyme_prioritization"
    assert len(pdata["steps"]) >= 8

    # 4. POST /api/manifest/test (Live type verification)
    sample_manifest = {
        "manifest": {
            "name": "mock_validator",
            "version": "1.0.0",
            "inputs": {"structure": "Structure"},
            "outputs": {"pockets": "Pocket"},
            "execution": {"type": "local", "command": "python -c \"pass\""},
            "adapter": "protflow.adapters.pockets:parse_p2rank"
        }
    }
    req_post = urllib.request.Request(
        f"{ui_test_server}/api/manifest/test",
        data=json.dumps(sample_manifest).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    res = urllib.request.urlopen(req_post)
    test_res = json.loads(res.read().decode("utf-8"))
    assert test_res["status"] == "PASS"
    assert "pockets" in test_res["outputs"]
    assert "logs" in test_res
    assert len(test_res["logs"]) >= 3

    # 5. POST /api/manifest/test on filter node (verify filter drop counters)
    filter_manifest = {
        "manifest": {
            "name": "mock_filter_node",
            "version": "1.0.0",
            "inputs": {"protein": "Sequence"},
            "outputs": {"protein": "Sequence", "filter_table": "Table"},
            "execution": {"type": "local", "command": "python -c \"pass\""},
            "adapter": "protflow.adapters.localization:parse_subcellular_filter"
        }
    }
    req_filter = urllib.request.Request(
        f"{ui_test_server}/api/manifest/test",
        data=json.dumps(filter_manifest).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    res_f = urllib.request.urlopen(req_filter)
    f_data = json.loads(res_f.read().decode("utf-8"))
    assert f_data["status"] == "PASS"
    assert f_data["filter_stats"] is not None
    assert "input_count" in f_data["filter_stats"]
    assert "passed_count" in f_data["filter_stats"]

    # 6. POST /api/manifest/save validation check
    invalid_save = {
        "filename": "invalid_test.yaml",
        "raw_yaml": "name: bad\ninputs:\n  p: NonExistentType\n"
    }
    req_bad_save = urllib.request.Request(
        f"{ui_test_server}/api/manifest/save",
        data=json.dumps(invalid_save).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with pytest.raises(urllib.error.HTTPError) as exc_info:
        urllib.request.urlopen(req_bad_save)
    assert exc_info.value.code == 400

    # 7. Verify UI features: emoji-free, interactive wire hitboxes, filter counters, terminal box, LOD
    import re
    assert not re.search(r"[\U00010000-\U0010ffff\u2600-\u27bf]", html), "UI contains emoji characters"
    assert "wire-hitbox" in html
    assert "deleteConnection" in html
    assert "deleteSelected" in html
    assert "overflow: visible !important;" in html
    assert "filter-stat-badge" in html
    assert "terminal-logs" in html
    assert "lod-compact" in html
    assert "swapToolForNode" in html


def test_ui_server_security_hardening(ui_test_server):
    import urllib.error

    # 1. Non-JSON Content-Type rejection (415)
    req_bad_ct = urllib.request.Request(
        f"{ui_test_server}/api/manifest/test",
        data=b"test data",
        headers={"Content-Type": "text/plain"}
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(req_bad_ct)
    assert exc.value.code == 415

    # 2. Path traversal attack on static files (403)
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(f"{ui_test_server}/../../../../Windows/win.ini")
    assert exc.value.code == 403

    # 3. Untrusted Origin header rejection (403)
    req_bad_origin = urllib.request.Request(
        f"{ui_test_server}/api/manifest/test",
        data=json.dumps({"manifest": {}}).encode("utf-8"),
        headers={"Content-Type": "application/json", "Origin": "http://evil-attacker.com"}
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(req_bad_origin)
    assert exc.value.code == 403


def test_ui_pipeline_run_with_prior_input(ui_test_server):
    """Verify that /api/pipeline/run accepts dynamic steps with manifest_path and custom genome input."""
    payload = {
        "genome_dna": "ATGGCTGGGCTGGTGGCCATCCTGGAGAAATTCGGCGACGAGCTGCGTTGGGCGCCGAAAACCTGCCACTATGGCGCGGTGATCCTGTCCTAA",
        "organism": "ecoli_k12",
        "pipeline": {
            "steps": [
                {
                    "name": "prior_input",
                    "manifest_path": "manifests/00_prior_input.yaml",
                    "inputs": {},
                    "parameters": {"organism": "ecoli_k12"}
                },
                {
                    "name": "gene_prediction",
                    "manifest_path": "manifests/01_prodigal.yaml",
                    "inputs": {"genome": "genome"},
                    "parameters": {}
                }
            ]
        }
    }
    req = urllib.request.Request(
        f"{ui_test_server}/api/pipeline/run",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    res = urllib.request.urlopen(req)
    data = json.loads(res.read().decode("utf-8"))
    assert data["status"] == "SUCCESS"
    assert len(data["step_results"]) == 2
    assert data["step_results"][0]["name"] == "prior_input"
    assert data["step_results"][0]["status"] in ("completed", "cached")
    assert "Genom:" in data["step_results"][0]["preview"]
    assert data["step_results"][1]["name"] == "gene_prediction"
    assert "Sekans:" in data["step_results"][1]["preview"]
    assert "node_outputs" in data
    assert "prior_input" in data["node_outputs"]
    assert "genome" in data["node_outputs"]["prior_input"]
    assert data["node_outputs"]["prior_input"]["genome"]["type"] == "Sequence"
    assert "gene_prediction" in data["node_outputs"]
    assert len(data["priority_table"]) >= 1
    assert len(data["alanine_mutants"]) >= 1


def test_ui_ncbi_fetch_endpoint(ui_test_server):
    """Verify that /api/ncbi/fetch validates queries and pulls genome data."""
    # 1. Empty query rejection (400)
    req_bad = urllib.request.Request(
        f"{ui_test_server}/api/ncbi/fetch",
        data=json.dumps({"query": ""}).encode("utf-8"),
        headers={"Content-Type": "application/json"}
    )
    with pytest.raises(urllib.error.HTTPError) as exc:
        urllib.request.urlopen(req_bad)
    assert exc.value.code == 400

    # 2. Live / Network-resilient query test
    try:
        req_ok = urllib.request.Request(
            f"{ui_test_server}/api/ncbi/fetch",
            data=json.dumps({"query": "GCF_000864885.1"}).encode("utf-8"),
            headers={"Content-Type": "application/json"}
        )
        res = urllib.request.urlopen(req_ok, timeout=12)
        n_data = json.loads(res.read().decode("utf-8"))
        assert n_data["status"] == "SUCCESS"
        assert "result" in n_data
        assert n_data["result"]["length"] > 1000
        assert n_data["result"]["gc_content"] > 0
    except (urllib.error.URLError, TimeoutError):
        pass  # Graceful fallback if offline


