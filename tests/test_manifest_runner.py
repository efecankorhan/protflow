"""
Phase 1 Milestone Test:
Verifies that the same tool defined in 3 different execution forms
(local path, container, and web service) runs uniformly through the engine
and produces identical typed outputs with provenance and privacy safeguards.
"""

from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path
import json
import threading
import pytest

from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.types import Pocket, Structure


class MockHTTPPocketServer(BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        payload = [
            {
                "pocket_id": 1,
                "score": 0.88,
                "center_x": 1.0,
                "center_y": 2.0,
                "center_z": 3.0,
                "volume": 350.0,
                "residues": "A:1,A:2"
            }
        ]
        self.wfile.write(json.dumps(payload).encode("utf-8"))

    def log_message(self, format, *args):
        # Silence standard HTTP request logs
        return


@pytest.fixture(scope="module")
def mock_server():
    server = HTTPServer(("127.0.0.1", 0), MockHTTPPocketServer)
    port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{port}"
    server.shutdown()


def test_three_execution_modes_produce_identical_typed_output(tmp_path, sample_data, mock_server):
    structure = sample_data["structure"]
    inputs = {"structure": structure}

    # Python command generating dummy pocket JSON
    dummy_cmd = (
        'python -c "import json; print(json.dumps([{'
        '\'pocket_id\': 1, \'score\': 0.88, \'center_x\': 1.0, '
        '\'center_y\': 2.0, \'center_z\': 3.0, \'volume\': 350.0, \'residues\': \'A:1,A:2\''
        '}]))"'
    )

    # 1. LOCAL MANIFEST
    manifest_local = Manifest.from_dict({
        "name": "pocket_predictor",
        "version": "1.0.0",
        "inputs": {"structure": "Structure"},
        "outputs": {"pockets": "Pocket"},
        "execution": {
            "type": "local",
            "command": dummy_cmd
        },
        "adapter": "protflow.adapters.parsers:parse_dummy_pocket_finder"
    })

    # 2. CONTAINER MANIFEST
    manifest_container = Manifest.from_dict({
        "name": "pocket_predictor",
        "version": "1.0.0",
        "inputs": {"structure": "Structure"},
        "outputs": {"pockets": "Pocket"},
        "execution": {
            "type": "container",
            "image": "bio/pocket_finder:1.0",
            "command": dummy_cmd
        },
        "adapter": "protflow.adapters.parsers:parse_dummy_pocket_finder"
    })

    # 3. HTTP / WEB SERVICE MANIFEST
    manifest_http = Manifest.from_dict({
        "name": "pocket_predictor",
        "version": "1.0.0",
        "inputs": {"structure": "Structure"},
        "outputs": {"pockets": "Pocket"},
        "execution": {
            "type": "http",
            "url": f"{mock_server}/pockets/{{inputs.structure.protein_id}}",
            "method": "GET"
        },
        "adapter": "protflow.adapters.parsers:parse_dummy_pocket_finder",
        "is_remote": True
    })

    privacy_warnings_logged = []
    def on_privacy_warn(tool_name, endpoint, protein_id):
        privacy_warnings_logged.append((tool_name, endpoint, protein_id))

    runner = Runner(
        cache_dir=tmp_path / "cache",
        enable_cache=True,
        privacy_warning_handler=on_privacy_warn
    )

    # Run Local
    out_local = runner.run(manifest_local, inputs, output_dir=tmp_path / "out_local")
    assert "pockets" in out_local
    assert isinstance(out_local["pockets"], Pocket)
    assert out_local["pockets"].df["score"].iloc[0] == 0.88
    assert out_local["pockets"].provenance.is_remote is False

    # Run Container
    out_container = runner.run(manifest_container, inputs, output_dir=tmp_path / "out_container")
    assert "pockets" in out_container
    assert isinstance(out_container["pockets"], Pocket)
    assert out_container["pockets"].df["score"].iloc[0] == 0.88

    # Run HTTP
    out_http = runner.run(manifest_http, inputs, output_dir=tmp_path / "out_http")
    assert "pockets" in out_http
    assert isinstance(out_http["pockets"], Pocket)
    assert out_http["pockets"].df["score"].iloc[0] == 0.88
    assert out_http["pockets"].provenance.is_remote is True
    assert mock_server in out_http["pockets"].provenance.endpoint

    # Verify Privacy Warning was triggered for HTTP mode
    assert len(privacy_warnings_logged) == 1
    assert privacy_warnings_logged[0][0] == "pocket_predictor"
    assert privacy_warnings_logged[0][2] == structure.protein_id

    # Verify identical data shape across all 3
    assert out_local["pockets"].df.equals(out_container["pockets"].df)
    assert out_local["pockets"].df.equals(out_http["pockets"].df)

    # Verify Caching on a second run
    cached_out = runner.run(manifest_http, inputs)
    assert cached_out["pockets"].provenance.parameters.get("cache_hit") is True


def test_runner_cache_roundtrip_collections(tmp_path):
    from protflow.types import Sequence, SequenceCollection, Annotation, Alignment, Provenance
    runner = Runner(cache_dir=tmp_path / "cache", enable_cache=True)
    prov = Provenance(tool_name="test_tool", parameters={"k": "v"})

    seqs = SequenceCollection(protein_id="COL1", sequences={
        "P1": Sequence(protein_id="P1", sequence="MKFLV", description="test1"),
        "P2": Sequence(protein_id="P2", sequence="AGGLL", description="test2")
    })
    annot = Annotation(protein_id="COL1", records=[{"seqid": "P1", "type": "CDS", "start": 1, "end": 5}])
    align = Alignment(protein_id="COL1", target_id="P2", identity=0.8, score=100.0)

    outputs = {"seqs": seqs, "annot": annot, "align": align}
    cache_file = tmp_path / "cache" / "test_cache.json"
    runner._save_to_cache(cache_file, outputs, prov)

    assert cache_file.exists()
    loaded = runner._load_cached_outputs(cache_file, None)

    assert "seqs" in loaded and isinstance(loaded["seqs"], SequenceCollection)
    assert len(loaded["seqs"]) == 2
    assert loaded["seqs"]["P1"].sequence == "MKFLV"
    assert "annot" in loaded and isinstance(loaded["annot"], Annotation)
    assert len(loaded["annot"].records) == 1
    assert "align" in loaded and isinstance(loaded["align"], Alignment)
    assert loaded["align"].identity == 0.8
    assert loaded["seqs"].provenance.parameters.get("cache_hit") is True

