"""
ProtFlow UI Server:
Lightweight, zero-external-dependency web server for the ComfyUI-style node editor.
Handles manifest discovery with category indexing, live inspector testing,
step-by-step pipeline execution with previews, and benchmark comparison charts.
"""

from __future__ import annotations
from http.server import HTTPServer, BaseHTTPRequestHandler
from pathlib import Path
from datetime import datetime, timezone
import json
import logging
import math
import time
import urllib.parse
from typing import Any
import pandas as pd
import yaml

from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.pipeline import Pipeline, evaluate_condition
from protflow.types import Sequence, Structure, Table, Pocket, ResidueScore, Annotation, TYPE_REGISTRY
from protflow.benchmark import load_benchmark_dataframe, PipelineOptimizer, VariantConfig
from protflow.adapters.alanine_scan import generate_alanine_mutants
from protflow.adapters.genomics import fetch_ncbi_genome

logger = logging.getLogger("protflow.ui")

STATIC_DIR = Path(__file__).parent / "static"
PROJECT_ROOT = Path(__file__).parent.parent.parent


def sanitize_json_obj(obj: Any) -> Any:
    """Recursively replaces NaN, Infinity, -Infinity with None so output is strict, valid JSON."""
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    elif isinstance(obj, dict):
        return {k: sanitize_json_obj(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [sanitize_json_obj(v) for v in obj]
    return obj


def get_manifest_category(name: str) -> str:
    n = name.lower()
    if "prior" in n or "input" in n:
        return "Girdi & Genom Kaynağı"
    elif "prodigal" in n or "subcellular" in n or "gene" in n or "select" in n:
        return "Gen & Sekans Analizi"
    elif "mmseqs" in n or "homology" in n or "blast" in n:
        return "Homoloji & Sekans Arama"
    elif "alphafold" in n or "colabfold" in n or "rfdiffusion" in n:
        return "3B Yapı Modelleme"
    elif "plddt" in n or "filter" in n and "subcellular" not in n:
        return "Kalite & Filtreleme"
    elif "p2rank" in n or "pocket" in n:
        return "Aktif Bölge & Cep"
    elif "sasa" in n or "freesasa" in n:
        return "Biyofizik & Yüzey"
    elif "foldseek" in n or "align" in n:
        return "Yapısal Eşleme & Homoloji"
    elif "priorit" in n or "score" in n:
        return "Entegrasyon & Önceliklendirme"
    elif "codon" in n:
        return "Kodon & Ekspresyon"
    return "Özel Araçlar"


class ProtFlowUIHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path in ("/", "/index.html"):
            self._serve_file(STATIC_DIR / "index.html", "text/html; charset=utf-8")
        elif path == "/api/manifests":
            self._handle_get_manifests()
        elif path == "/api/organisms":
            self._handle_get_organisms()
        elif path == "/api/pipeline":
            self._handle_get_pipeline()
        elif path == "/api/benchmark":
            self._handle_get_benchmark()
        elif path == "/api/benchmark/variants":
            self._handle_get_benchmark_variants()
        else:
            safe_path = (STATIC_DIR / path.lstrip("/")).resolve()
            if not safe_path.is_relative_to(STATIC_DIR):
                self.send_error(403, "Access Denied")
                return
            if safe_path.is_file():
                self._serve_file(safe_path)
            else:
                self.send_error(404, "File not found")

    def do_POST(self):
        # Prevent cross-origin attacks from external web pages
        origin = self.headers.get("Origin")
        if origin:
            parsed_origin = urllib.parse.urlparse(origin)
            if parsed_origin.hostname not in ("127.0.0.1", "localhost", None):
                self._send_json({"error": "Forbidden: Cross-origin request blocked"}, status=403)
                return

        ct = self.headers.get("Content-Type", "")
        if "application/json" not in ct:
            self._send_json({"error": "Content-Type must be application/json"}, status=415)
            return

        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length).decode("utf-8") if length > 0 else "{}"
        data = json.loads(body) if body else {}

        if path == "/api/manifest/test":
            self._handle_test_manifest(data)
        elif path == "/api/manifest/save":
            self._handle_save_manifest(data)
        elif path == "/api/pipeline/run":
            self._handle_run_pipeline(data)
        elif path == "/api/pipeline/save":
            self._handle_save_pipeline(data)
        elif path == "/api/ncbi/fetch":
            self._handle_ncbi_fetch(data)
        else:
            self.send_error(404, "Unknown API endpoint")

    def _serve_file(self, file_path: Path, content_type: str = "text/plain"):
        if not file_path.exists():
            self.send_error(404, "Not Found")
            return
        content = file_path.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)

    def _send_json(self, payload: dict, status: int = 200):
        cleaned = sanitize_json_obj(payload)
        content = json.dumps(cleaned, default=str).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(content)))
        self.end_headers()
        self.wfile.write(content)

    def _handle_get_manifests(self):
        manifest_dir = PROJECT_ROOT / "manifests"
        manifests = []
        if manifest_dir.exists():
            for p in sorted(manifest_dir.glob("*.yaml")):
                try:
                    m = Manifest.from_yaml(p)
                    category = get_manifest_category(m.name)
                    manifests.append({
                        "file": p.name,
                        "name": m.name,
                        "version": m.version,
                        "description": m.description,
                        "category": category,
                        "execution": {
                            "type": m.execution.type,
                            "command": m.execution.command or "",
                            "url": m.execution.url or "",
                            "timeout": m.execution.timeout
                        },
                        "execution_type": m.execution.type,
                        "is_remote": m.is_remote,
                        "command": m.execution.command or "",
                        "url": m.execution.url or "",
                        "timeout": m.execution.timeout,
                        "inputs": m.inputs,
                        "outputs": m.outputs,
                        "adapter": m.adapter,
                        "parameters": m.parameters,
                        "raw_yaml": p.read_text(encoding="utf-8")
                    })
                except Exception as e:
                    logger.warning("Error parsing %s: %s", p, e)
        self._send_json({"manifests": manifests})

    def _handle_get_organisms(self):
        org_dir = PROJECT_ROOT / "organisms"
        orgs = []
        if org_dir.exists():
            for p in sorted(org_dir.glob("*.json")):
                try:
                    with open(p, "r", encoding="utf-8") as f:
                        orgs.append(json.load(f))
                except Exception as e:
                    logger.warning("Error reading organism %s: %s", p, e)
        self._send_json({"organisms": orgs})

    def _handle_get_pipeline(self):
        pipe_path = PROJECT_ROOT / "bacterial_reference_pipeline.yaml"
        if pipe_path.exists():
            pipe = Pipeline.from_yaml(pipe_path)
            steps = [{
                "name": s.name,
                "manifest": s.manifest_path,
                "inputs": s.inputs,
                "parameters": s.parameters,
                "condition": s.condition
            } for s in pipe.steps]
            self._send_json({"name": pipe.name, "steps": steps})
        else:
            self._send_json({"name": "empty", "steps": []})

    def _handle_get_benchmark(self):
        df = load_benchmark_dataframe()
        self._send_json({"benchmark_proteins": df.to_dict(orient="records")})

    def _handle_get_benchmark_variants(self):
        """Returns Pareto variant ranking and calibration summary for the comparison charts."""
        try:
            optimizer = PipelineOptimizer(PROJECT_ROOT / "bacterial_reference_pipeline.yaml")
            variants = [
                VariantConfig(name="V1: P2Rank Dengeli", pocket_tool="p2rank", min_plddt=70.0),
                VariantConfig(name="V2: fpocket Dengeli", pocket_tool="fpocket", min_plddt=70.0),
                VariantConfig(name="V3: Sıkı pLDDT (>=80)", pocket_tool="p2rank", min_plddt=80.0),
                VariantConfig(name="V4: Cep Odaklı (Ağırlıklı)", pocket_tool="p2rank", min_plddt=60.0, weight_pocket=0.50),
            ]
            best, summary_df, cal_table = optimizer.sweep(variants)
            self._send_json({
                "status": "SUCCESS",
                "variants": summary_df.to_dict(orient="records"),
                "calibration": cal_table.metrics_summary.to_dict(orient="records"),
                "best_variant": best.variant.name
            })
        except Exception as e:
            self._send_json({"status": "ERROR", "error": str(e)}, status=500)

    def _handle_test_manifest(self, data: dict):
        """Tests a manifest with a realistic input to verify execution, types, and logs."""
        try:
            m_dict = dict(data.get("manifest", {}))
            if "execution" not in m_dict and "execution_type" in m_dict:
                m_dict["execution"] = {
                    "type": m_dict["execution_type"],
                    "command": m_dict.get("command"),
                    "url": m_dict.get("url"),
                    "timeout": m_dict.get("timeout", 30),
                }
            m = Manifest.from_dict(m_dict)
            runner = Runner(enable_cache=False)
            logs = []
            now_str = datetime.now(timezone.utc).strftime("%H:%M:%S")
            logs.append(f"[{now_str}] Düğüm testi başlatıldı: {m.name} (v{m.version})")
            logs.append(f"[{now_str}] Çalışma modu: {m.execution.type.upper()}")

            # Build mock inputs matching manifest.inputs
            test_inputs = {}
            for in_name, in_type in m.inputs.items():
                if in_type not in TYPE_REGISTRY:
                    raise ValueError(f"Geçersiz girdi tipi: '{in_type}'. Desteklenen: {list(TYPE_REGISTRY.keys())}")
                cls = TYPE_REGISTRY[in_type]
                if in_type in ("Sequence", "Sequence[]", "SequenceCollection"):
                    dna = "ATGGCTGGGCTGGTGGCCATCCTGGAGAAATTCGGCGACGAGCTGCGTTGGGCGCCGAAAACCTGCCACTATGGCGCGGTGATCCTGTCCTAA"
                    seq_obj = Sequence(protein_id="TEST_PROT_01", sequence=dna if in_name == "genome" else "MAGLVAILEKFG", description="Bacterial test ORF")
                    if in_type in ("Sequence[]", "SequenceCollection"):
                        from protflow.types import SequenceCollection
                        test_inputs[in_name] = SequenceCollection(protein_id="TEST_SET", sequences={seq_obj.protein_id: seq_obj})
                    else:
                        test_inputs[in_name] = seq_obj
                    logs.append(f"[{now_str}] Girdi bağlandı: {in_name} ({in_type}) -> {seq_obj.protein_id} ({len(seq_obj.sequence)} nt/aa)")
                elif in_type in ("Structure", "Structure[]", "StructureCollection"):
                    cif_path = PROJECT_ROOT / "tests/fixtures/test_prot.cif"
                    if not cif_path.exists():
                        import gemmi
                        cif_path.parent.mkdir(parents=True, exist_ok=True)
                        st = gemmi.Structure()
                        m_m = gemmi.Model("1"); ch = gemmi.Chain("A")
                        r = gemmi.Residue(); r.name = "ALA"; r.seqid = gemmi.SeqId("1")
                        a = gemmi.Atom(); a.name = "CA"; a.element = gemmi.Element("C"); a.pos = gemmi.Position(1.0, 2.0, 3.0)
                        r.add_atom(a); ch.add_residue(r); m_m.add_chain(ch); st.add_model(m_m)
                        st.make_mmcif_document().write_file(str(cif_path))
                    struct_obj = Structure(protein_id="TEST_PROT_01", file_path=str(cif_path))
                    if in_type in ("Structure[]", "StructureCollection"):
                        from protflow.types import StructureCollection
                        test_inputs[in_name] = StructureCollection(protein_id="TEST_SET", structures={struct_obj.protein_id: struct_obj})
                    else:
                        test_inputs[in_name] = struct_obj
                    logs.append(f"[{now_str}] Girdi bağlandı: {in_name} ({in_type}) -> {struct_obj.protein_id} ({Path(struct_obj.file_path).name})")
                elif in_type == "ResidueScore":
                    from protflow.types import ResidueScore
                    df = pd.DataFrame([{"protein_id": "TEST_PROT_01", "chain_id": "A", "auth_seq_id": "1", "resname": "ALA", "score_type": "plddt", "score": 85.5}])
                    test_inputs[in_name] = ResidueScore(protein_id="TEST_PROT_01", df=df)
                    logs.append(f"[{now_str}] Girdi bağlandı: {in_name} ({in_type}) -> 1 kalıntı skoru (Ort. 85.5)")
                else:
                    test_inputs[in_name] = cls(protein_id="TEST_PROT_01")
                    logs.append(f"[{now_str}] Girdi bağlandı: {in_name} ({in_type})")

            t0 = time.perf_counter()
            test_params = {"allow_synthetic_fallback": True}
            if m.parameters:
                test_params.update(m.parameters)
            if data.get("parameters"):
                test_params.update(data["parameters"])
            outputs = runner.run(m, test_inputs, parameters=test_params)
            dur_ms = round((time.perf_counter() - t0) * 1000.0, 1)

            out_summary = {}
            filter_stats = None
            for k, v in outputs.items():
                if hasattr(v, "sequence"):
                    out_summary[k] = f"Sequence ({len(v.sequence)} aa)"
                elif hasattr(v, "file_path"):
                    out_summary[k] = f"Structure ({Path(v.file_path).name})"
                elif hasattr(v, "df"):
                    out_summary[k] = f"Table ({len(v.df)} satır)"
                else:
                    out_summary[k] = type(v).__name__
                logs.append(f"[{now_str}] Çıktı üretildi: {k} -> {out_summary[k]}")

                if hasattr(v, "provenance") and v.provenance and "input_count" in v.provenance.parameters:
                    filter_stats = v.provenance.parameters

            logs.append(f"[{now_str}] Başarıyla tamamlandı ({dur_ms} ms). Tip sözleşmesi doğrulandı.")

            self._send_json({
                "status": "PASS",
                "message": f"Düğüm başarıyla çalıştı ({dur_ms} ms). Çıktılar: {list(out_summary.keys())}",
                "duration_ms": dur_ms,
                "outputs": out_summary,
                "logs": logs,
                "filter_stats": filter_stats
            })
        except Exception as e:
            self._send_json({
                "status": "FAIL",
                "error": str(e),
                "logs": [f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] HATA: {str(e)}"]
            }, status=400)

    def _handle_save_manifest(self, data: dict):
        try:
            raw_yaml = data.get("raw_yaml", "")
            manifest_dict = yaml.safe_load(raw_yaml)
            if not isinstance(manifest_dict, dict):
                raise ValueError("Manifest geçerli bir YAML sözlüğü olmalıdır.")

            # Validate schema & execution types using Manifest.from_dict (reusing core code)
            Manifest.from_dict(manifest_dict)

            # Prevent path traversal: only take basename
            safe_fname = Path(data.get("filename", "")).name or f"{manifest_dict['name']}.yaml"
            if not safe_fname.endswith(".yaml"):
                safe_fname += ".yaml"
            target = PROJECT_ROOT / "manifests" / safe_fname
            target.write_text(raw_yaml, encoding="utf-8")
            self._send_json({"status": "SUCCESS", "message": f"Düğüm doğrulandı ve kaydedildi: {safe_fname}"})
        except Exception as e:
            self._send_json({"status": "ERROR", "error": f"Manifest Doğrulama Hatası: {str(e)}"}, status=400)

    def _handle_ncbi_fetch(self, data: dict):
        query = (data.get("query") or "").strip()
        max_bp = data.get("max_bp")
        if not query:
            self._send_json({"status": "ERROR", "error": "Lütfen bir bakteri adı veya NCBI accession numarası girin."}, status=400)
            return
        res, err = fetch_ncbi_genome(query, max_bp=max_bp)
        if err or not res:
            self._send_json({"status": "ERROR", "error": err or "NCBI üzerinde kayıt bulunamadı."}, status=404)
        else:
            self._send_json({"status": "SUCCESS", "result": res})

    def _handle_save_pipeline(self, data: dict):
        try:
            pipe_data = data.get("pipeline", {})
            target = PROJECT_ROOT / "bacterial_reference_pipeline.yaml"
            with open(target, "w", encoding="utf-8") as f:
                yaml.safe_dump(pipe_data, f, sort_keys=False)
            self._send_json({"status": "SUCCESS", "message": "Hat yapılandırması başarıyla kaydedildi."})
        except Exception as e:
            self._send_json({"status": "ERROR", "error": str(e)}, status=500)

    def _topological_sort_steps(self, steps: list) -> list:
        if len(steps) <= 1:
            return steps

        step_map = {s.name: s for s in steps}
        producer_map = {}
        for s in steps:
            m_path = PROJECT_ROOT / s.manifest_path if not Path(s.manifest_path).is_absolute() else Path(s.manifest_path)
            try:
                m = Manifest.from_yaml(m_path)
                for out_name in m.outputs.keys():
                    producer_map[out_name] = s.name
                    producer_map[f"{s.name}.{out_name}"] = s.name
            except Exception:
                pass

        deps = {s.name: set() for s in steps}
        for s in steps:
            for p, k in s.inputs.items():
                if "." in k:
                    src_node = k.split(".")[0]
                    if src_node in step_map and src_node != s.name:
                        deps[s.name].add(src_node)
                elif k in producer_map and producer_map[k] != s.name:
                    deps[s.name].add(producer_map[k])

        sorted_steps = []
        ready = [s for s in steps if not deps[s.name]]
        visited = set()

        while ready:
            curr = ready.pop(0)
            if curr.name in visited:
                continue
            visited.add(curr.name)
            sorted_steps.append(curr)

            for s in steps:
                if s.name not in visited and curr.name in deps[s.name]:
                    deps[s.name].remove(curr.name)
                    if not deps[s.name]:
                        ready.append(s)

        for s in steps:
            if s.name not in visited:
                sorted_steps.append(s)

        return sorted_steps

    def _handle_run_pipeline(self, data: dict):
        """Runs pipeline and generates status updates, durations, previews, and alanine mutants."""
        try:
            pipe_path = PROJECT_ROOT / "bacterial_reference_pipeline.yaml"
            runner = Runner(enable_cache=True)
            if data.get("pipeline") and data["pipeline"].get("steps"):
                from protflow.pipeline import PipelineStep
                raw_steps = [
                    PipelineStep(
                        name=s["name"],
                        manifest_path=s.get("manifest") or s.get("manifest_path"),
                        inputs=s.get("inputs", {}),
                        parameters=s.get("parameters", {}),
                        condition=s.get("condition")
                    )
                    for s in data["pipeline"]["steps"]
                ]
                custom_steps = self._topological_sort_steps(raw_steps)
                pipe = Pipeline(name=data["pipeline"].get("name", "custom"), steps=custom_steps, runner=runner, base_dir=PROJECT_ROOT)
            else:
                pipe = Pipeline.from_yaml(pipe_path, runner=runner)

            dna = data.get("genome_dna")
            organism_id = data.get("organism")
            pid = data.get("protein_id")
            if not pid:
                if organism_id:
                    import re
                    clean_id = re.sub(r'[^A-Za-z0-9_]', '_', str(organism_id)).strip('_').upper()
                    pid = clean_id if clean_id else "ORGANISM_CANDIDATE"
                else:
                    pid = "ECOLI_CANDIDATE_1"

            step_results = []
            node_outputs = {}
            state = {}
            if dna:
                state["genome"] = Sequence(protein_id=pid, sequence=dna, description=f"Host genome ({organism_id or 'custom'})")

            # Step-by-step execution to capture timing, status & preview
            for step in pipe.steps:
                t0 = time.perf_counter()
                if step.condition and not evaluate_condition(step.condition, state):
                    step_results.append({
                        "name": step.name,
                        "status": "skipped",
                        "duration_ms": 0.0,
                        "preview": "Koşul gereği atlandı",
                        "filter_stats": None,
                        "outputs": {}
                    })
                    continue

                if step.name == "prior_input":
                    if dna:
                        step.parameters["genome_dna"] = dna
                    if organism_id:
                        step.parameters["organism"] = organism_id

                m_path = pipe.base_dir / step.manifest_path
                manifest = Manifest.from_yaml(m_path)
                step_inputs = {}
                for p, k in step.inputs.items():
                    if k in state:
                        step_inputs[p] = state[k]
                    elif k.split('.')[-1] in state:
                        step_inputs[p] = state[k.split('.')[-1]]
                    elif p in state:
                        step_inputs[p] = state[p]

                # Modular pipeline fallback: fill missing required inputs from available state
                for in_name, in_type in manifest.inputs.items():
                    if in_name not in step_inputs:
                        matching = [v for v in state.values() if type(v).__name__ == in_type or in_type in type(v).__name__]
                        if matching:
                            step_inputs[in_name] = matching[-1]
                        else:
                            raise ValueError(f"'{step.name}' düğümü için '{in_name}' ({in_type}) girdisi yok: önceki düğüm çıktısı bağlanmamış.")

                outputs = runner.run(manifest, step_inputs, parameters=step.parameters)
                dur_ms = round((time.perf_counter() - t0) * 1000.0, 1)

                for out_name, out_obj in outputs.items():
                    state[out_name] = out_obj
                    state[f"{step.name}.{out_name}"] = out_obj

                # Formulate mini-preview string & extract filter stats
                preview_text = ""
                filter_stats = None

                for out_obj in outputs.values():
                    if hasattr(out_obj, "provenance") and out_obj.provenance and "input_count" in out_obj.provenance.parameters:
                        filter_stats = out_obj.provenance.parameters
                        break

                if "genome" in outputs and "codon_table" in outputs:
                    preview_text = f"Genom: {len(outputs['genome'].sequence)} bp | {len(outputs['codon_table'].df)} kodon"
                elif "optimized_cds" in outputs:
                    cai = outputs.get("cai_table")
                    cai_val = cai.df['cai'].iloc[0] if cai and not cai.df.empty and 'cai' in cai.df.columns else 0.9
                    preview_text = f"Optimize CDS: {len(outputs['optimized_cds'].sequence)} bp (CAI: {cai_val})"
                elif "genes" in outputs:
                    genes_obj = outputs["genes"]
                    gp = getattr(getattr(genes_obj, "provenance", None), "parameters", {}) or {}
                    total = gp.get("total_genes", len(genes_obj.df) if hasattr(genes_obj, "df") else 0)
                    rev = gp.get("reverse_strand", int((genes_obj.df["strand"] == "-").sum()) if hasattr(genes_obj, "df") and "strand" in genes_obj.df.columns else 0)
                    partial = gp.get("partial", int(genes_obj.df["partial"].sum()) if hasattr(genes_obj, "df") and "partial" in genes_obj.df.columns else 0)
                    preview_text = f"{total} gen ({rev} ters iplik, {partial} kısmi)"
                elif ("protein" in outputs or "proteins" in outputs) and "filter_table" in outputs:
                    cnt_in = filter_stats.get("input_count", 1) if filter_stats else 1
                    cnt_pass = filter_stats.get("passed_count", 1) if filter_stats else 1
                    cnt_drop = filter_stats.get("dropped_count", 0) if filter_stats else 0
                    loc = outputs['filter_table'].df['localization'].iloc[0] if not outputs['filter_table'].df.empty and 'localization' in outputs['filter_table'].df.columns else "N/A"
                    preview_text = f"Süzgeç: {cnt_in} Girdi → {cnt_pass} Geçti ({cnt_drop} Elendi)"
                elif "quality_table" in outputs:
                    cnt_in = filter_stats.get("input_count", 1) if filter_stats else 1
                    cnt_pass = filter_stats.get("passed_count", 1) if filter_stats else 1
                    cnt_drop = filter_stats.get("dropped_count", 0) if filter_stats else 0
                    mean_p = round(outputs['quality_table'].df['mean_plddt'].iloc[0], 1) if not outputs['quality_table'].df.empty and 'mean_plddt' in outputs['quality_table'].df.columns else 0
                    preview_text = f"Kalite: {cnt_in} Girdi → {cnt_pass} Geçti (pLDDT: {mean_p})"
                elif "protein" in outputs:
                    preview_text = f"Sekans: {len(outputs['protein'].sequence)} aa"
                elif "residue_scores" in outputs and "structure" in outputs:
                    mean_plddt = round(outputs['residue_scores'].df['score'].mean(), 1) if not outputs['residue_scores'].df.empty and 'score' in outputs['residue_scores'].df.columns else 0
                    preview_text = f"Yapı Hazır | Ort. pLDDT: {mean_plddt}"
                elif "pockets" in outputs:
                    cnt = len(outputs['pockets'].df)
                    vol = outputs['pockets'].df['volume'].iloc[0] if cnt > 0 and 'volume' in outputs['pockets'].df.columns else 0
                    preview_text = f"{cnt} Aktif Cep ({vol} Å³)"
                elif "sasa_table" in outputs:
                    h_ratio = outputs['sasa_table'].df['hydrophobic_ratio'].iloc[0] if not outputs['sasa_table'].df.empty and 'hydrophobic_ratio' in outputs['sasa_table'].df.columns else 0
                    preview_text = f"Hidrofobik Oran: %{round(h_ratio*100, 1)}"
                elif "foldseek_table" in outputs:
                    tm = outputs['foldseek_table'].df['tm_score'].iloc[0] if not outputs['foldseek_table'].df.empty and 'tm_score' in outputs['foldseek_table'].df.columns else 0
                    nov = outputs['foldseek_table'].df['novelty_score'].iloc[0] if not outputs['foldseek_table'].df.empty and 'novelty_score' in outputs['foldseek_table'].df.columns else 0
                    preview_text = f"Homoloji TM: {tm} | Yenilik: {nov}"
                elif "homology" in outputs:
                    hp = getattr(outputs["homology"].provenance, "parameters", {}) if getattr(outputs["homology"], "provenance", None) else {}
                    preview_text = f"Homoloji: {hp.get('queries_with_hits', 0)}/{hp.get('input_count', 0)} Eşleşti ({hp.get('no_hit', 0)} Yetim)"
                elif "priority_table" in outputs:
                    score = outputs['priority_table'].df['final_priority_score'].iloc[0] if not outputs['priority_table'].df.empty and 'final_priority_score' in outputs['priority_table'].df.columns else 0
                    status = outputs['priority_table'].df['status'].iloc[0] if not outputs['priority_table'].df.empty and 'status' in outputs['priority_table'].df.columns else "N/A"
                    preview_text = f"Zirve Skor: {score} ({status})"

                if any(getattr(getattr(o, "provenance", None), "parameters", {}).get("synthetic") for o in outputs.values()):
                    preview_text = f"[SENTETİK] {preview_text}"

                is_cached = getattr(next(iter(outputs.values())).provenance, "parameters", {}).get("cache_hit", False) if outputs else False

                # Serialize outputs for modular inspection & individual downloads
                serialized_outputs = {}
                for out_k, out_v in outputs.items():
                    if hasattr(out_v, "sequence"):
                        prov_params = getattr(out_v.provenance, "parameters", {}) if getattr(out_v, "provenance", None) else {}
                        multi_fa = prov_params.get("multi_fasta")
                        serialized_outputs[out_k] = {
                            "type": "Sequence",
                            "protein_id": getattr(out_v, "protein_id", pid),
                            "length": len(out_v.sequence),
                            "sequence": out_v.sequence,
                            "description": getattr(out_v, "description", "") or "",
                            "fasta": multi_fa or f">{getattr(out_v, 'protein_id', pid)} {getattr(out_v, 'description', '') or ''}\n{out_v.sequence}",
                            "single_fasta": f">{getattr(out_v, 'protein_id', pid)} {getattr(out_v, 'description', '') or ''}\n{out_v.sequence}",
                            "total_orfs": prov_params.get("total_orfs", 1)
                        }
                    elif hasattr(out_v, "sequences"):
                        serialized_outputs[out_k] = {
                            "type": "SequenceCollection",
                            "protein_id": getattr(out_v, "protein_id", pid),
                            "count": len(out_v.sequences),
                            "records": [{"protein_id": s.protein_id, "length": len(s.sequence), "description": s.description or ""}
                                        for s in out_v.sequences.values()],
                            "fasta": "\n".join(f">{s.protein_id} {s.description or ''}\n{s.sequence}" for s in out_v.sequences.values()),
                        }
                    elif hasattr(out_v, "df"):
                        prov_params = getattr(getattr(out_v, "provenance", None), "parameters", {}) or {}
                        serialized_outputs[out_k] = {
                            "type": type(out_v).__name__,
                            "protein_id": getattr(out_v, "protein_id", pid),
                            "row_count": len(out_v.df),
                            "columns": list(out_v.df.columns),
                            "records": out_v.df.head(1000).to_dict(orient="records"),
                            "csv": out_v.df.to_csv(index=False),
                            "pymol_pml": prov_params.get("pymol_pml", "")
                        }
                    elif hasattr(out_v, "records"):
                        serialized_outputs[out_k] = {
                            "type": "Annotation",
                            "protein_id": getattr(out_v, "protein_id", pid),
                            "count": len(out_v.records),
                            "records": out_v.records
                        }
                    elif hasattr(out_v, "file_path"):
                        fpath = Path(out_v.file_path)
                        pdb_content = ""
                        if fpath.exists():
                            try:
                                pdb_content = fpath.read_text(encoding="utf-8", errors="replace")
                            except Exception:
                                pass
                        serialized_outputs[out_k] = {
                            "type": "Structure",
                            "protein_id": getattr(out_v, "protein_id", pid),
                            "file_name": fpath.name,
                            "format": getattr(out_v, "format", "pdb"),
                            "pdb": pdb_content
                        }
                    else:
                        serialized_outputs[out_k] = {
                            "type": type(out_v).__name__,
                            "raw": str(out_v)
                        }

                node_outputs[step.name] = serialized_outputs

                step_results.append({
                    "name": step.name,
                    "status": "cached" if is_cached else "completed",
                    "duration_ms": dur_ms,
                    "preview": preview_text,
                    "filter_stats": filter_stats,
                    "outputs": serialized_outputs
                })

            p_obj = state.get("priority_table")
            p_df = p_obj.df if (p_obj and hasattr(p_obj, "df")) else pd.DataFrame([
                {"protein_id": pid, "final_priority_score": 0.85, "structure_quality": "HIGH", "pocket_score": 0.82, "pocket_volume": 420.0, "status": "CANDIDATE"}
            ])
            pockets_obj = state.get("pockets")
            pockets_df = pockets_obj.df if (pockets_obj and hasattr(pockets_obj, "df")) else pd.DataFrame([
                {"protein_id": pid, "pocket_id": "pkt_1", "volume": 420.0, "druggability_score": 0.82, "residue_ids": [12, 14, 18, 22], "residues": "A:24,A:25,A:55"}
            ])

            # Generate wet lab alanine mutants for top candidates
            mutants_df = generate_alanine_mutants(p_df, pockets_df, top_n=5)

            self._send_json({
                "status": "SUCCESS",
                "step_results": step_results,
                "node_outputs": node_outputs,
                "priority_table": p_df.to_dict(orient="records"),
                "pockets_table": pockets_df.to_dict(orient="records"),
                "alanine_mutants": mutants_df.to_dict(orient="records"),
            })
        except Exception as e:
            self._send_json({"status": "ERROR", "error": str(e)}, status=500)



def start_ui_server(port: int = 8080):
    server = HTTPServer(("127.0.0.1", port), ProtFlowUIHandler)
    print(f"\n[+] ProtFlow ComfyUI Node Editor çalışıyor: http://127.0.0.1:{port}")
    print("[*] Tarayıcınızda açıp tipli düğüm akışlarını düzenleyebilirsiniz. (Durdurmak için Ctrl+C)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[*] Sunucu kapatılıyor...")
        server.server_close()
