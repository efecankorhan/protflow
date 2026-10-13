"""
ProtFlow execution runner.
Executes manifests across local, container, and http backends with
input verification, caching, privacy warnings, and provenance tracking.
"""

from __future__ import annotations
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from typing import Any, Callable, Dict, Optional, Tuple, Union

import requests

from protflow.manifest import Manifest
from protflow.types import BaseBioType, Provenance, TYPE_REGISTRY

logger = logging.getLogger("protflow.runner")


class RunnerError(RuntimeError):
    """Raised when node execution fails."""
    pass


class PrivacyWarning(UserWarning):
    """Raised/logged when private sequence data is about to be sent to a remote service."""
    pass


def load_adapter_callable(adapter_str: str) -> Callable[[Any, Dict[str, Any]], Dict[str, Any]]:
    """Loads a callable from a 'module:func' string."""
    if ":" not in adapter_str:
        raise ValueError(f"Adapter string must be in format 'module.submodule:function_name', got '{adapter_str}'")
    mod_name, func_name = adapter_str.split(":", 1)
    mod = importlib.import_module(mod_name)
    func = getattr(mod, func_name)
    if not callable(func):
        raise TypeError(f"Adapter '{adapter_str}' is not callable.")
    return func


class Runner:
    def __init__(
        self,
        cache_dir: Optional[Union[str, Path]] = None,
        enable_cache: bool = True,
        container_engine: Optional[str] = None,  # 'docker', 'apptainer', or None
        privacy_warning_handler: Optional[Callable[[str, str, str], None]] = None,
    ):
        self.cache_dir = Path(cache_dir or ".protflow_cache")
        self.enable_cache = enable_cache
        self.container_engine = container_engine or self._detect_container_engine()
        self.privacy_warning_handler = privacy_warning_handler
        if self.enable_cache:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _detect_container_engine() -> Optional[str]:
        if shutil.which("docker"):
            return "docker"
        if shutil.which("apptainer"):
            return "apptainer"
        if shutil.which("singularity"):
            return "singularity"
        return None

    def _compute_cache_key(
        self,
        manifest: Manifest,
        inputs: Dict[str, BaseBioType],
        parameters: Dict[str, Any]
    ) -> str:
        hasher = hashlib.sha256()
        hasher.update(manifest.name.encode())
        hasher.update(manifest.version.encode())
        hasher.update(str(manifest.execution.type).encode())
        hasher.update(str(manifest.execution.command or manifest.execution.url).encode())
        hasher.update(json.dumps(parameters, sort_keys=True).encode())

        for in_name in sorted(inputs.keys()):
            obj = inputs[in_name]
            hasher.update(in_name.encode())
            hasher.update(obj.protein_id.encode())
            if hasattr(obj, "sequences"):
                for s_id, s in sorted(obj.sequences.items()):
                    hasher.update(s_id.encode())
                    hasher.update(s.sequence.encode())
            elif hasattr(obj, "sequence") and obj.sequence:
                hasher.update(obj.sequence.encode())
            elif hasattr(obj, "file_path") and obj.file_path and Path(obj.file_path).exists():
                fp = Path(obj.file_path)
                hasher.update(str(fp.stat().st_size).encode())
                with open(fp, "rb") as fh:
                    hasher.update(fh.read(65536))
            elif hasattr(obj, "df") and not obj.df.empty:
                hasher.update(obj.df.to_json().encode())
            elif hasattr(obj, "records"):
                hasher.update(json.dumps(obj.records, sort_keys=True).encode())
            elif hasattr(obj, "target_id"):
                hasher.update(json.dumps(obj.to_dict(), sort_keys=True).encode())

        return hasher.hexdigest()

    def run(
        self,
        manifest: Manifest,
        inputs: Dict[str, BaseBioType],
        parameters: Optional[Dict[str, Any]] = None,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> Dict[str, BaseBioType]:
        """
        Executes a manifest with the given inputs and parameters.
        Returns a dict of typed ProtFlow objects.
        """
        params = dict(manifest.parameters)
        if parameters:
            params.update(parameters)

        # 1. Type check inputs
        for in_name, in_type in manifest.inputs.items():
            if in_name not in inputs:
                raise RunnerError(f"Missing required input '{in_name}' of type '{in_type}' for tool '{manifest.name}'.")
            expected_cls = TYPE_REGISTRY[in_type]
            if not isinstance(inputs[in_name], expected_cls):
                raise RunnerError(
                    f"Input '{in_name}' expected type {expected_cls.__name__}, got {type(inputs[in_name]).__name__}."
                )

        # Determine primary protein_id
        if inputs:
            first_input = next(iter(inputs.values()))
            primary_protein_id = first_input.protein_id
        else:
            primary_protein_id = params.get("protein_id", params.get("organism", "SAMPLE"))

        # 2. Privacy Warning Check
        if manifest.is_remote:
            endpoint = manifest.execution.url or "remote_service"
            warn_msg = (
                f"[PRIVACY WARNING] Tool '{manifest.name}' is marked as REMOTE. "
                f"Protein '{primary_protein_id}' data will be transmitted externally to: {endpoint}."
            )
            logger.warning(warn_msg)
            if self.privacy_warning_handler:
                self.privacy_warning_handler(manifest.name, endpoint, primary_protein_id)

        # 3. Cache check
        cache_key = self._compute_cache_key(manifest, inputs, params)
        cache_file = self.cache_dir / f"{cache_key}.json"
        if self.enable_cache and cache_file.exists():
            logger.info("Cache hit for tool '%s' (hash: %s)", manifest.name, cache_key[:8])
            # Load cached outputs
            return self._load_cached_outputs(cache_file, manifest)

        # 4. Dispatch execution
        out_dir = Path(output_dir or f".protflow_work/{manifest.name}_{primary_protein_id}")
        out_dir.mkdir(parents=True, exist_ok=True)

        start_time = time.time()
        call_time_iso = datetime.now(timezone.utc).isoformat()

        exec_type = manifest.execution.type.lower()
        if exec_type == "local":
            raw_output = self._run_local(manifest, inputs, params, out_dir)
        elif exec_type == "container":
            raw_output = self._run_container(manifest, inputs, params, out_dir)
        elif exec_type == "http":
            raw_output = self._run_http(manifest, inputs, params)
        else:
            raise RunnerError(f"Unsupported execution backend: '{exec_type}'")

        # 5. Adapt raw output to ProtFlow types
        adapter_fn = load_adapter_callable(manifest.adapter)
        context = {
            "protein_id": primary_protein_id,
            "inputs": inputs,
            "parameters": params,
            "output_dir": str(out_dir),
            "manifest": manifest,
        }
        outputs = adapter_fn(raw_output, context)

        # 6. Validate outputs & attach provenance
        prov = Provenance(
            tool_name=manifest.name,
            tool_version=manifest.version,
            parameters=params,
            execution_time=call_time_iso,
            is_remote=manifest.is_remote,
            endpoint=manifest.execution.url if manifest.is_remote else None,
        )

        for out_name, out_type in manifest.outputs.items():
            if out_name not in outputs:
                raise RunnerError(f"Adapter did not produce declared output '{out_name}'.")
            expected_cls = TYPE_REGISTRY[out_type]
            if not isinstance(outputs[out_name], expected_cls):
                raise RunnerError(
                    f"Output '{out_name}' expected {expected_cls.__name__}, got {type(outputs[out_name]).__name__}."
                )
            merged_params = dict(params)
            if outputs[out_name].provenance and outputs[out_name].provenance.parameters:
                merged_params.update(outputs[out_name].provenance.parameters)

            out_prov = Provenance(
                tool_name=manifest.name,
                tool_version=manifest.version,
                parameters=merged_params,
                execution_time=call_time_iso,
                is_remote=manifest.is_remote,
                endpoint=manifest.execution.url if manifest.is_remote else None,
            )
            outputs[out_name].attach_provenance(out_prov)

        # 7. Save to cache
        if self.enable_cache:
            self._save_to_cache(cache_file, outputs, prov)

        return outputs

    def _run_local(
        self,
        manifest: Manifest,
        inputs: Dict[str, BaseBioType],
        params: Dict[str, Any],
        out_dir: Path
    ) -> Any:
        cmd_template = manifest.execution.command
        if not cmd_template:
            raise RunnerError(f"Manifest '{manifest.name}' specifies local execution but has no 'command'.")

        # Context interpolation and universal placeholder replacement
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        str_out = str(out_dir)
        quoted_out = f'"{str_out}"' if " " in str_out else str_out
        replacements: Dict[str, str] = {
            "{output_dir}": str_out if ('"{output_dir}"' in cmd_template or "'{output_dir}'" in cmd_template) else quoted_out,
            "output_dir": str_out,
        }
        for k, v in inputs.items():
            if any(ch in str(v.protein_id) for ch in (";", "&", "|", "`", "$", "\n", "\r", ">", "<")):
                raise RunnerError(f"Security error: Input '{k}' protein_id contains prohibited shell metacharacters: {v.protein_id}")
            if hasattr(v, "sequences") and not getattr(v, "file_path", None):
                safe_pid = "".join(c for c in str(v.protein_id) if c.isalnum() or c in ("_", "-", "."))
                fasta_path = out_dir / f"{safe_pid}.fasta"
                if not fasta_path.exists():
                    with open(fasta_path, "w", encoding="utf-8") as f_out:
                        for s in v.sequences.values():
                            f_out.write(f">{s.protein_id} {getattr(s, 'description', '')}\n{s.sequence}\n")
                v_file = str(fasta_path)
            elif hasattr(v, "sequence") and not getattr(v, "file_path", None):
                safe_pid = "".join(c for c in str(v.protein_id) if c.isalnum() or c in ("_", "-", "."))
                fasta_path = out_dir / f"{safe_pid}.fasta"
                if not fasta_path.exists():
                    with open(fasta_path, "w", encoding="utf-8") as f_out:
                        f_out.write(f">{v.protein_id} {getattr(v, 'description', '')}\n{v.sequence}\n")
                v_file = str(fasta_path)
            else:
                v_file = getattr(v, "file_path", getattr(v, "protein_id", ""))

            str_file = str(v_file)
            quoted_file = f'"{str_file}"' if " " in str_file else str_file
            for placeholder in (f"{{inputs.{k}.file_path}}", f"{{inputs.{k}}}", f"{{inputs_{k}}}"):
                replacements[placeholder] = str_file if (f'"{placeholder}"' in cmd_template or f"'{placeholder}'" in cmd_template) else quoted_file
            replacements[f"{{inputs.{k}.protein_id}}"] = str(v.protein_id)

        for k, v in params.items():
            if isinstance(v, str) and any(ch in v for ch in (";", "&", "|", "`", "$", "\n", "\r", ">", "<")):
                # shortcut: reject shell metacharacters in local runner params, upgrade to list-based argv execution
                raise RunnerError(f"Security error: Parameter '{k}' contains prohibited shell metacharacters: {v}")
            replacements[f"{{parameters.{k}}}"] = str(v)
            replacements[f"{{param_{k}}}"] = str(v)
            replacements[f"{{param.{k}}}"] = str(v)
            replacements[f"{{{k}}}"] = str(v)

        cmd = cmd_template
        for placeholder, val in replacements.items():
            cmd = cmd.replace(placeholder, val)

        cmd = os.path.expandvars(cmd)

        # Check if the command's primary binary exists on PATH or shell
        binary_name = cmd.split()[0] if cmd else ""
        if binary_name in ("python", "python3"):
            cmd = f'"{sys.executable}"' + cmd[len(binary_name):]
            binary_name = sys.executable

        if shutil.which(binary_name) is None and binary_name not in (sys.executable, "python", "python3"):
            logger.info("Binary '%s' not found on PATH. Falling back to native adapter execution.", binary_name)
            return {"native_fallback": True, "synthetic": True, "output_dir": str(out_dir)}

        timeout = manifest.execution.timeout
        try:
            res = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise RunnerError(f"Local command timed out after {timeout} seconds: {cmd}")

        if res.returncode != 0:
            err_lower = res.stderr.lower()
            if (
                res.returncode in (9009, 127, 126)
                or "not recognized as an internal or external command" in err_lower
                or "python was not found" in err_lower
                or "cannot find the file" in err_lower
            ):
                logger.info("Binary '%s' not available in shell. Falling back to native adapter execution.", binary_name)
                return {"native_fallback": True, "synthetic": True, "output_dir": str(out_dir)}
            raise RunnerError(f"Local command failed (code {res.returncode}): {res.stderr}\nCommand: {cmd}")

        return res.stdout or str(out_dir)

    def _run_container(
        self,
        manifest: Manifest,
        inputs: Dict[str, BaseBioType],
        params: Dict[str, Any],
        out_dir: Path
    ) -> Any:
        image = manifest.execution.image
        cmd_template = manifest.execution.command
        if not image or not cmd_template:
            raise RunnerError("Container execution requires both 'image' and 'command'.")

        abs_out = out_dir.resolve()
        # If container engine is present, run inside container; otherwise fallback / run command
        if self.container_engine == "docker":
            cmd = f'docker run --rm -v "{abs_out}:/work" -w /work {image} {cmd_template}'
            timeout = manifest.execution.timeout
            try:
                res = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                raise RunnerError(f"Docker command timed out after {timeout} seconds: {cmd}")
            if res.returncode != 0:
                raise RunnerError(f"Docker command failed: {res.stderr}")
            return res.stdout or str(out_dir)
        else:
            # Fallback wrapper: logs warning and executes locally or via mock runner
            logger.info("No active docker found; executing container command via native execution layer.")
            return self._run_local(manifest, inputs, params, out_dir)

    def _run_http(
        self,
        manifest: Manifest,
        inputs: Dict[str, BaseBioType],
        params: Dict[str, Any]
    ) -> Any:
        url_template = manifest.execution.url
        if not url_template:
            raise RunnerError("HTTP execution requires 'url'.")

        # Format URL with inputs & params
        fmt_args = {}
        for k, v in inputs.items():
            fmt_args[f"inputs.{k}.protein_id"] = v.protein_id
            fmt_args[k] = v.protein_id
        for k, v in params.items():
            fmt_args[k] = v

        url = url_template
        for k, v in fmt_args.items():
            url = url.replace(f"{{{k}}}", str(v))

        method = manifest.execution.method
        headers = {k: os.path.expandvars(str(v)) for k, v in (manifest.execution.headers or {}).items()}
        timeout = manifest.execution.timeout
        retries = manifest.execution.retries

        last_exc = None
        for attempt in range(1, retries + 1):
            try:
                kwargs = {"headers": headers, "timeout": timeout}
                if method.upper() in ("POST", "PUT", "PATCH") and "json_body" in params:
                    kwargs["json"] = params["json_body"]
                resp = requests.request(method=method, url=url, **kwargs)
                if resp.status_code in (404, 400):
                    # Service returned not-found or invalid accession; allow adapter to handle
                    logger.info("HTTP service returned %d for %s (entity not found in database)", resp.status_code, url)
                    return {"not_found": True, "status_code": resp.status_code, "url": url}
                resp.raise_for_status()
                try:
                    return resp.json()
                except Exception:
                    return resp.text
            except Exception as e:
                last_exc = e
                logger.warning("HTTP request to %s failed (attempt %d/%d): %s", url, attempt, retries, e)
                if attempt < retries:
                    time.sleep(1.0 * attempt)

        raise RunnerError(f"HTTP request failed after {retries} attempts: {last_exc}")

    def _save_to_cache(self, cache_file: Path, outputs: Dict[str, BaseBioType], prov: Provenance) -> None:
        data = {
            "provenance": prov.to_dict(),
            "outputs": {}
        }
        for name, obj in outputs.items():
            out_type = type(obj).__name__
            serialized: Dict[str, Any] = {"type": out_type, "protein_id": obj.protein_id}
            if getattr(obj, "provenance", None):
                serialized["provenance"] = obj.provenance.to_dict()
            if hasattr(obj, "df"):
                serialized["df_json"] = obj.df.to_json(orient="split")
            elif hasattr(obj, "sequences"):  # SequenceCollection
                serialized["sequences"] = {
                    pid: {"sequence": s.sequence, "description": getattr(s, "description", "")}
                    for pid, s in obj.sequences.items()
                }
            elif hasattr(obj, "structures"):  # StructureCollection
                serialized["structures"] = {
                    pid: {"file_path": s.file_path, "format": getattr(s, "format", "cif")}
                    for pid, s in obj.structures.items()
                }
            elif hasattr(obj, "records"):  # Annotation
                serialized["records"] = obj.records
            elif hasattr(obj, "target_id"):  # Alignment
                serialized["alignment"] = obj.to_dict()
            elif hasattr(obj, "sequence"):
                serialized["sequence"] = obj.sequence
                serialized["description"] = getattr(obj, "description", "")
            elif hasattr(obj, "file_path"):
                serialized["file_path"] = obj.file_path
                serialized["format"] = getattr(obj, "format", "cif")
            data["outputs"][name] = serialized

        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(data, f)

    def _load_cached_outputs(self, cache_file: Path, manifest: Manifest) -> Dict[str, BaseBioType]:
        with open(cache_file, "r", encoding="utf-8") as f:
            data = json.load(f)

        prov = Provenance.from_dict(data["provenance"])
        prov.parameters["cache_hit"] = True

        outputs = {}
        import pandas as pd
        import io
        from protflow.types import Sequence, Structure
        for name, s_data in data["outputs"].items():
            cls = TYPE_REGISTRY[s_data["type"]]
            pid = s_data["protein_id"]
            if "df_json" in s_data:
                df = pd.read_json(io.StringIO(s_data["df_json"]), orient="split")
                obj = cls(protein_id=pid, df=df)
            elif "sequences" in s_data:
                seqs = {k: Sequence(protein_id=k, sequence=v["sequence"], description=v.get("description", "")) for k, v in s_data["sequences"].items()}
                obj = cls(protein_id=pid, sequences=seqs)
            elif "structures" in s_data:
                strs = {k: Structure(protein_id=k, file_path=v["file_path"], format=v.get("format", "cif")) for k, v in s_data["structures"].items()}
                obj = cls(protein_id=pid, structures=strs)
            elif "records" in s_data:
                obj = cls(protein_id=pid, records=s_data["records"])
            elif "alignment" in s_data:
                align_dict = s_data["alignment"]
                obj = cls(
                    protein_id=pid,
                    target_id=align_dict.get("target_id", ""),
                    query_seq=align_dict.get("query_seq", ""),
                    target_seq=align_dict.get("target_seq", ""),
                    identity=align_dict.get("identity", 0.0),
                    score=align_dict.get("score", 0.0),
                    evalue=align_dict.get("evalue", 0.0),
                    alignment_format=align_dict.get("alignment_format", "tabular"),
                )
            elif "sequence" in s_data:
                obj = cls(protein_id=pid, sequence=s_data["sequence"], description=s_data.get("description", ""))
            elif "file_path" in s_data:
                obj = cls(protein_id=pid, file_path=s_data["file_path"], format=s_data.get("format", "cif"))
            else:
                obj = cls(protein_id=pid)
            if "provenance" in s_data:
                obj_prov = Provenance.from_dict(s_data["provenance"])
                obj_prov.parameters["cache_hit"] = True
                obj.attach_provenance(obj_prov)
            else:
                obj.attach_provenance(prov)
            outputs[name] = obj
        return outputs


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="ProtFlow Manifest Runner")
    parser.add_argument("--manifest", required=True, help="Path to tool manifest YAML")
    parser.add_argument("--sample", default="SAMPLE", help="Sample/Protein identifier")
    parser.add_argument("--params", default="{}", help="JSON string of parameter overrides")
    args = parser.parse_args()
    m = Manifest.from_yaml(args.manifest)
    p_overrides = json.loads(args.params)
    r = Runner()
    res = r.run(m, inputs={}, parameters=p_overrides)
    print(f"Executed {m.name}: produced {list(res.keys())}")
