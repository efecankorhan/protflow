"""
Manifest parser and validator for ProtFlow nodes.
Loads tool manifests from YAML files and validates type contracts and execution modes.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import yaml
from protflow.types import TYPE_REGISTRY


class ManifestValidationError(ValueError):
    """Raised when a manifest fails validation."""
    pass


@dataclass
class ExecutionConfig:
    type: str  # local, container, http
    command: Optional[str] = None
    image: Optional[str] = None
    url: Optional[str] = None
    method: str = "GET"
    headers: Dict[str, str] = field(default_factory=dict)
    timeout: int = 30
    retries: int = 3


@dataclass
class Manifest:
    name: str
    version: str
    inputs: Dict[str, str]
    outputs: Dict[str, str]
    execution: ExecutionConfig
    adapter: str
    description: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)
    is_remote: bool = False

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Manifest:
        # Validate required fields
        for req in ("name", "version", "inputs", "outputs", "execution", "adapter"):
            if req not in data:
                raise ManifestValidationError(f"Missing required field in manifest: '{req}'")

        # Validate types against TYPE_REGISTRY
        for in_name, in_type in data["inputs"].items():
            if in_type not in TYPE_REGISTRY:
                raise ManifestValidationError(
                    f"Unknown input type '{in_type}' for input '{in_name}'. "
                    f"Allowed types: {list(TYPE_REGISTRY.keys())}"
                )

        for out_name, out_type in data["outputs"].items():
            if out_type not in TYPE_REGISTRY:
                raise ManifestValidationError(
                    f"Unknown output type '{out_type}' for output '{out_name}'. "
                    f"Allowed types: {list(TYPE_REGISTRY.keys())}"
                )

        exec_data = data["execution"]
        exec_type = exec_data.get("type", "").lower()
        if exec_type not in ("local", "container", "http"):
            raise ManifestValidationError(
                f"Invalid execution type '{exec_type}'. Must be 'local', 'container', or 'http'."
            )

        is_remote = data.get("is_remote", False) or (exec_type == "http")

        exec_cfg = ExecutionConfig(
            type=exec_type,
            command=exec_data.get("command"),
            image=exec_data.get("image"),
            url=exec_data.get("url"),
            method=exec_data.get("method", "GET").upper(),
            headers=exec_data.get("headers", {}),
            timeout=exec_data.get("timeout", 30),
            retries=exec_data.get("retries", 3),
        )

        return cls(
            name=data["name"],
            version=str(data["version"]),
            inputs=data["inputs"],
            outputs=data["outputs"],
            execution=exec_cfg,
            adapter=data["adapter"],
            description=data.get("description", ""),
            parameters=data.get("parameters", {}),
            is_remote=is_remote,
        )

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> Manifest:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Manifest file not found: {path}")
        with open(p, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)

    def to_yaml(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "is_remote": self.is_remote,
            "inputs": self.inputs,
            "outputs": self.outputs,
            "parameters": self.parameters,
            "execution": {
                "type": self.execution.type,
                "command": self.execution.command,
                "image": self.execution.image,
                "url": self.execution.url,
                "method": self.execution.method,
                "headers": self.execution.headers,
                "timeout": self.execution.timeout,
                "retries": self.execution.retries,
            },
            "adapter": self.adapter,
        }
        with open(p, "w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, sort_keys=False)
        return p
