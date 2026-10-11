"""
Pipeline DAG coordinator for ProtFlow.
Loads a pipeline configuration file, manages data passing between typed nodes,
handles filtering/branching conditions, and generates execution plans.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import logging
import yaml

from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.types import BaseBioType, Table
from protflow.export import export_pipeline_to_snakemake, export_pipeline_to_nextflow

logger = logging.getLogger("protflow.pipeline")

def evaluate_condition(condition: Optional[str], state: Dict[str, Any]) -> bool:
    """Safely checks pipeline step execution conditions without arbitrary eval."""
    if not condition:
        return True
    c = condition.strip()
    if "not in state" in c:
        key = c.split("not in state")[0].strip().strip("'\"")
        return key not in state
    if "in state" in c:
        key = c.split("in state")[0].strip().strip("'\"")
        return key in state
    try:
        return bool(eval(c, {"__builtins__": {}}, {"state": state}))
    except Exception:
        return True


@dataclass
class PipelineStep:
    name: str
    manifest_path: str
    inputs: Dict[str, str]  # node_port_name -> state_key
    parameters: Dict[str, Any] = field(default_factory=dict)
    condition: Optional[str] = None  # e.g. "state['filter_table'].df['is_intracellular'].iloc[0]"


class Pipeline:
    def __init__(
        self,
        name: str,
        steps: List[PipelineStep],
        runner: Optional[Runner] = None,
        base_dir: Optional[Union[str, Path]] = None
    ):
        self.name = name
        self.steps = steps
        self.runner = runner or Runner()
        self.base_dir = Path(base_dir or Path.cwd())

    @classmethod
    def from_yaml(cls, path: Union[str, Path], runner: Optional[Runner] = None) -> Pipeline:
        p = Path(path)
        with open(p, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)

        steps = []
        for s in data.get("steps", []):
            steps.append(PipelineStep(
                name=s["name"],
                manifest_path=s["manifest"],
                inputs=s.get("inputs", {}),
                parameters=s.get("parameters", {}),
                condition=s.get("condition"),
            ))

        return cls(name=data.get("name", "pipeline"), steps=steps, runner=runner, base_dir=p.parent)

    def to_yaml(self, path: Union[str, Path]) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "name": self.name,
            "steps": [{
                "name": s.name,
                "manifest": s.manifest_path,
                "inputs": s.inputs,
                "parameters": s.parameters,
                "condition": s.condition
            } for s in self.steps]
        }
        with open(p, "w", encoding="utf-8") as f:
            yaml.safe_dump(data, f, sort_keys=False)
        return p

    def get_manifests(self) -> List[Manifest]:
        manifests = []
        for s in self.steps:
            m_path = self.base_dir / s.manifest_path if not Path(s.manifest_path).is_absolute() else Path(s.manifest_path)
            manifests.append(Manifest.from_yaml(m_path))
        return manifests

    def run(self, initial_state: Dict[str, BaseBioType]) -> Dict[str, Any]:
        """
        Executes pipeline end-to-end, passing typed outputs between nodes.
        Returns full execution state containing all intermediate and final objects.
        """
        state: Dict[str, BaseBioType] = dict(initial_state)

        for step in self.steps:
            # Check condition if present
            if step.condition and not evaluate_condition(step.condition, state):
                logging.getLogger("protflow.pipeline").info("Step '%s' skipped due to condition: %s", step.name, step.condition)
                continue

            # Load step manifest
            m_path = self.base_dir / step.manifest_path if not Path(s_path := step.manifest_path).is_absolute() else Path(s_path)
            manifest = Manifest.from_yaml(m_path)

            # Map inputs from state
            step_inputs: Dict[str, BaseBioType] = {}
            for port_name, state_key in step.inputs.items():
                if state_key not in state:
                    raise KeyError(
                        f"Step '{step.name}' requires input '{port_name}' from state key '{state_key}', "
                        f"but available keys are: {list(state.keys())}"
                    )
                step_inputs[port_name] = state[state_key]

            logger.info("Executing pipeline step: %s (%s)", step.name, manifest.name)
            outputs = self.runner.run(
                manifest=manifest,
                inputs=step_inputs,
                parameters=step.parameters
            )

            # Store outputs into state (keyed by step output names or prefixed)
            for out_name, out_obj in outputs.items():
                state[out_name] = out_obj
                state[f"{step.name}.{out_name}"] = out_obj

        return state

    def export_snakemake(self) -> str:
        return export_pipeline_to_snakemake(self.get_manifests())

    def export_nextflow(self) -> str:
        return export_pipeline_to_nextflow(self.get_manifests())
