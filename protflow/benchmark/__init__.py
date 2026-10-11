"""
ProtFlow benchmark, calibration, and optimization suite.
"""

from protflow.benchmark.dataset import (
    BenchmarkEntry,
    SWISS_PROT_BENCHMARK,
    load_benchmark_dataframe,
    get_benchmark_sequences
)
from protflow.benchmark.calibration import CalibrationTable
from protflow.benchmark.optimizer import (
    VariantConfig,
    VariantResult,
    PipelineOptimizer
)
from protflow.benchmark.evaluator import (
    ToolBenchmarkReport,
    NodeEvaluator
)

__all__ = [
    "BenchmarkEntry",
    "SWISS_PROT_BENCHMARK",
    "load_benchmark_dataframe",
    "get_benchmark_sequences",
    "CalibrationTable",
    "VariantConfig",
    "VariantResult",
    "PipelineOptimizer",
    "ToolBenchmarkReport",
    "NodeEvaluator",
]
