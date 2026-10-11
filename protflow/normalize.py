"""
Node-level score normalization utilities for ProtFlow.
Ensures heterogeneous tool metrics (probabilities, geometric scores, SASA, pLDDT)
are mapped to compatible [0.0, 1.0] scales without distortion.
"""

from __future__ import annotations
import math
from typing import Union, List
import numpy as np


def minmax_norm(val: float, vmin: float, vmax: float) -> float:
    """Clips and scales value linearly between 0.0 and 1.0."""
    if vmax <= vmin:
        return 0.5
    clipped = max(vmin, min(vmax, val))
    return float((clipped - vmin) / (vmax - vmin))


def sigmoid_norm(val: float, midpoint: float = 0.5, slope: float = 10.0) -> float:
    """Logistic sigmoid transformation."""
    try:
        return float(1.0 / (1.0 + math.exp(-slope * (val - midpoint))))
    except OverflowError:
        return 0.0 if val < midpoint else 1.0


def normalize_plddt(mean_plddt: float) -> float:
    """
    Normalizes AlphaFold pLDDT score (0 to 100):
    < 50: very low confidence (< 0.2)
    50-70: low (0.2 - 0.5)
    70-90: confident (0.5 - 0.85)
    > 90: very high (> 0.85)
    """
    if mean_plddt >= 90.0:
        return minmax_norm(mean_plddt, 90.0, 100.0) * 0.15 + 0.85
    elif mean_plddt >= 70.0:
        return minmax_norm(mean_plddt, 70.0, 90.0) * 0.35 + 0.50
    elif mean_plddt >= 50.0:
        return minmax_norm(mean_plddt, 50.0, 70.0) * 0.30 + 0.20
    else:
        return minmax_norm(mean_plddt, 0.0, 50.0) * 0.20


def normalize_pocket_score(raw_score: float, tool: str = "p2rank") -> float:
    """
    Normalizes pocket detector scores into standard 0.0 - 1.0 scale:
    - P2Rank: outputs calibrated probability [0.0, 1.0]
    - fpocket: raw druggability or score [0.0, 1.0] or [-5.0, 50.0]
    """
    t = tool.lower()
    if "p2rank" in t:
        return minmax_norm(raw_score, 0.0, 1.0)
    elif "fpocket" in t:
        if raw_score <= 1.0:
            return minmax_norm(raw_score, 0.0, 1.0)
        return sigmoid_norm(raw_score, midpoint=20.0, slope=0.2)
    return minmax_norm(raw_score, 0.0, 1.0)


def normalize_solubility_from_sasa(hydrophobic_ratio: float, ideal_target: float = 0.40) -> float:
    """
    Estimates solubility from surface hydrophobic exposure:
    Penalizes deviations from physiological globular ratio (~0.35-0.45).
    """
    dev = abs(hydrophobic_ratio - ideal_target)
    return max(0.0, min(1.0, float(1.0 - dev * 2.5)))
