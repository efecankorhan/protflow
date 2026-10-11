"""
Curated Swiss-Prot enzyme benchmark dataset for pipeline calibration and evaluation.
Includes characterized enzymes across 6 EC classes and inactive/non-catalytic controls.
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd
from protflow.types import Sequence


@dataclass
class BenchmarkEntry:
    protein_id: str
    uniprot_acc: str
    name: str
    ec_class: str  # EC 1-6 or "NON_ENZYME"
    is_active_enzyme: bool
    sequence: str
    known_pocket_residues: str
    description: str = ""


# Representative Swiss-Prot benchmark reference set
SWISS_PROT_BENCHMARK: List[BenchmarkEntry] = [
    # EC 1: Oxidoreductase
    BenchmarkEntry(
        protein_id="BENCH_EC1_ADH",
        uniprot_acc="P00325",
        name="Alcohol dehydrogenase 1A",
        ec_class="EC 1 (Oxidoreductase)",
        is_active_enzyme=True,
        sequence="MSTAGKVIKCKAAVLWEVKKPFSIEDVEVAPPKAHEVRIKMVATGICRSDDHVVSGTLVTP",
        known_pocket_residues="A:15,A:16,A:40,A:41",
        description="Zinc-dependent alcohol dehydrogenase"
    ),
    # EC 2: Transferase
    BenchmarkEntry(
        protein_id="BENCH_EC2_PGK",
        uniprot_acc="P0A799",
        name="Phosphoglycerate kinase",
        ec_class="EC 2 (Transferase)",
        is_active_enzyme=True,
        sequence="MSLSSKLSVQDLDLKDKRVFIRVDFNVPLDGKKITSNQRIVAALPTIKYVLEHHPRYVVLAS",
        known_pocket_residues="A:21,A:22,A:23,A:38",
        description="Glycolytic transferase with ATP-binding pocket"
    ),
    # EC 3: Hydrolase
    BenchmarkEntry(
        protein_id="BENCH_EC3_LYZ",
        uniprot_acc="P00698",
        name="Lysozyme C",
        ec_class="EC 3 (Hydrolase)",
        is_active_enzyme=True,
        sequence="KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWW",
        known_pocket_residues="A:35,A:52,A:53,A:62",
        description="Catalytic cleft with Glu35 and Asp52"
    ),
    # EC 4: Lyase
    BenchmarkEntry(
        protein_id="BENCH_EC4_TNAA",
        uniprot_acc="P0A855",
        name="Tryptophanase",
        ec_class="EC 4 (Lyase)",
        is_active_enzyme=True,
        sequence="MENFKHLPEPFRIRVIEPVKRTTRAYREEAIIKSGMNPFLLDSEDVFIDLLTDSGTGAVTQSM",
        known_pocket_residues="A:24,A:25,A:26,A:55",
        description="Pyridoxal phosphate-dependent lyase"
    ),
    # EC 5: Isomerase
    BenchmarkEntry(
        protein_id="BENCH_EC5_TIM",
        uniprot_acc="P00939",
        name="Triosephosphate isomerase",
        ec_class="EC 5 (Isomerase)",
        is_active_enzyme=True,
        sequence="MRHPLVMGNNKLAEMGAKLLGELGAAQVVLCHAPTLEQLTELAAEGADVWVIGEPGEL",
        known_pocket_residues="A:10,A:12,A:45,A:46",
        description="TIM barrel catalytic pocket"
    ),
    # EC 6: Ligase
    BenchmarkEntry(
        protein_id="BENCH_EC6_GLNA",
        uniprot_acc="P0A9C9",
        name="Glutamine synthetase",
        ec_class="EC 6 (Ligase)",
        is_active_enzyme=True,
        sequence="MSAEHVLTMLNEHEVKFVDLRFTDTKGKEQHVTIPAHQVNSEFEEKDGKTYKIRITFTSP",
        known_pocket_residues="A:18,A:19,A:32,A:33",
        description="ATP-dependent glutamine ligase"
    ),
    # Negative Control 1: Collagen/Fibrous protein (non-catalytic)
    BenchmarkEntry(
        protein_id="BENCH_NEG_COL",
        uniprot_acc="P02452",
        name="Collagen alpha-1(I) chain",
        ec_class="NON_ENZYME",
        is_active_enzyme=False,
        sequence="GPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPPGPP",
        known_pocket_residues="",
        description="Triple-helical structural fibrous protein with no catalytic pocket"
    ),
    # Negative Control 2: Bacterial Flagellin (structural motor filament)
    BenchmarkEntry(
        protein_id="BENCH_NEG_FLA",
        uniprot_acc="P04949",
        name="Flagellin FliC",
        ec_class="NON_ENZYME",
        is_active_enzyme=False,
        sequence="MAQVINTNSLSLITQNNINKNQSALSSSIERLSSGLRINSAKDDAAGQAIANRFTSNIKGL",
        known_pocket_residues="",
        description="Filamentous structural protein, no enzymatic active pocket"
    ),
    # Negative Control 3: Disordered peptide / Decoy
    BenchmarkEntry(
        protein_id="BENCH_NEG_DIS",
        uniprot_acc="P99999",
        name="Intrinsically Disordered Decoy",
        ec_class="NON_ENZYME",
        is_active_enzyme=False,
        sequence="GSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGSGS",
        known_pocket_residues="",
        description="Low complexity unstructured repeat with zero structural pockets"
    ),
]


def load_benchmark_dataframe() -> pd.DataFrame:
    """Loads Swiss-Prot benchmark set as a structured pandas DataFrame."""
    return pd.DataFrame([
        {
            "protein_id": e.protein_id,
            "uniprot_acc": e.uniprot_acc,
            "name": e.name,
            "ec_class": e.ec_class,
            "is_active_enzyme": e.is_active_enzyme,
            "sequence": e.sequence,
            "known_pocket_residues": e.known_pocket_residues,
            "description": e.description
        }
        for e in SWISS_PROT_BENCHMARK
    ])


def get_benchmark_sequences() -> List[Sequence]:
    """Returns benchmark entries as ProtFlow Sequence entities."""
    return [
        Sequence(protein_id=e.protein_id, sequence=e.sequence, description=f"{e.name} ({e.ec_class})")
        for e in SWISS_PROT_BENCHMARK
    ]
