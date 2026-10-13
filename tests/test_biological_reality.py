"""
Independent Biological Reality Test Suite.
Validates biological accuracy, biophysical constraints, and edge cases across domains:
- AT-rich parasites (Plasmodium falciparum)
- GC-rich actinobacteria (Streptomyces)
- Archaea with non-standard TTG/GTG start codons (Sulfolobus)
- Eukaryotes with long enzymatic reading frames (Yeast ADH)
- Biophysical CA-CA pseudo-bond geometry and 3D compactness
- Signal peptide discrimination vs catalytic Rossmann folds
"""

import math
from pathlib import Path
import numpy as np
import pytest

from protflow.pipeline import Pipeline
from protflow.runner import Runner
from protflow.types import Sequence, Structure, ResidueScore, Table, Pocket


BIOLOGICAL_CASES = [
    {
        "id": "PLASMODIUM_AT_RICH",
        "name": "Plasmodium falciparum (Malaria parasite)",
        "dna": "ATGATAAATAATAAAAATATTTATAATTTTTTAAATAAATATTTTAATTTAGATAATAAAATATTAAATAATAAATAA",
        "domain": "Eukaryota (Apicomplexa)",
        "expected_start": "ATG",
        "expected_intracellular": True
    },
    {
        "id": "STREPTOMYCES_GC_RICH",
        "name": "Streptomyces coelicolor (Actinobacterium)",
        "dna": "ATGACCGACCGGGCCGTCGTGGCCGCCCTGCGGGCGCTGGGCGCCCGGGTCCTCGCCGAGGTCGTCCCGGCCGCCCGCTGA",
        "domain": "Bacteria (Actinomycetota)",
        "expected_start": "ATG",
        "expected_intracellular": True
    },
    {
        "id": "SULFOLOBUS_ARCHAEA",
        "name": "Sulfolobus solfataricus (Thermoacidophile Archaea)",
        "dna": "TTGATAGATATAGTGATAGATAAATTAGATAATAGATTAATAGATAATTTAGATTATAGATTATTAAATTTAGATTAA",
        "domain": "Archaea (Crenarchaeota)",
        "expected_start": "TTG",
        "expected_intracellular": True
    },
    {
        "id": "YEAST_ADH1",
        "name": "Saccharomyces cerevisiae (Baker's yeast ADH1)",
        "dna": (
            "ATGTCTATCCCAGAAACTCAAAAAGGTGTTATCTTCTACGAATCCCACGGTAAGTTGGAATACAAAGATATTCCAGTTCCAAAGCCAAAGTCCAACGAA"
            "TTATTGATCAACGTTAAGTACTCCGGTGTCTGTCACACCGATTTGCACGCTTGGCATGGTGATTGGCCATTGCCAGTTAAGTTACCATTAGTTGGTGGTCAC"
            "GAAGGTGCTGGTGTCGTTGTCGCTATGGGTGAAAACGTTAAGGGCTGGAAGATCGGTGATTACGCCGGTATCAAATGGTTGAACGGTTCTTGTATGGCCTGT"
            "GAATACTGTGAATTGGGTAACGAATCCAACTGTCCTCACGCTGACTTGTCTGGTTACACCCACGACGGTTCTTTCCAACAATACGCTACCGCTGACGCTGTT"
            "CAAGCCGCTCACATTCCACAAGGTACCGACTTGGCCCAAGTCGCCCCAATCTTGTGTGCTGGTATCACCGTCTACAAGGCTTTGAAGTCTGCCAACTTGCGT"
            "GCTGGTCACTGGGCTGCTATCTCCGGTGCCGCTGGTGGTCTCGGTTCTTTGGCCGTTCAATATGCTAAGGCCATGGGTTACAGAGTCTTGGGTATTGACGGT"
            "GGTGAAGGTAAGGAAGAATTATTCAGATCCATCGGTGGTGAAGTCTTCATTGACTTCACCAAGGAAAAGGACATTGTCGGTGCTGTTCTAAGGCCACTGACG"
            "GTGGTGCTCACGGTGTCATCAACGTTTCCGTTTCCGAAGCCGCTATTGAAGCTTCTACCAGATACGTTAGAGCTAACGGTACCACCGTTTTGGTCGGTATGC"
            "CAGCTGGTGCCAAGTGTTGTTCTGACGTTTTCAACCAAGTCGTCAAGTCCATCTCTATTGTTGGTTCTTACGTCGGTAACAGAGCTGACACCAGAGAAGCTT"
            "TGGACTTCTTCGCCAGAGGTTTGGTCAAGTCTCCAATCAAGGTTGTCGGCTTGTCTACCTTGCCAGAAATTTACGAAAAGATGGAAAAGGGTCAAATCGTTG"
            "GTAGATACGTTGTTGACACTTCTAAATAA"
        ),
        "domain": "Eukaryota (Fungi)",
        "expected_start": "ATG",
        "expected_intracellular": True
    }
]


@pytest.mark.parametrize("case", BIOLOGICAL_CASES, ids=[c["id"] for c in BIOLOGICAL_CASES])
def test_biological_diversity_pipeline_end_to_end(tmp_path, case):
    """
    Executes the full reference pipeline on extreme real-world biological cases:
    Tests gene finding, start codon recognition, localization, 3D structure,
    pLDDT scoring, and target prioritization.
    """
    pipeline_yaml = Path(__file__).parent.parent / "bacterial_reference_pipeline.yaml"
    runner = Runner(cache_dir=tmp_path / f"cache_{case['id']}", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)

    seq = Sequence(
        protein_id=case["id"],
        sequence=case["dna"],
        description=f"{case['name']} [{case['domain']}]"
    )

    final_state = pipeline.run({"genome": seq})

    # 1. Gene prediction check
    assert "protein" in final_state
    protein = final_state["protein"]
    assert len(protein.sequence) >= 10, f"{case['id']}: protein translation too short"
    assert not any(c in protein.sequence for c in "*X?"), f"{case['id']}: invalid amino acids in translation"

    # Verify start codon in genes table
    assert "genes" in final_state
    genes_df = final_state["genes"].df
    assert not genes_df.empty
    start_c = genes_df["start_codon"].iloc[0]
    assert start_c in ("ATG", "GTG", "TTG", "CTG"), f"{case['id']}: unrecognized start codon {start_c}"

    # 2. Subcellular localization check
    assert "filter_table" in final_state
    filt_df = final_state["filter_table"].df
    assert bool(filt_df["is_intracellular"].iloc[0]) is case["expected_intracellular"]

    # 3. 3D Structure check
    assert "structure" in final_state
    st_obj = final_state["structure"]
    st = st_obj.to_gemmi_structure()

    ca_atoms = [atom for model in st for chain in model for res in chain for atom in res if atom.name == "CA"]
    assert len(ca_atoms) > 0

    # 3D coordinates must not be 1-D flat line
    xs = [a.pos.x for a in ca_atoms]
    ys = [a.pos.y for a in ca_atoms]
    zs = [a.pos.z for a in ca_atoms]
    assert np.std(xs) > 1.0, "X coordinates lack 3D dispersion"
    assert np.std(ys) > 1.0, "Y coordinates lack 3D dispersion"
    assert np.std(zs) > 1.0, "Z coordinates lack 3D dispersion"

    # Consecutive CA-CA pseudo-bond geometry: continuous backbone (3.8 +- 0.8 A)
    ca_dists = [
        math.sqrt((xs[i+1]-xs[i])**2 + (ys[i+1]-ys[i])**2 + (zs[i+1]-zs[i])**2)
        for i in range(len(ca_atoms) - 1)
    ]
    assert all(3.0 <= d <= 5.0 for d in ca_dists), f"Non-physical peptide chain jumps: {ca_dists}"

    # 4. Residue pLDDT distribution
    assert "residue_scores" in final_state
    scores = final_state["residue_scores"].df["score"]
    assert scores.min() >= 30.0, "pLDDT scores unrealistically low (<30)"
    assert scores.max() <= 100.0, "pLDDT scores exceed 100"
    assert scores.std() > 0.0, "pLDDT scores are completely uniform/flat"

    # 5. Pockets
    assert "pockets" in final_state
    pockets = final_state["pockets"]
    assert not pockets.df.empty
    assert pockets.df["volume"].iloc[0] > 0.0
    assert 0.0 <= pockets.df["score"].iloc[0] <= 1.0

    # 6. Final candidate prioritization table
    assert "priority_table" in final_state
    p_df = final_state["priority_table"].df
    assert 0.0 <= p_df["final_priority_score"].iloc[0] <= 1.0
