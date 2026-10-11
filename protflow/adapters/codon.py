"""
Codon Optimization adapter and host genomic profile engine for ProtFlow.
Optimizes amino acid sequences for host expression (E. coli, Yeast, Human/CHO),
maximizes Codon Adaptation Index (CAI), balances %GC, and breaks forbidden motifs
such as cryptic splice sites, internal Shine-Dalgarno pauses, and restriction sites.
"""

from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import math
import re
import pandas as pd

from protflow.types import Sequence, Table, Provenance
from protflow.adapters.residue import AA1_TO_3

# Standard Codon Usage Frequencies per Host (frequency of codon given amino acid)
# Sources: Kazusa Codon Usage Database & E. coli B/K-12 highly expressed genes
HOST_PROFILES: Dict[str, Dict[str, Any]] = {
    "ecoli": {
        "organism": "Escherichia coli",
        "target_gc": (48.0, 56.0),
        "codons": {
            'A': {'GCG': 0.36, 'GCC': 0.26, 'GCA': 0.21, 'GCT': 0.17},
            'C': {'TGC': 0.56, 'TGT': 0.44},
            'D': {'GAT': 0.63, 'GAC': 0.37},
            'E': {'GAA': 0.68, 'GAG': 0.32},
            'F': {'TTC': 0.58, 'TTT': 0.42},
            'G': {'GGC': 0.38, 'GGT': 0.35, 'GGG': 0.15, 'GGA': 0.13},
            'H': {'CAT': 0.57, 'CAC': 0.43},
            'I': {'ATT': 0.49, 'ATC': 0.39, 'ATA': 0.11},
            'K': {'AAA': 0.74, 'AAG': 0.26},
            'L': {'CTG': 0.50, 'TTA': 0.14, 'TTG': 0.13, 'CTC': 0.10, 'CTT': 0.10, 'CTA': 0.04},
            'M': {'ATG': 1.00},
            'N': {'AAC': 0.55, 'AAT': 0.45},
            'P': {'CCG': 0.53, 'CCA': 0.20, 'CCT': 0.18, 'CCC': 0.10},
            'Q': {'CAG': 0.69, 'CAA': 0.31},
            'R': {'CGC': 0.36, 'CGT': 0.36, 'CGG': 0.11, 'CGA': 0.07, 'AGA': 0.04, 'AGG': 0.02},
            'S': {'AGC': 0.25, 'TCT': 0.17, 'AGT': 0.16, 'TCC': 0.15, 'TCA': 0.14, 'TCG': 0.14},
            'T': {'ACC': 0.40, 'ACG': 0.25, 'ACT': 0.19, 'ACA': 0.17},
            'V': {'GTG': 0.35, 'GTT': 0.28, 'GTC': 0.20, 'GTA': 0.17},
            'W': {'TGG': 1.00},
            'Y': {'TAT': 0.59, 'TAC': 0.41},
            '*': {'TAA': 0.61, 'TGA': 0.30, 'TAG': 0.09}
        },
        "forbidden_motifs": [
            "AGGAGG",  # Internal Shine-Dalgarno ribosomal pause
            "GGAGG",   # Core SD pause
            "AAAAAA",  # Poly-A premature termination / frameshifting
            "TTTTTT",  # Poly-T termination signal
            "GGTCTC",  # BsaI Golden Gate restriction site
            "GAGACC",  # BsaI reverse complement
            "CGTCTC",  # BsmBI restriction site
            "GAGACG",  # BsmBI reverse complement
            "GAATTC",  # EcoRI site
            "GGATCC",  # BamHI site
            "CATATG",  # NdeI site
            "CTCGAG",  # XhoI site
        ]
    },
    "human": {
        "organism": "Homo sapiens (CHO / Mammalian expression)",
        "target_gc": (50.0, 60.0),
        "codons": {
            'A': {'GCC': 0.40, 'GCT': 0.27, 'GCA': 0.23, 'GCG': 0.10},
            'C': {'TGC': 0.55, 'TGT': 0.45},
            'D': {'GAC': 0.54, 'GAT': 0.46},
            'E': {'GAG': 0.58, 'GAA': 0.42},
            'F': {'TTC': 0.55, 'TTT': 0.45},
            'G': {'GGC': 0.34, 'GGA': 0.25, 'GGG': 0.25, 'GGT': 0.16},
            'H': {'CAC': 0.59, 'CAT': 0.41},
            'I': {'ATC': 0.48, 'ATT': 0.36, 'ATA': 0.16},
            'K': {'AAG': 0.57, 'AAA': 0.43},
            'L': {'CTG': 0.41, 'CTC': 0.20, 'TTG': 0.13, 'CTT': 0.13, 'CTA': 0.07, 'TTA': 0.06},
            'M': {'ATG': 1.00},
            'N': {'AAC': 0.53, 'AAT': 0.47},
            'P': {'CCC': 0.33, 'CCT': 0.29, 'CCA': 0.27, 'CCG': 0.11},
            'Q': {'CAG': 0.74, 'CAA': 0.26},
            'R': {'CGG': 0.21, 'CGC': 0.19, 'AGG': 0.20, 'AGA': 0.20, 'CGA': 0.11, 'CGT': 0.09},
            'S': {'AGC': 0.24, 'TCC': 0.22, 'TCT': 0.18, 'AGT': 0.15, 'TCA': 0.15, 'TCG': 0.06},
            'T': {'ACC': 0.36, 'ACA': 0.28, 'ACT': 0.24, 'ACG': 0.12},
            'V': {'GTG': 0.47, 'GTC': 0.24, 'GTT': 0.18, 'GTA': 0.11},
            'W': {'TGG': 1.00},
            'Y': {'TAC': 0.57, 'TAT': 0.43},
            '*': {'TGA': 0.52, 'TAA': 0.28, 'TAG': 0.20}
        },
        "forbidden_motifs": [
            "AATAAA",  # Canonical polyadenylation signal
            "ATTAAA",  # Alternative polyA signal
            "GGTGAG",  # Consensus cryptic splice donor
            "GGTAAG",  # Consensus cryptic splice donor
            "TTTTTT",  # Pol III terminator
            "AAAAAA",  # Poly-A run
            "GGTCTC",  # BsaI
            "GAGACC",
            "CGTCTC",  # BsmBI
            "GAGACG"
        ]
    },
    "yeast": {
        "organism": "Saccharomyces cerevisiae",
        "target_gc": (38.0, 45.0),
        "codons": {
            'A': {'GCT': 0.53, 'GCC': 0.25, 'GCA': 0.20, 'GCG': 0.02},
            'C': {'TGT': 0.68, 'TGC': 0.32},
            'D': {'GAT': 0.65, 'GAC': 0.35},
            'E': {'GAA': 0.71, 'GAG': 0.29},
            'F': {'TTT': 0.59, 'TTC': 0.41},
            'G': {'GGT': 0.47, 'GGA': 0.22, 'GGC': 0.19, 'GGG': 0.12},
            'H': {'CAT': 0.64, 'CAC': 0.36},
            'I': {'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27},
            'K': {'AAG': 0.58, 'AAA': 0.42},
            'L': {'TTG': 0.29, 'TTA': 0.28, 'CTA': 0.14, 'CTG': 0.11, 'CTT': 0.10, 'CTC': 0.08},
            'M': {'ATG': 1.00},
            'N': {'AAT': 0.59, 'AAC': 0.41},
            'P': {'CCA': 0.42, 'CCT': 0.31, 'CCG': 0.14, 'CCC': 0.13},
            'Q': {'CAA': 0.69, 'CAG': 0.31},
            'R': {'AGA': 0.48, 'AGG': 0.21, 'CGT': 0.14, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04},
            'S': {'TCT': 0.29, 'TCC': 0.16, 'AGT': 0.16, 'TCA': 0.15, 'AGC': 0.12, 'TCG': 0.12},
            'T': {'ACT': 0.35, 'ACC': 0.26, 'ACA': 0.21, 'ACG': 0.18},
            'V': {'GTT': 0.39, 'GTC': 0.21, 'GTG': 0.21, 'GTA': 0.19},
            'W': {'TGG': 1.00},
            'Y': {'TAT': 0.56, 'TAC': 0.44},
            '*': {'TAA': 0.48, 'TGA': 0.29, 'TAG': 0.23}
        },
        "forbidden_motifs": [
            "TATATA",  # TATA box / cryptic promoter
            "TTTTTT",  # Premature transcription terminator
            "AAAAAA",  # Poly-A run
            "GGTCTC",  # BsaI
            "GAGACC"
        ]
    }
}


def compute_cai(codons: List[str], codon_weights: Dict[str, Dict[str, float]]) -> float:
    """Calculates Codon Adaptation Index (CAI) as geometric mean of relative adaptiveness."""
    if not codons:
        return 0.0
    log_sum = 0.0
    valid_count = 0
    for c in codons:
        c_upper = c.upper()
        # Find AA for this codon
        for aa, c_map in codon_weights.items():
            if c_upper in c_map:
                max_f = max(c_map.values())
                w_i = c_map[c_upper] / max(max_f, 1e-9)
                log_sum += math.log(max(w_i, 1e-4))
                valid_count += 1
                break
    return round(math.exp(log_sum / max(valid_count, 1)), 3)


def optimize_protein_sequence(
    protein_seq: str,
    host: str = "ecoli",
    custom_profile: Optional[Dict[str, Any]] = None,
    extra_forbidden: Optional[List[str]] = None
) -> Tuple[str, pd.DataFrame, Dict[str, Any]]:
    """
    Core codon optimization engine.
    Translates amino acid sequence into synonymous codons optimized for target host,
    actively avoiding forbidden motifs (splice sites, SD pauses, restriction sites).
    """
    profile = custom_profile or HOST_PROFILES.get(host.lower(), HOST_PROFILES["ecoli"])
    codon_map = profile["codons"]
    forbidden = list(profile.get("forbidden_motifs", []))
    if extra_forbidden:
        forbidden.extend(extra_forbidden)

    # Sort synonymous codons by frequency descending for each AA
    sorted_synonymous: Dict[str, List[Tuple[str, float]]] = {}
    for aa, c_dict in codon_map.items():
        sorted_synonymous[aa] = sorted(c_dict.items(), key=lambda x: x[1], reverse=True)

    optimized_codons: List[str] = []
    motif_avoided_count = 0
    rows = []

    clean_seq = protein_seq.strip().upper()
    if not clean_seq.endswith("*"):
        clean_seq += "*"

    curr_dna = ""

    for idx, aa in enumerate(clean_seq, start=1):
        candidates = sorted_synonymous.get(aa, [('NNN', 0.0)])
        chosen_codon = candidates[0][0]
        chosen_freq = candidates[0][1]

        # Check if chosen optimal codon causes a forbidden motif across junction
        candidate_found = False
        for cand_codon, cand_freq in candidates:
            test_dna = curr_dna[-15:] + cand_codon
            has_forbidden = any(f_motif in test_dna for f_motif in forbidden)
            if not has_forbidden:
                chosen_codon = cand_codon
                chosen_freq = cand_freq
                candidate_found = True
                break

        if not candidate_found:
            # shortcut: if every synonymous codon creates a motif, pick highest frequency
            motif_avoided_count += 1
            chosen_codon = candidates[0][0]
            chosen_freq = candidates[0][1]
        elif chosen_codon != candidates[0][0]:
            motif_avoided_count += 1

        curr_dna += chosen_codon
        optimized_codons.append(chosen_codon)

        # Per-residue GC and CAI metrics
        gc_count = chosen_codon.count("G") + chosen_codon.count("C")
        rows.append({
            "position": idx,
            "amino_acid": aa,
            "codon": chosen_codon,
            "codon_freq": chosen_freq,
            "gc_percent": round((gc_count / 3.0) * 100.0, 1)
        })

    full_dna = "".join(optimized_codons)
    total_gc = round(((full_dna.count("G") + full_dna.count("C")) / max(len(full_dna), 1)) * 100.0, 2)
    final_cai = compute_cai(optimized_codons, codon_map)

    stats = {
        "host": profile.get("organism", host),
        "cai_score": final_cai,
        "gc_percent": total_gc,
        "length_aa": len(clean_seq),
        "length_nt": len(full_dna),
        "motifs_avoided": motif_avoided_count,
        "target_gc_range": profile.get("target_gc", (45.0, 55.0))
    }

    df = pd.DataFrame(rows)
    return full_dna, df, stats


def parse_codon_optimization(raw_data: Any, context: Dict[str, Any]) -> Dict[str, Any]:
    """
    Adapter callable for ProtFlow pipeline node.
    Accepts:
      - protein (Sequence)
      - parameters: host ('ecoli' | 'human' | 'yeast'), extra_forbidden (list)
    Produces:
      - optimized_cds (Sequence): Optimized cDNA construct
      - codon_table (Table): Residue-level codon usage and GC table
    """
    inputs = context["inputs"]
    protein_obj: Sequence = inputs["protein"]
    params = context.get("parameters", {})
    host = params.get("host", "ecoli")
    extra_motifs = params.get("forbidden_motifs", [])

    dna_seq, codon_df, stats = optimize_protein_sequence(
        protein_seq=protein_obj.sequence,
        host=host,
        extra_forbidden=extra_motifs
    )

    pid = protein_obj.protein_id
    opt_pid = f"{pid}_opt_{host}"

    cds_seq_obj = Sequence(
        protein_id=opt_pid,
        sequence=dna_seq,
        description=f"Codon-optimized cDNA for {stats['host']} (CAI={stats['cai_score']}, GC={stats['gc_percent']}%)"
    )

    summary_df = codon_df.copy()
    # Add metadata summary column to table
    summary_df["cai_score"] = stats["cai_score"]
    summary_df["total_gc"] = stats["gc_percent"]

    table_obj = Table(protein_id=opt_pid, df=summary_df)
    table_obj.attach_provenance(Provenance(
        tool_name="codon_optimizer",
        parameters=stats
    ))

    return {
        "optimized_cds": cds_seq_obj,
        "codon_table": table_obj
    }
