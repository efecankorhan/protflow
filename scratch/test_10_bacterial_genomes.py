"""
Comprehensive Bacterial Genome Scale (150,000+ bp) Benchmark and Reality Audit.
Evaluates 10 diverse bacterial species across the tree of life up to localization filter,
compares biological statistics against empirical literature, and audits 1 randomly selected
protein through AlphaFold 3D structure prediction and P2Rank pocket finding.
"""

import json
import math
import random
from pathlib import Path
import numpy as np
import pandas as pd
import gemmi

from protflow.pipeline import Pipeline, PipelineStep
from protflow.runner import Runner
from protflow.types import Sequence, SequenceCollection, Table
from protflow.adapters.genomics import fetch_ncbi_genome
from protflow.adapters.structure_fetch import parse_alphafold_fetcher
from protflow.adapters.pockets import parse_p2rank


ORGANISMS = [
    {
        "id": "NC_000913.3",
        "name": "Escherichia coli K-12 MG1655",
        "phylum": "Pseudomonadota (Gammaproteobacteria)",
        "expected_gc": 50.8,
        "lit_sec_range": (4.0, 12.0)  # empirical % secreted/periplasmic in E. coli
    },
    {
        "id": "NC_000964.3",
        "name": "Bacillus subtilis 168",
        "phylum": "Bacillota (Firmicutes)",
        "expected_gc": 43.5,
        "lit_sec_range": (3.0, 12.0)
    },
    {
        "id": "NC_002516.2",
        "name": "Pseudomonas aeruginosa PAO1",
        "phylum": "Pseudomonadota (Gammaproteobacteria)",
        "expected_gc": 66.6,
        "lit_sec_range": (5.0, 15.0)
    },
    {
        "id": "NC_007795.1",
        "name": "Staphylococcus aureus USA300",
        "phylum": "Bacillota (Staphylococcaceae)",
        "expected_gc": 32.8,
        "lit_sec_range": (4.0, 14.0)
    },
    {
        "id": "NC_000962.3",
        "name": "Mycobacterium tuberculosis H37Rv",
        "phylum": "Actinomycetota (Mycobacteriaceae)",
        "expected_gc": 65.6,
        "lit_sec_range": (3.0, 10.0)
    },
    {
        "id": "NC_000915.1",
        "name": "Helicobacter pylori 26695",
        "phylum": "Campylobacterota",
        "expected_gc": 38.9,
        "lit_sec_range": (3.0, 10.0)
    },
    {
        "id": "NC_003028.3",
        "name": "Streptococcus pneumoniae TIGR4",
        "phylum": "Bacillota (Streptococcaceae)",
        "expected_gc": 39.7,
        "lit_sec_range": (4.0, 12.0)
    },
    {
        "id": "NC_006461.1",
        "name": "Thermus thermophilus HB8",
        "phylum": "Deinococcota",
        "expected_gc": 69.5,
        "lit_sec_range": (3.0, 10.0)
    },
    {
        "id": "NC_001318.1",
        "name": "Borrelia burgdorferi B31",
        "phylum": "Spirochaetota",
        "expected_gc": 28.6,
        "lit_sec_range": (4.0, 15.0)
    },
    {
        "id": "NC_003888.3",
        "name": "Streptomyces coelicolor A3(2)",
        "phylum": "Actinomycetota",
        "expected_gc": 72.1,
        "lit_sec_range": (4.0, 14.0)
    },
    {
        "id": "NC_002505.1",
        "name": "Vibrio cholerae O1 El Tor",
        "phylum": "Pseudomonadota (Vibrionaceae)",
        "expected_gc": 47.5,
        "lit_sec_range": (4.0, 12.0)
    },
    {
        "id": "NC_003210.1",
        "name": "Listeria monocytogenes EGD-e",
        "phylum": "Bacillota (Listeriaceae)",
        "expected_gc": 38.0,
        "lit_sec_range": (3.0, 11.0)
    }
]


def run_benchmark():
    work_dir = Path("scratch/benchmark_run")
    work_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(cache_dir=work_dir / "cache", enable_cache=True)

    results = []
    all_candidate_proteins = []

    print(f"=== Starting Bacterial Genomes Benchmark ({len(ORGANISMS)} organisms, 150,000+ bp each) ===\n")

    for idx, org in enumerate(ORGANISMS, 1):
        acc = org["id"]
        print(f"[{idx}/{len(ORGANISMS)}] Fetching & Analyzing {org['name']} ({acc})...")
        res, err = fetch_ncbi_genome(acc, max_bp=150000)
        if not res or err:
            print(f"  Warning: NCBI fetch issue for {acc}: {err}. Skipping...")
            continue

        seq_dna = res["sequence"]
        actual_len = len(seq_dna)
        gc_pct = res["gc_content"]
        org_title = res["organism"] or org["name"]

        # Run pipeline up to localization filter
        genome_seq = Sequence(
            protein_id=acc.replace(".", "_"),
            sequence=seq_dna,
            description=f"{org_title} | {org['phylum']}"
        )

        steps = [
            PipelineStep(
                name="gene_prediction",
                manifest_path="manifests/01_augustus.yaml",
                inputs={"genome": "genome"},
                parameters={"species": "generic"}
            ),
            PipelineStep(
                name="localization_filter",
                manifest_path="manifests/02_subcellular_filter.yaml",
                inputs={"proteins": "proteins"}
            )
        ]
        pipeline = Pipeline(name=f"pipeline_{acc}", steps=steps, runner=runner)
        state = pipeline.run({"genome": genome_seq})

        genes_df = state["genes"].df
        filter_df = state["filter_table"].df
        intracellular_coll = state["proteins"]

        # Biological metrics
        num_genes = len(genes_df)
        gene_density_per_100kb = round((num_genes / actual_len) * 100_000, 1)

        # Total coding nucleotides (sum of gene lengths * 3)
        total_coding_nt = genes_df["length_aa"].sum() * 3
        coding_density_pct = round((total_coding_nt / actual_len) * 100.0, 1)

        # Start codon frequencies
        start_counts = genes_df["start_codon"].value_counts().to_dict()
        atg_pct = round((start_counts.get("ATG", 0) / max(num_genes, 1)) * 100.0, 1)
        gtg_pct = round((start_counts.get("GTG", 0) / max(num_genes, 1)) * 100.0, 1)
        ttg_pct = round((start_counts.get("TTG", 0) / max(num_genes, 1)) * 100.0, 1)

        # Localization metrics
        loc_counts = filter_df["localization"].value_counts().to_dict()
        cyto_count = loc_counts.get("Cytoplasm", 0)
        sec_count = loc_counts.get("Extracellular", 0)
        sec_pct = round((sec_count / max(num_genes, 1)) * 100.0, 1)
        cyto_pct = round((cyto_count / max(num_genes, 1)) * 100.0, 1)

        mean_len_aa = round(float(genes_df["length_aa"].mean()), 1)

        # Add passed proteins to candidate pool
        for pid, p_seq in intracellular_coll.sequences.items():
            if len(p_seq.sequence) >= 50:
                all_candidate_proteins.append({
                    "organism": org_title,
                    "acc": acc,
                    "protein": p_seq
                })

        entry = {
            "organism": org_title,
            "accession": acc,
            "phylum": org["phylum"],
            "bp_length": actual_len,
            "gc_content_pct": gc_pct,
            "expected_gc": org["expected_gc"],
            "num_genes": num_genes,
            "gene_density_per_100kb": gene_density_per_100kb,
            "coding_density_pct": coding_density_pct,
            "mean_protein_len_aa": mean_len_aa,
            "start_codon_atg_pct": atg_pct,
            "start_codon_gtg_pct": gtg_pct,
            "start_codon_ttg_pct": ttg_pct,
            "intracellular_pct": cyto_pct,
            "secreted_pct": sec_pct,
            "lit_sec_range": org["lit_sec_range"],
            "sec_in_lit_range": org["lit_sec_range"][0] <= sec_pct <= org["lit_sec_range"][1]
        }
        results.append(entry)
        print(f"  -> Genes: {num_genes} | Coding: {coding_density_pct}% | Mean len: {mean_len_aa} aa | Secreted: {sec_pct}% (Lit: {org['lit_sec_range']}%)")

    # Save summary table
    df_results = pd.DataFrame(results)
    csv_path = work_dir / "bacterial_10_genome_benchmark.csv"
    df_results.to_csv(csv_path, index=False)
    print(f"\n[Saved benchmark summary to {csv_path}]")

    # Pick 1 random protein for AlphaFold 3D structure prediction and audit
    random.seed(42)  # reproducible random selection
    chosen = random.choice(all_candidate_proteins)
    chosen_prot = chosen["protein"]
    print(f"\n=== Selected 1 Random Protein for AlphaFold Audit ===")
    print(f"Organism: {chosen['organism']} ({chosen['acc']})")
    print(f"Protein ID: {chosen_prot.protein_id}")
    print(f"Length: {len(chosen_prot.sequence)} aa")
    print(f"Description: {chosen_prot.description}")

    # Run structure prediction via parse_alphafold_fetcher
    af_context = {
        "inputs": {"protein": chosen_prot},
        "output_dir": str(work_dir / "alphafold_output"),
        "parameters": {"allow_synthetic_fallback": True}
    }
    af_res = parse_alphafold_fetcher(None, af_context)
    st_obj = af_res["structure"]
    scores_obj = af_res["residue_scores"]

    st = st_obj.to_gemmi_structure()
    ca_atoms = [atom for model in st for chain in model for res in chain for atom in res if atom.name == "CA"]

    # Coordinates
    xs = np.array([a.pos.x for a in ca_atoms])
    ys = np.array([a.pos.y for a in ca_atoms])
    zs = np.array([a.pos.z for a in ca_atoms])

    ca_dists = [
        math.sqrt((xs[i+1]-xs[i])**2 + (ys[i+1]-ys[i])**2 + (zs[i+1]-zs[i])**2)
        for i in range(len(ca_atoms)-1)
    ]

    plddt_df = scores_obj.df
    plddt_scores = plddt_df["score"]

    # Run P2Rank pocket finder on this structure
    p2rank_context = {
        "inputs": {"structure": st_obj},
        "output_dir": str(work_dir / "p2rank_output"),
        "parameters": {}
    }
    pocket_res = parse_p2rank(None, p2rank_context)
    pockets = pocket_res["pockets"]

    audit_report = {
        "protein_id": chosen_prot.protein_id,
        "organism": chosen["organism"],
        "sequence_len": len(chosen_prot.sequence),
        "structure_path": st_obj.file_path,
        "format": st_obj.format,
        "model_source": st_obj.provenance.parameters.get("model_source"),
        "num_ca_atoms": len(ca_atoms),
        "mean_ca_ca_bond_A": round(float(np.mean(ca_dists)), 3),
        "min_ca_ca_bond_A": round(float(min(ca_dists)), 3),
        "max_ca_ca_bond_A": round(float(max(ca_dists)), 3),
        "bond_length_conforms_to_physics": (3.0 <= min(ca_dists) and max(ca_dists) <= 5.0),
        "std_x_dispersion": round(float(np.std(xs)), 2),
        "std_y_dispersion": round(float(np.std(ys)), 2),
        "std_z_dispersion": round(float(np.std(zs)), 2),
        "plddt_mean": round(float(plddt_scores.mean()), 2),
        "plddt_std": round(float(plddt_scores.std()), 2),
        "plddt_min": round(float(plddt_scores.min()), 2),
        "plddt_max": round(float(plddt_scores.max()), 2),
        "num_pockets": len(pockets.df),
        "top_pocket_volume": float(pockets.df["volume"].iloc[0]) if not pockets.df.empty else 0.0,
        "top_pocket_score": float(pockets.df["score"].iloc[0]) if not pockets.df.empty else 0.0,
    }

    audit_json_path = work_dir / "alphafold_audit_report.json"
    audit_json_path.write_text(json.dumps(audit_report, indent=2), encoding="utf-8")
    print(f"\n[Saved AlphaFold audit report to {audit_json_path}]")
    print(json.dumps(audit_report, indent=2))

    return df_results, audit_report


if __name__ == "__main__":
    run_benchmark()
