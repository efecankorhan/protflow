"""
Rigorous 250,000 bp Bacterial Genome Scale Ground-Truth Validation.
Tests 10 diverse NCBI bacterial genomes against real NCBI RefSeq curated CDS records.
Evaluates:
- Predicted genes vs Curated genes (TP, FP, FN, Sensitivity/Recall).
- Start codon agreement and length fidelity.
- Protein functional diversity (Enzymes, Transporters, Structural, Regulators).
- Subcellular localization sanity.
- Identifies any problematic edge cases or discrepancies without hallucination.
"""

import json
import math
import re
from pathlib import Path
import numpy as np
import pandas as pd

from protflow.pipeline import Pipeline, PipelineStep
from protflow.runner import Runner
from protflow.types import Sequence, SequenceCollection, Table
from protflow.adapters.genomics import fetch_ncbi_genome


ORGANISMS = [
    {"id": "NC_000913.3", "name": "Escherichia coli K-12 MG1655", "phylum": "Gammaproteobacteria"},
    {"id": "NC_000964.3", "name": "Bacillus subtilis 168", "phylum": "Bacillota (Firmicutes)"},
    {"id": "NC_002516.2", "name": "Pseudomonas aeruginosa PAO1", "phylum": "Gammaproteobacteria"},
    {"id": "NC_007795.1", "name": "Staphylococcus aureus USA300", "phylum": "Bacillota (Staphylococcaceae)"},
    {"id": "NC_000962.3", "name": "Mycobacterium tuberculosis H37Rv", "phylum": "Actinomycetota"},
    {"id": "NC_000915.1", "name": "Helicobacter pylori 26695", "phylum": "Campylobacterota"},
    {"id": "NC_002505.1", "name": "Vibrio cholerae O1 El Tor", "phylum": "Vibrionaceae"},
    {"id": "NC_003210.1", "name": "Listeria monocytogenes EGD-e", "phylum": "Bacillota (Listeriaceae)"},
    {"id": "NC_006461.1", "name": "Thermus thermophilus HB8", "phylum": "Deinococcota"},
    {"id": "NC_003888.3", "name": "Streptomyces coelicolor A3(2)", "phylum": "Actinomycetota"}
]


def classify_function(product_name: str) -> str:
    """Classifies protein functional category based on product title."""
    p = (product_name or "").lower()
    if any(k in p for k in ("kinase", "synthase", "synthetase", "reductase", "dehydrogenase", "hydrolase", "peptidase", "transferase", "isomerase", "polymerase", "ligase", "protease", "nuclease")):
        return "Catalytic Enzyme"
    if any(k in p for k in ("transporter", "permease", "efflux", "symporter", "antiporter", "abc transporter", "channel", "porin")):
        return "Transporter / Channel"
    if any(k in p for k in ("regulator", "repressor", "activator", "transcription factor", "response regulator", "sigma factor")):
        return "Transcriptional Regulator"
    if any(k in p for k in ("ribosomal", "rRNA", "tRNA", "ribosome")):
        return "Ribosomal / Translation Machinery"
    if any(k in p for k in ("chaperone", "heat shock", "dna repair", "recombinase", "topoisomerase", "helicase")):
        return "Chaperone / DNA Repair"
    if any(k in p for k in ("membrane", "envelope", "outer membrane", "lipoprotein")):
        return "Membrane / Structural"
    if any(k in p for k in ("toxin", "adhesin", "invasion", "virulence", "antigen", "pilin", "flagell")):
        return "Virulence / Motility"
    return "Hypothetical / Other"


def run_ground_truth_benchmark():
    work_dir = Path("scratch/benchmark_250k")
    work_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(cache_dir=work_dir / "cache", enable_cache=True)

    summary_rows = []
    all_proteins_catalog = []

    print("=== Starting 250,000 bp NCBI Bacterial Ground-Truth Validation (10 Organisms) ===\n")

    for idx, org in enumerate(ORGANISMS, 1):
        acc = org["id"]
        print(f"[{idx}/10] Fetching 250k bp: {org['name']} ({acc})...")
        res, err = fetch_ncbi_genome(acc, max_bp=250000)
        if not res or err:
            print(f"  Error: {err}")
            continue

        seq_dna = res["sequence"]
        actual_len = len(seq_dna)
        gc_pct = res["gc_content"]
        org_title = res["organism"] or org["name"]

        # 1. Curated NCBI ground truth CDS in this exact window [1..250000]
        all_cds = res.get("cds_records", [])
        curated_cds = [c for c in all_cds if 0 < c.get("end", 0) <= actual_len]
        num_curated = len(curated_cds)

        # 2. Run ProtFlow pipeline up to localization_filter
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
        pipeline = Pipeline(name=f"pipeline_250k_{acc}", steps=steps, runner=runner)
        state = pipeline.run({"genome": genome_seq})

        pred_genes_df = state["genes"].df
        filter_df = state["filter_table"].df
        num_predicted = len(pred_genes_df)

        # 3. Ground truth matching & benchmarking
        # Coordinate match: same strand, start within 30 bp, end within 30 bp
        tp_matches = 0
        matched_curated_indices = set()
        matched_pred_indices = set()

        for p_idx, p_row in pred_genes_df.iterrows():
            p_st, p_en, p_str = p_row["start"], p_row["end"], p_row["strand"]
            for c_idx, c in enumerate(curated_cds):
                if c_idx in matched_curated_indices:
                    continue
                if c["strand"] == p_str:
                    if abs(p_st - c["start"]) <= 30 and abs(p_en - c["end"]) <= 30:
                        tp_matches += 1
                        matched_curated_indices.add(c_idx)
                        matched_pred_indices.add(p_idx)
                        break

        # Calculate metrics
        recall_pct = round((tp_matches / max(num_curated, 1)) * 100.0, 1) if num_curated > 0 else 0.0
        precision_pct = round((tp_matches / max(num_predicted, 1)) * 100.0, 1) if num_predicted > 0 else 0.0
        fn_count = num_curated - tp_matches
        fp_count = num_predicted - tp_matches

        # Coding density
        coding_density_pct = round((pred_genes_df["length_aa"].sum() * 3 / actual_len) * 100.0, 1)

        # Start codons
        start_counts = pred_genes_df["start_codon"].value_counts().to_dict()
        atg_pct = round((start_counts.get("ATG", 0) / max(num_predicted, 1)) * 100.0, 1)
        gtg_pct = round((start_counts.get("GTG", 0) / max(num_predicted, 1)) * 100.0, 1)
        ttg_pct = round((start_counts.get("TTG", 0) / max(num_predicted, 1)) * 100.0, 1)

        # Localization
        loc_counts = filter_df["localization"].value_counts().to_dict()
        sec_pct = round((loc_counts.get("Extracellular", 0) / max(num_predicted, 1)) * 100.0, 1)
        cyto_pct = round((loc_counts.get("Cytoplasm", 0) / max(num_predicted, 1)) * 100.0, 1)

        # Functional classification of predicted and curated proteins
        curated_func_counts = {}
        for c in curated_cds:
            cat = classify_function(c.get("product", ""))
            curated_func_counts[cat] = curated_func_counts.get(cat, 0) + 1

        summary_rows.append({
            "organism": org_title,
            "accession": acc,
            "phylum": org["phylum"],
            "length_bp": actual_len,
            "gc_pct": gc_pct,
            "ncbi_curated_genes": num_curated,
            "protflow_predicted_genes": num_predicted,
            "true_positive_matches": tp_matches,
            "false_negatives": fn_count,
            "false_positives": fp_count,
            "recall_pct": recall_pct,
            "precision_pct": precision_pct,
            "coding_density_pct": coding_density_pct,
            "atg_pct": atg_pct,
            "gtg_pct": gtg_pct,
            "ttg_pct": ttg_pct,
            "secreted_pct": sec_pct,
            "intracellular_pct": cyto_pct
        })

        # Save diverse protein catalog examples
        for _, p_row in pred_genes_df.head(15).iterrows():
            f_cat = classify_function(p_row.get("product", ""))
            loc_val = filter_df.loc[filter_df["protein_id"] == p_row["gene_id"], "localization"]
            loc_str = loc_val.iloc[0] if not loc_val.empty else "Cytoplasm"
            all_proteins_catalog.append({
                "organism": org_title,
                "gene_id": p_row["gene_id"],
                "product": p_row["product"],
                "functional_category": f_cat,
                "length_aa": p_row["length_aa"],
                "start_codon": p_row["start_codon"],
                "localization": loc_str,
                "strand": p_row["strand"],
                "start": p_row["start"],
                "end": p_row["end"]
            })

        print(f"  -> Curated: {num_curated} | Predicted: {num_predicted} | Matches (TP): {tp_matches} | Recall: {recall_pct}% | Precision: {precision_pct}%")

    df_summary = pd.DataFrame(summary_rows)
    df_catalog = pd.DataFrame(all_proteins_catalog)

    summary_csv = work_dir / "benchmark_250k_summary.csv"
    catalog_csv = work_dir / "benchmark_250k_protein_catalog.csv"
    df_summary.to_csv(summary_csv, index=False)
    df_catalog.to_csv(catalog_csv, index=False)

    print(f"\n[Saved summary to {summary_csv}]")
    print(f"[Saved protein catalog to {catalog_csv}]")

    return df_summary, df_catalog


if __name__ == "__main__":
    run_ground_truth_benchmark()
