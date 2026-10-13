"""
Comprehensive 400,000 bp Genome-Scale Benchmark on 10 Under-Studied / Extremophile Bacteria.
Validates de novo gene prediction and localization against NCBI RefSeq curated ground truth.
Assesses GC-stratified start codon behavior across High GC (>62%), Moderate GC (45-62%), and Low GC (<45%).
"""

import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

from protflow.pipeline import Pipeline, PipelineStep
from protflow.runner import Runner
from protflow.types import Sequence, Table
from protflow.adapters.genomics import fetch_ncbi_genome


ORGANISMS_400K = [
    {
        "id": "NC_005823.1",
        "name": "Deinococcus geothermalis DSM 11300",
        "phylum": "Deinococcota",
        "ecological_niche": "Hot springs / Radiation resistant",
        "gc_regime": "High GC"
    },
    {
        "id": "NC_009440.1",
        "name": "Caldicellulosiruptor bescii DSM 6725",
        "phylum": "Bacillota",
        "ecological_niche": "Hyperthermophile / Cellulolytic anaerobe",
        "gc_regime": "Low GC"
    },
    {
        "id": "NC_002935.2",
        "name": "Chlorobium tepidum TLS",
        "phylum": "Chlorobiota",
        "ecological_niche": "Green sulfur anoxygenic phototroph",
        "gc_regime": "Moderate GC"
    },
    {
        "id": "NC_003030.1",
        "name": "Clostridium acetobutylicum ATCC 824",
        "phylum": "Bacillota",
        "ecological_niche": "Solventogenic obligate anaerobe",
        "gc_regime": "Low GC"
    },
    {
        "id": "NC_003272.1",
        "name": "Nostoc sp. PCC 7120",
        "phylum": "Cyanobacteriota",
        "ecological_niche": "Filamentous nitrogen-fixing phototroph",
        "gc_regime": "Low/Moderate GC"
    },
    {
        "id": "NC_004342.2",
        "name": "Porphyromonas gingivalis W83",
        "phylum": "Bacteroidota",
        "ecological_niche": "Anaerobic asaccharolytic mucosal pathogen",
        "gc_regime": "Moderate GC"
    },
    {
        "id": "NC_004757.1",
        "name": "Rhodopseudomonas palustris CGA009",
        "phylum": "Pseudomonadota (Alphaproteobacteria)",
        "ecological_niche": "Purple non-sulfur phototroph / Aromatic degrader",
        "gc_regime": "High GC"
    },
    {
        "id": "NC_006840.2",
        "name": "Symbiobacterium thermophilum IAM 14863",
        "phylum": "Bacillota",
        "ecological_niche": "Obligate commensal thermophile",
        "gc_regime": "High GC"
    },
    {
        "id": "NC_007164.1",
        "name": "Campylobacter lari RM2100",
        "phylum": "Campylobacterota",
        "ecological_niche": "Avian coastal thermotolerant pathogen",
        "gc_regime": "Low GC"
    },
    {
        "id": "NC_007297.1",
        "name": "Acidiphilium cryptum JF-5",
        "phylum": "Pseudomonadota (Alphaproteobacteria)",
        "ecological_niche": "Acid mine drainage / Iron-reducing acidophile",
        "gc_regime": "High GC"
    }
]


def classify_product(product_name: str) -> str:
    p = (product_name or "").lower()
    if any(k in p for k in ("kinase", "synthase", "synthetase", "reductase", "dehydrogenase", "hydrolase", "peptidase", "transferase", "isomerase", "polymerase", "ligase", "protease", "nuclease")):
        return "Catalytic Enzyme"
    if any(k in p for k in ("transporter", "permease", "efflux", "symporter", "antiporter", "abc transporter", "channel", "porin")):
        return "Transporter / Channel"
    if any(k in p for k in ("regulator", "repressor", "activator", "transcription factor", "response regulator", "sigma factor")):
        return "Transcriptional Regulator"
    if any(k in p for k in ("ribosomal", "rrna", "trna", "ribosome")):
        return "Ribosomal / Translation"
    if any(k in p for k in ("chaperone", "heat shock", "dna repair", "recombinase", "topoisomerase", "helicase", "cas", "crispr")):
        return "Chaperone / DNA Repair / Defense"
    if any(k in p for k in ("membrane", "envelope", "outer membrane", "lipoprotein")):
        return "Membrane / Structural"
    if any(k in p for k in ("toxin", "adhesin", "invasion", "virulence", "antigen", "pilin", "flagell", "nitrogenase", "photosynth", "rubisco", "cellul")):
        return "Specialized Niche Function"
    return "Hypothetical / Other"


def run_400k_benchmark():
    work_dir = Path("scratch/benchmark_400k")
    work_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(cache_dir=work_dir / "cache", enable_cache=True)

    summary_rows = []
    catalog_rows = []

    print("=== Starting 400,000 bp Benchmark on 10 Under-Studied Bacterial Genomes ===\n")

    for idx, org in enumerate(ORGANISMS_400K, 1):
        acc = org["id"]
        print(f"[{idx}/10] Fetching 400k bp: {org['name']} ({acc})...")
        res, err = fetch_ncbi_genome(acc, max_bp=400000)
        if not res or err:
            print(f"  Error: {err}")
            continue

        seq_dna = res["sequence"]
        actual_len = len(seq_dna)
        gc_pct = res["gc_content"]
        org_title = res["organism"] or org["name"]

        # Curated ground truth CDS records in [1..400000]
        all_cds = res.get("cds_records", [])
        curated_cds = [c for c in all_cds if 0 < c.get("end", 0) <= actual_len]
        num_curated = len(curated_cds)

        # Run ProtFlow pipeline up to localization filter
        genome_seq = Sequence(
            protein_id=acc.replace(".", "_"),
            sequence=seq_dna,
            description=f"{org_title} | {org['phylum']} | {org['ecological_niche']}"
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
        pipeline = Pipeline(name=f"pipeline_400k_{acc}", steps=steps, runner=runner)
        state = pipeline.run({"genome": genome_seq})

        pred_genes_df = state["genes"].df
        filter_df = state["filter_table"].df
        num_predicted = len(pred_genes_df)

        # Matching metrics:
        # 1. Exact start/stop match (within 30 bp)
        # 2. Stop codon match (same stop within 6 bp, measuring full CDS locus recovery)
        exact_tp = 0
        locus_tp = 0
        matched_c_exact = set()
        matched_c_locus = set()

        for p_idx, p in pred_genes_df.iterrows():
            p_st, p_en, p_str = p["start"], p["end"], p["strand"]
            p_stop = p_en if p_str == "+" else p_st

            for c_idx, c in enumerate(curated_cds):
                if c["strand"] == p_str:
                    c_stop = c["end"] if c["strand"] == "+" else c["start"]
                    # Check exact match
                    if c_idx not in matched_c_exact:
                        if abs(p_st - c["start"]) <= 30 and abs(p_en - c["end"]) <= 30:
                            exact_tp += 1
                            matched_c_exact.add(c_idx)
                    # Check locus stop match
                    if c_idx not in matched_c_locus:
                        if abs(p_stop - c_stop) <= 6:
                            locus_tp += 1
                            matched_c_locus.add(c_idx)

        exact_recall_pct = round((exact_tp / max(num_curated, 1)) * 100.0, 1)
        exact_precision_pct = round((exact_tp / max(num_predicted, 1)) * 100.0, 1)
        locus_recall_pct = round((locus_tp / max(num_curated, 1)) * 100.0, 1)

        # Coding density
        coding_density_pct = round((pred_genes_df["length_aa"].sum() * 3 / actual_len) * 100.0, 1)

        # Start codon distribution
        start_counts = pred_genes_df["start_codon"].value_counts().to_dict()
        atg_pct = round((start_counts.get("ATG", 0) / max(num_predicted, 1)) * 100.0, 1)
        gtg_pct = round((start_counts.get("GTG", 0) / max(num_predicted, 1)) * 100.0, 1)
        ttg_pct = round((start_counts.get("TTG", 0) / max(num_predicted, 1)) * 100.0, 1)

        # Subcellular localization
        loc_counts = filter_df["localization"].value_counts().to_dict()
        sec_pct = round((loc_counts.get("Extracellular", 0) / max(num_predicted, 1)) * 100.0, 1)
        cyto_pct = round((loc_counts.get("Cytoplasm", 0) / max(num_predicted, 1)) * 100.0, 1)

        summary_rows.append({
            "organism": org_title,
            "accession": acc,
            "phylum": org["phylum"],
            "ecological_niche": org["ecological_niche"],
            "gc_regime": org["gc_regime"],
            "actual_len_bp": actual_len,
            "gc_pct": gc_pct,
            "curated_ncbi_genes": num_curated,
            "protflow_predicted_genes": num_predicted,
            "exact_matches_tp": exact_tp,
            "exact_recall_pct": exact_recall_pct,
            "exact_precision_pct": exact_precision_pct,
            "locus_stop_matches_tp": locus_tp,
            "locus_recovery_pct": locus_recall_pct,
            "coding_density_pct": coding_density_pct,
            "start_atg_pct": atg_pct,
            "start_gtg_pct": gtg_pct,
            "start_ttg_pct": ttg_pct,
            "secreted_pct": sec_pct,
            "intracellular_pct": cyto_pct
        })

        # Capture representative proteins with diverse functions
        for _, p_row in pred_genes_df.head(12).iterrows():
            loc_val = filter_df.loc[filter_df["protein_id"] == p_row["gene_id"], "localization"]
            loc_str = loc_val.iloc[0] if not loc_val.empty else "Cytoplasm"
            catalog_rows.append({
                "organism": org_title,
                "niche": org["ecological_niche"],
                "gene_id": p_row["gene_id"],
                "length_aa": p_row["length_aa"],
                "start_codon": p_row["start_codon"],
                "strand": p_row["strand"],
                "start": p_row["start"],
                "end": p_row["end"],
                "localization": loc_str
            })

        print(f"  -> Curated: {num_curated} | Predicted: {num_predicted} | Exact Recall: {exact_recall_pct}% | Locus Recall: {locus_recall_pct}% | Start ATG/GTG/TTG: {atg_pct}%/{gtg_pct}%/{ttg_pct}%")

    df_summary = pd.DataFrame(summary_rows)
    df_catalog = pd.DataFrame(catalog_rows)

    summary_csv = work_dir / "benchmark_400k_summary.csv"
    catalog_csv = work_dir / "benchmark_400k_protein_catalog.csv"
    df_summary.to_csv(summary_csv, index=False)
    df_catalog.to_csv(catalog_csv, index=False)

    print(f"\n[Saved summary to {summary_csv}]")
    print(f"[Saved catalog to {catalog_csv}]")

    return df_summary, df_catalog


if __name__ == "__main__":
    run_400k_benchmark()
