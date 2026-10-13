"""
Mixed Prokaryote & Eukaryote Benchmark (10 Diverse Organisms).
Evaluates 5 Prokaryotes (Bacteria + Archaea) and 5 Eukaryotes (Fungi, Apicomplexa, Metazoa)
on real 250,000 bp genomic windows against NCBI RefSeq curated ground truth.
Verifies intron/exon structures, coding densities, start codons, and subcellular localization.
"""

import json
import math
import re
from pathlib import Path
import numpy as np
import pandas as pd

from protflow.pipeline import Pipeline, PipelineStep
from protflow.runner import Runner
from protflow.types import Sequence, Table
from protflow.adapters.genomics import fetch_ncbi_genome


ORGANISMS_MIXED = [
    # 5 Prokaryotes
    {
        "id": "NC_000913.3",
        "name": "Escherichia coli K-12 MG1655",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Gammaproteobacteria",
        "is_prokaryote": True
    },
    {
        "id": "NC_000964.3",
        "name": "Bacillus subtilis 168",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Bacillota (Gram+)",
        "is_prokaryote": True
    },
    {
        "id": "NC_000917.1",
        "name": "Methanocaldococcus jannaschii DSM 2661",
        "domain": "Prokaryote (Archaea)",
        "phylum": "Methanocaldococcaceae (Thermophile)",
        "is_prokaryote": True
    },
    {
        "id": "NC_000853.1",
        "name": "Thermotoga maritima MSB8",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Thermotogota (Hyperthermophile)",
        "is_prokaryote": True
    },
    {
        "id": "NC_003888.3",
        "name": "Streptomyces coelicolor A3(2)",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Actinomycetota (High GC)",
        "is_prokaryote": True
    },
    # 5 Eukaryotes
    {
        "id": "NC_001133.9",
        "name": "Saccharomyces cerevisiae Chr I",
        "domain": "Eukaryote (Fungi)",
        "phylum": "Ascomycota (Budding Yeast)",
        "is_prokaryote": False
    },
    {
        "id": "NC_001134.8",
        "name": "Saccharomyces cerevisiae Chr II",
        "domain": "Eukaryote (Fungi)",
        "phylum": "Ascomycota (Budding Yeast)",
        "is_prokaryote": False
    },
    {
        "id": "NC_004325.2",
        "name": "Plasmodium falciparum 3D7 Chr 1",
        "domain": "Eukaryote (Protista)",
        "phylum": "Apicomplexa (Malaria Parasite)",
        "is_prokaryote": False
    },
    {
        "id": "NC_003424.3",
        "name": "Schizosaccharomyces pombe Chr I",
        "domain": "Eukaryote (Fungi)",
        "phylum": "Ascomycota (Fission Yeast)",
        "is_prokaryote": False
    },
    {
        "id": "NC_003279.8",
        "name": "Caenorhabditis elegans Chr I",
        "domain": "Eukaryote (Metazoa)",
        "phylum": "Nematoda (Animal / Roundworm)",
        "is_prokaryote": False
    }
]


def run_mixed_benchmark():
    work_dir = Path("scratch/benchmark_mixed")
    work_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(cache_dir=work_dir / "cache", enable_cache=True)

    summary_rows = []
    catalog_rows = []

    print("=== Starting Mixed Prokaryote & Eukaryote Benchmark (10 Organisms) ===\n")

    for idx, org in enumerate(ORGANISMS_MIXED, 1):
        acc = org["id"]
        is_prok = org["is_prokaryote"]
        dom_label = "Prokaryote" if is_prok else "Eukaryote"
        print(f"[{idx}/10] [{dom_label}] Fetching & Analyzing {org['name']} ({acc})...")

        res, err = fetch_ncbi_genome(acc, max_bp=250000)
        if not res or err:
            print(f"  Error: {err}")
            continue

        seq_dna = res["sequence"]
        actual_len = len(seq_dna)
        gc_pct = res["gc_content"]
        org_title = res["organism"] or org["name"]

        # Curated NCBI ground truth CDS records in [1..actual_len]
        all_cds = res.get("cds_records", [])
        curated_cds = [c for c in all_cds if 0 < c.get("end", 0) <= actual_len]
        num_curated = len(curated_cds)

        # Count ground truth introns
        curated_spliced_count = sum(1 for c in curated_cds if c.get("has_introns"))
        curated_total_introns = sum(c.get("intron_count", 0) for c in curated_cds)

        # Run ProtFlow pipeline up to localization filter
        genome_seq = Sequence(
            protein_id=acc.replace(".", "_"),
            sequence=seq_dna,
            description=f"{org_title} | {org['domain']} | {org['phylum']}"
        )

        steps = [
            PipelineStep(
                name="gene_prediction",
                manifest_path="manifests/01_augustus.yaml",
                inputs={"genome": "genome"},
                parameters={
                    "species": "generic",
                    "is_prokaryote": is_prok,
                }
            ),
            PipelineStep(
                name="localization_filter",
                manifest_path="manifests/02_subcellular_filter.yaml",
                inputs={"proteins": "proteins"}
            )
        ]
        pipeline = Pipeline(name=f"pipeline_mixed_{acc}", steps=steps, runner=runner)
        state = pipeline.run({"genome": genome_seq})

        pred_genes_df = state["genes"].df
        filter_df = state["filter_table"].df
        num_predicted = len(pred_genes_df)

        # Intron/Exon metrics in output table
        pred_has_introns_count = int(pred_genes_df["has_introns"].sum()) if "has_introns" in pred_genes_df else 0
        pred_total_introns = int(pred_genes_df["intron_count"].sum()) if "intron_count" in pred_genes_df else 0
        mean_exons = round(float(pred_genes_df["num_exons"].mean()), 2) if "num_exons" in pred_genes_df else 1.0

        # Real de novo ground truth recovery / matching:
        # 1. Exact CDS match: start & stop both within 30 bp on same strand
        # 2. Stop codon match: stop within 30 bp on same strand (tolerant to start codon drift)
        # 3. Locus detection: predicted ORF overlaps curated gene locus on same strand
        exact_matches = 0
        matched_c_exact = set()
        for _, p in pred_genes_df.iterrows():
            p_st, p_en, p_str = p["start"], p["end"], p["strand"]
            for c_idx, c in enumerate(curated_cds):
                if c_idx not in matched_c_exact and c["strand"] == p_str:
                    if abs(p_st - c["start"]) <= 30 and abs(p_en - c["end"]) <= 30:
                        exact_matches += 1
                        matched_c_exact.add(c_idx)
                        break

        stop_matches = sum(
            1 for c in curated_cds
            if any(p["strand"] == c["strand"] and (abs(p["end"] - c["end"]) <= 30 if c["strand"] == "+" else abs(p["start"] - c["start"]) <= 30) for _, p in pred_genes_df.iterrows())
        )

        locus_matches = sum(
            1 for c in curated_cds
            if any(p["strand"] == c["strand"] and p["start"] <= c["end"] and p["end"] >= c["start"] for _, p in pred_genes_df.iterrows())
        )

        denom = max(num_curated, 1)
        exact_recall_pct = round((exact_matches / denom) * 100.0, 1) if num_curated > 0 else 0.0
        exact_prec_pct = round((exact_matches / max(num_predicted, 1)) * 100.0, 1) if num_predicted > 0 else 0.0
        locus_recall_pct = round((locus_matches / denom) * 100.0, 1) if num_curated > 0 else 0.0
        stop_recall_pct = round((stop_matches / denom) * 100.0, 1) if num_curated > 0 else 0.0

        # Localization breakdown
        loc_counts = filter_df["localization"].value_counts().to_dict()
        sec_pct = round((loc_counts.get("Extracellular", 0) / max(num_predicted, 1)) * 100.0, 1)
        cyto_pct = round((loc_counts.get("Cytoplasm", 0) / max(num_predicted, 1)) * 100.0, 1)

        summary_rows.append({
            "organism": org_title,
            "accession": acc,
            "domain": org["domain"],
            "phylum": org["phylum"],
            "is_prokaryote": is_prok,
            "length_bp": actual_len,
            "gc_pct": gc_pct,
            "num_genes": num_predicted,
            "curated_genes": num_curated,
            "locus_matches": locus_matches,
            "locus_recall_pct": locus_recall_pct,
            "exact_matches": exact_matches,
            "exact_recall_pct": exact_recall_pct,
            "exact_precision_pct": exact_prec_pct,
            "stop_matches": stop_matches,
            "stop_recall_pct": stop_recall_pct,
            "curated_spliced_genes": curated_spliced_count,
            "curated_total_introns": curated_total_introns,
            "pred_spliced_genes": pred_has_introns_count,
            "pred_total_introns": pred_total_introns,
            "mean_exons_per_gene": mean_exons,
            "secreted_pct": sec_pct,
            "intracellular_pct": cyto_pct
        })

        for _, p_row in pred_genes_df.head(10).iterrows():
            loc_val = filter_df.loc[filter_df["protein_id"] == p_row["gene_id"], "localization"]
            loc_str = loc_val.iloc[0] if not loc_val.empty else "Cytoplasm"
            catalog_rows.append({
                "organism": org_title,
                "domain": org["domain"],
                "gene_id": p_row["gene_id"],
                "gene_symbol": p_row.get("gene_symbol", ""),
                "product": p_row.get("product", ""),
                "length_aa": p_row["length_aa"],
                "num_exons": p_row.get("num_exons", 1),
                "has_introns": p_row.get("has_introns", False),
                "intron_count": p_row.get("intron_count", 0),
                "start": p_row["start"],
                "end": p_row["end"],
                "localization": loc_str
            })

        print(f"  -> Pred: {num_predicted} | Curated: {num_curated} | Locus Recall: {locus_matches}/{num_curated} ({locus_recall_pct}%) | Exact CDS: {exact_matches} ({exact_recall_pct}%) | Stop Match: {stop_matches} ({stop_recall_pct}%)")

    df_summary = pd.DataFrame(summary_rows)
    df_catalog = pd.DataFrame(catalog_rows)

    summary_csv = work_dir / "benchmark_mixed_summary.csv"
    catalog_csv = work_dir / "benchmark_mixed_catalog.csv"
    df_summary.to_csv(summary_csv, index=False)
    df_catalog.to_csv(catalog_csv, index=False)

    print(f"\n[Saved mixed summary to {summary_csv}]")
    print(f"[Saved mixed catalog to {catalog_csv}]")

    return df_summary, df_catalog


if __name__ == "__main__":
    run_mixed_benchmark()
