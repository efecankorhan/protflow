"""
10 New Diverse Organisms Benchmark (250,000 bp windows).
Includes 1 large DNA Virus, 3 Bacteria, 2 Archaea, 1 Fungus, 1 Plant, 1 Insect, 1 Vertebrate Animal.
Strictly Zero Data Leakage: Pipeline receives ONLY raw genomic DNA sequence.
Evaluates de novo gene prediction against NCBI RefSeq ground truth.
"""

from pathlib import Path
import numpy as np
import pandas as pd

from protflow.pipeline import Pipeline, PipelineStep
from protflow.runner import Runner
from protflow.types import Sequence
from protflow.adapters.genomics import fetch_ncbi_genome


NEW_ORGANISMS = [
    # 1. Large DNA Virus
    {
        "id": "NC_006273.2",
        "name": "Human cytomegalovirus / HCMV",
        "domain": "Virus (Betaherpesvirinae)",
        "phylum": "Duplodnaviria (Large dsDNA Virus)",
        "is_prokaryote": True  # Continuous viral CDS models
    },
    # 2. Bacteria (Gammaproteobacteria)
    {
        "id": "NC_002516.2",
        "name": "Pseudomonas aeruginosa PAO1",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Pseudomonadota (Opportunistic Pathogen)",
        "is_prokaryote": True
    },
    # 3. Bacteria (High-GC Actinobacteria)
    {
        "id": "NC_000962.3",
        "name": "Mycobacterium tuberculosis H37Rv",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Actinomycetota (High GC, Slow Growing)",
        "is_prokaryote": True
    },
    # 4. Bacteria (Low-GC Microaerophile)
    {
        "id": "NC_000915.1",
        "name": "Helicobacter pylori 26695",
        "domain": "Prokaryote (Bacteria)",
        "phylum": "Campylobacterota (Low GC)",
        "is_prokaryote": True
    },
    # 5. Archaea (Thermoacidophile)
    {
        "id": "NC_007181.1",
        "name": "Sulfolobus acidocaldarius DSM 639",
        "domain": "Prokaryote (Archaea)",
        "phylum": "Thermoproteota (Acidophile / Thermophile)",
        "is_prokaryote": True
    },
    # 6. Archaea (Extreme Halophile, High GC)
    {
        "id": "NC_013967.1",
        "name": "Haloferax volcanii DS2",
        "domain": "Prokaryote (Archaea)",
        "phylum": "Halobacteriota (High Salt / High GC)",
        "is_prokaryote": True
    },
    # 7. Eukaryote (Fungi - Opportunistic Yeast)
    {
        "id": "NC_032089.1",
        "name": "Candida albicans SC5314 Chr 1",
        "domain": "Eukaryote (Fungi)",
        "phylum": "Ascomycota (Diploid Yeast)",
        "is_prokaryote": False
    },
    # 8. Eukaryote (Plant - Model Dicot)
    {
        "id": "NC_003070.9",
        "name": "Arabidopsis thaliana Chr 1",
        "domain": "Eukaryote (Plantae)",
        "phylum": "Brassicaceae (Model Flowering Plant)",
        "is_prokaryote": False
    },
    # 9. Eukaryote (Animal - Invertebrate Model)
    {
        "id": "NC_004354.4",
        "name": "Drosophila melanogaster 2L",
        "domain": "Eukaryote (Metazoa)",
        "phylum": "Arthropoda (Fruit Fly)",
        "is_prokaryote": False
    },
    # 10. Eukaryote (Animal - Vertebrate Fish)
    {
        "id": "NC_007112.7",
        "name": "Danio rerio Chr 1",
        "domain": "Eukaryote (Metazoa)",
        "phylum": "Chordata (Zebrafish / Vertebrate)",
        "is_prokaryote": False
    }
]


def run_benchmark():
    work_dir = Path("scratch/benchmark_new_10")
    work_dir.mkdir(parents=True, exist_ok=True)
    runner = Runner(cache_dir=work_dir / "cache", enable_cache=True)

    summary_rows = []
    catalog_rows = []

    print("=== Starting Benchmark on 10 NEW Organisms (250,000 bp windows) ===")
    print("Strict Zero-Leakage: De novo prediction from raw DNA sequence only.\n")

    for idx, org in enumerate(NEW_ORGANISMS, 1):
        acc = org["id"]
        is_prok = org["is_prokaryote"]
        dom_label = org["domain"]
        print(f"[{idx}/10] [{dom_label}] Fetching & Predicting {org['name']} ({acc})...")

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

        # Distinct biological gene loci (merging alternative splicing isoforms of the same gene)
        seen_syms = set()
        distinct_loci = []
        for c in curated_cds:
            sym = (c.get("gene_symbol") or "").strip()
            if sym and sym in seen_syms:
                continue
            if sym:
                seen_syms.add(sym)
            distinct_loci.append(c)
        num_unique_genes = len(distinct_loci)

        curated_spliced_count = sum(1 for c in curated_cds if c.get("has_introns"))
        curated_total_introns = sum(c.get("intron_count", 0) for c in curated_cds)

        # STRICT ZERO LEAKAGE: Feed ONLY raw DNA into pipeline
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
                    "min_len_aa": 25
                }
            ),
            PipelineStep(
                name="localization_filter",
                manifest_path="manifests/02_subcellular_filter.yaml",
                inputs={"proteins": "proteins"}
            )
        ]
        pipeline = Pipeline(name=f"pipe_new_{acc}", steps=steps, runner=runner)
        state = pipeline.run({"genome": genome_seq})

        pred_genes_df = state["genes"].df
        filter_df = state["filter_table"].df
        num_predicted = len(pred_genes_df)

        pred_has_introns_count = int(pred_genes_df["has_introns"].sum()) if "has_introns" in pred_genes_df else 0
        pred_total_introns = int(pred_genes_df["intron_count"].sum()) if "intron_count" in pred_genes_df else 0
        mean_exons = round(float(pred_genes_df["num_exons"].mean()), 2) if "num_exons" in pred_genes_df else 1.0

        # Ground truth recovery:
        # 1. Exact CDS match: both start and stop within 30 bp on same strand
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

        # 2. Stop codon match (tolerant to alternative start codon selection)
        stop_matches = sum(
            1 for c in curated_cds
            if any(p["strand"] == c["strand"] and (abs(p["end"] - c["end"]) <= 30 if c["strand"] == "+" else abs(p["start"] - c["start"]) <= 30) for _, p in pred_genes_df.iterrows())
        )

        # 3. Locus detection (predicted ORF overlaps curated gene locus on same strand)
        locus_matches = sum(
            1 for c in curated_cds
            if any(p["strand"] == c["strand"] and p["start"] <= c["end"] and p["end"] >= c["start"] for _, p in pred_genes_df.iterrows())
        )

        denom = max(num_curated, 1)
        locus_recall_pct = round((locus_matches / denom) * 100.0, 1) if num_curated > 0 else 0.0
        exact_recall_pct = round((exact_matches / denom) * 100.0, 1) if num_curated > 0 else 0.0
        exact_prec_pct = round((exact_matches / max(num_predicted, 1)) * 100.0, 1) if num_predicted > 0 else 0.0
        stop_recall_pct = round((stop_matches / denom) * 100.0, 1) if num_curated > 0 else 0.0

        # 4. Exon-level recovery (individual coding exon boundary detection)
        all_exons = [e for c in curated_cds for e in c.get("exons", [(c["start"], c["end"])])]
        total_exons = len(all_exons)
        exon_matches = sum(
            1 for e_st, e_en in all_exons
            if any(abs(p["start"] - e_st) <= 30 or abs(p["end"] - e_en) <= 30 for _, p in pred_genes_df.iterrows())
        )
        exon_recall_pct = round((exon_matches / max(total_exons, 1)) * 100.0, 1) if total_exons > 0 else 100.0

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
            "unique_loci": num_unique_genes,
            "total_exons": total_exons,
            "exon_matches": exon_matches,
            "exon_recall_pct": exon_recall_pct,
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

        print(f"  -> Pred: {num_predicted} | Unique Loci: {num_unique_genes} (CDS: {num_curated}) | Locus: {locus_matches}/{num_curated} ({locus_recall_pct}%) | Exon: {exon_matches}/{total_exons} ({exon_recall_pct}%) | Exact CDS: {exact_matches} ({exact_recall_pct}%) | Stop: {stop_matches} ({stop_recall_pct}%)")

    df_summary = pd.DataFrame(summary_rows)
    df_catalog = pd.DataFrame(catalog_rows)

    summary_csv = work_dir / "benchmark_new_10_summary.csv"
    catalog_csv = work_dir / "benchmark_new_10_catalog.csv"
    df_summary.to_csv(summary_csv, index=False)
    df_catalog.to_csv(catalog_csv, index=False)

    print(f"\n[Saved summary to {summary_csv}]")
    print(f"[Saved catalog to {catalog_csv}]")
    return df_summary, df_catalog


if __name__ == "__main__":
    run_benchmark()
