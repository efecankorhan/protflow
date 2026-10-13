from protflow.adapters.genomics import fetch_ncbi_genome, parse_augustus
from protflow.types import Sequence

def check_unmatched(acc, is_prok):
    res, _ = fetch_ncbi_genome(acc, max_bp=250000)
    seq = Sequence(protein_id=acc.replace(".", "_"), sequence=res["sequence"])
    out = parse_augustus(None, {"inputs": {"genome": seq}, "parameters": {"is_prokaryote": is_prok}})
    pred = out["genes"].df
    curated = [c for c in res["cds_records"] if 0 < c.get("end", 0) <= 250000]

    unmatched = []
    for c in curated:
        found = any(p["strand"] == c["strand"] and abs(p["start"] - c["start"]) <= 30 and abs(p["end"] - c["end"]) <= 30 for _, p in pred.iterrows())
        if not found:
            unmatched.append(c)

    print(f"=== {acc} (is_prok={is_prok}) ===")
    print(f"Curated: {len(curated)} | Pred: {len(pred)} | Unmatched: {len(unmatched)}")
    
    no_overlap = 0
    start_diff_only = 0
    stop_diff_only = 0
    both_diff = 0
    short_genes = 0

    for u in unmatched:
        ovs = pred[(pred["strand"] == u["strand"]) & (pred["start"] <= u["end"]) & (pred["end"] >= u["start"])]
        length = u["end"] - u["start"] + 1
        if length < 90:
            short_genes += 1
        if ovs.empty:
            no_overlap += 1
        else:
            best_ov = ovs.iloc[0]
            d_st = abs(best_ov["start"] - u["start"])
            d_en = abs(best_ov["end"] - u["end"])
            if d_en <= 30 and d_st > 30:
                start_diff_only += 1
            elif d_st <= 30 and d_en > 30:
                stop_diff_only += 1
            else:
                both_diff += 1

    print(f"Analysis of {len(unmatched)} unmatched:")
    print(f"  - No overlap at all: {no_overlap} (likely small peptides / ignored by threshold)")
    print(f"  - Length < 90 bp (<30 aa): {short_genes}")
    print(f"  - Stop codon matched (<=30 bp) but Start codon shifted (>30 bp): {start_diff_only}")
    print(f"  - Start codon matched but Stop shifted: {stop_diff_only}")
    print(f"  - Both shifted: {both_diff}")

    print("\nSample 5 unmatched details:")
    for u in unmatched[:5]:
        ovs = pred[(pred["strand"] == u["strand"]) & (pred["start"] <= u["end"]) & (pred["end"] >= u["start"])]
        gname = u.get("gene_symbol") or u.get("gene_id")
        prod = u.get("product")
        l = u["end"] - u["start"] + 1
        if not ovs.empty:
            ov = ovs.iloc[0]
            print(f"  Curated {gname} ({prod}) [{u['start']}..{u['end']}] len={l} vs Pred [{ov['start']}..{ov['end']}] d_start={ov['start']-u['start']} d_stop={ov['end']-u['end']}")
        else:
            print(f"  Curated {gname} ({prod}) [{u['start']}..{u['end']}] len={l} -> NO OVERLAP")

if __name__ == "__main__":
    check_unmatched("NC_000913.3", True)   # E. coli
    check_unmatched("NC_003888.3", True)   # S. coelicolor
    check_unmatched("NC_001133.9", False)  # S. cerevisiae
    check_unmatched("NC_003279.8", False)  # C. elegans
