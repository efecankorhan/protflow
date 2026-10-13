from protflow.adapters.genomics import fetch_ncbi_genome
from protflow.pipeline import Pipeline, PipelineStep
from protflow.runner import Runner
from protflow.types import Sequence

def inspect_organism(acc, name):
    print(f"\n=================== {name} ({acc}) ===================")
    res, _ = fetch_ncbi_genome(acc, max_bp=250000)
    seq = Sequence(protein_id=acc.replace(".", "_"), sequence=res["sequence"])
    runner = Runner(enable_cache=False)
    steps = [
        PipelineStep(name="gp", manifest_path="manifests/01_augustus.yaml", inputs={"genome": "genome"}, parameters={"is_prokaryote": False})
    ]
    p = Pipeline(name=f"test_{acc}", steps=steps, runner=runner)
    state = p.run({"genome": seq})
    pred = state["genes"].df
    curated = [c for c in res["cds_records"] if 0 < c.get("end", 0) <= 250000]

    print(f"Curated CDS: {len(curated)} | Pred Genes: {len(pred)}")
    for idx, c in enumerate(curated[:5], 1):
        gname = c.get("gene_symbol") or f"gene_{idx}"
        c_st, c_en = c["start"], c["end"]
        exons = c.get("exons", [])
        print(f"\n[{idx}] Curated {gname} [{c_st}..{c_en}] strand={c['strand']} ({len(exons)} exons: {exons})")
        ovs = pred[(pred["strand"] == c["strand"]) & (pred["start"] <= c_en) & (pred["end"] >= c_st)]
        if ovs.empty:
            print("   -> NO OVERLAP")
        else:
            for _, ov in ovs.iterrows():
                print(f"   -> Pred {ov['gene_id']} [{ov['start']}..{ov['end']}] n_exons={ov.get('num_exons')} d_start={ov['start']-c_st} d_end={ov['end']-c_en}")

if __name__ == "__main__":
    inspect_organism("NC_004354.4", "Drosophila melanogaster")
    inspect_organism("NC_007112.7", "Danio rerio")
