"""
Reference Pipeline Runner:
Executes Genome -> Prodigal -> Subcellular -> AFDB -> pLDDT filter -> Pocket -> SASA -> Foldseek -> Prioritization.
Outputs formatted candidate priority table.
"""

from pathlib import Path
import json
from protflow.pipeline import Pipeline
from protflow.runner import Runner
from protflow.types import Sequence

def main():
    print("=" * 70)
    print("ProtFlow: Bacterial Reference Pipeline Execution (Aşama 2)")
    print("=" * 70)

    # 1. Load pipeline configuration
    pipeline_yaml = Path("bacterial_reference_pipeline.yaml")
    runner = Runner(cache_dir=".protflow_cache", enable_cache=True)
    pipeline = Pipeline.from_yaml(pipeline_yaml, runner=runner)
    print(f"[*] Pipeline loaded: {pipeline.name} ({len(pipeline.steps)} adim)")

    # 2. Input Genome (E. coli K-12 fragment coding for enzyme)
    genome_dna = (
        "ATGGCTGGGCTGGTGGCCATCCTGGAGAAATTCGGCGACGAGCTGCGTTGGGCGCCGAAA"
        "ACCTGCCACTATGGCGCGGTGATCCTGTCCTAA"
    )
    genome_seq = Sequence(
        protein_id="ECOLI_MG1655_EXP1",
        sequence=genome_dna,
        description="E. coli K-12 representative enzyme locus"
    )
    print(f"[*] Input Genome: {genome_seq.protein_id} ({len(genome_dna)} bp)")

    # 3. Run Pipeline
    print("\n[*] Calistirma baslatiliyor...")
    final_state = pipeline.run({"genome": genome_seq})

    # 4. Results
    print("\n" + "=" * 70)
    print("SONUC: Protein Oncelik Tablosu (Candidate Prioritization)")
    print("=" * 70)
    priority_table = final_state["priority_table"].df
    print(priority_table.to_string(index=False))

    print("\n" + "-" * 70)
    print("Tespit Edilen Baglanma Cebleri (P2Rank / Cavity):")
    print("-" * 70)
    pockets_table = final_state["pockets"].df
    print(pockets_table.to_string(index=False))

    print("\n" + "-" * 70)
    print("Provenance (AlphaFold Yapisi & Guven Kaydi):")
    print("-" * 70)
    struct_prov = final_state["structure"].provenance
    print(f"Tool: {struct_prov.tool_name} v{struct_prov.tool_version}")
    print(f"Remote: {struct_prov.is_remote} | Endpoint: {struct_prov.endpoint}")
    print(f"Execution Time: {struct_prov.execution_time}")
    print("=" * 70)

if __name__ == "__main__":
    main()
