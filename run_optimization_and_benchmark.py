"""
Phase 3 Benchmark & Optimization Runner:
Runs multi-variant pipeline grid sweep, generates empirical calibration tables,
and evaluates node additions against Swiss-Prot ground-truth.
"""

from pathlib import Path
import pandas as pd
from protflow.manifest import Manifest
from protflow.runner import Runner
from protflow.benchmark import (
    load_benchmark_dataframe,
    PipelineOptimizer,
    VariantConfig,
    NodeEvaluator
)

def main():
    print("=" * 80)
    print("ProtFlow: Pipeline Optimization, Calibration & Benchmarking (Aşama 3)")
    print("=" * 80)

    # 1. Dataset
    bench_df = load_benchmark_dataframe()
    print(f"[*] Swiss-Prot Kıyas Seti Yüklendi: {len(bench_df)} protein")
    print(f"    - Aktif Enzimler (EC 1-6): {len(bench_df[bench_df['is_active_enzyme']==True])}")
    print(f"    - Negatif Kontroller: {len(bench_df[bench_df['is_active_enzyme']==False])}")

    # 2. Pipeline Sweep across Variants
    print("\n[*] Hat Varyantları Kıyaslama Taraması Başlatılıyor (Grid Sweep)...")
    optimizer = PipelineOptimizer(Path("bacterial_reference_pipeline.yaml"))

    variants = [
        VariantConfig(name="V1_P2Rank_Balanced", pocket_tool="p2rank", min_plddt=70.0),
        VariantConfig(name="V2_fpocket_Balanced", pocket_tool="fpocket", min_plddt=70.0),
        VariantConfig(name="V3_P2Rank_Strict_pLDDT", pocket_tool="p2rank", min_plddt=80.0),
        VariantConfig(name="V4_P2Rank_Pocket_Heavy", pocket_tool="p2rank", min_plddt=60.0, weight_pocket=0.50, weight_quality=0.20),
    ]

    best, summary_df, cal_table = optimizer.sweep(variants)

    print("\n" + "=" * 80)
    print("HAT VARYANTLARI KARŞILAŞTIRMA TABLOSU (Pareto Ranking):")
    print("=" * 80)
    print(summary_df.to_string(index=False))

    print(f"\n[+] En Başarılı Hat: '{best.variant.name}' (F1={best.f1_score}, Precision={best.precision})")

    # 3. Calibration Table
    print("\n" + "=" * 80)
    print("EMPİRİK KALİBRASYON TABLOSU (Enzim Sınıflarına Göre Güvenilirlik):")
    print("=" * 80)
    print(cal_table.metrics_summary.to_string(index=False))

    # 4. Automated Node Evaluator
    print("\n" + "=" * 80)
    print("OTOMATİK DÜĞÜM KIYASLAMA (Yeni Araç vs Mevcut Standart):")
    print("=" * 80)
    evaluator = NodeEvaluator()
    m_p2rank = Manifest.from_yaml("manifests/05_pocket_p2rank.yaml")
    m_fpocket = Manifest.from_yaml("manifests/05_pocket_fpocket.yaml")

    # Evaluate on a benchmark structure
    from protflow.types import Structure
    import gemmi
    cif_file = Path("tests/fixtures/test_prot.cif")
    if not cif_file.exists():
        cif_file.parent.mkdir(parents=True, exist_ok=True)
        st = gemmi.Structure()
        m = gemmi.Model("1")
        ch = gemmi.Chain("A")
        r = gemmi.Residue()
        r.name = "ALA"
        r.seqid = gemmi.SeqId("1")
        a = gemmi.Atom()
        a.name = "CA"
        a.element = gemmi.Element("C")
        a.pos = gemmi.Position(1.0, 2.0, 3.0)
        r.add_atom(a)
        ch.add_residue(r)
        m.add_chain(ch)
        st.add_model(m)
        st.make_mmcif_document().write_file(str(cif_file))

    test_struct = Structure(protein_id="BENCH_EC4_TNAA", file_path=str(cif_file))
    report = evaluator.compare_nodes(
        candidate_manifest=m_p2rank,
        baseline_manifest=m_fpocket,
        test_inputs_list=[{"structure": test_struct}]
    )

    print(f"Aday Araç: {report.candidate_name} | Mevcut Araç: {report.baseline_name}")
    print(f"Hızlanma (Speedup): {report.speedup}x (Aday: {report.candidate_avg_ms}ms, Mevcut: {report.baseline_avg_ms}ms)")
    print(f"Ortalama Skor: Aday={report.mean_candidate_score} vs Mevcut={report.mean_baseline_score}")
    print(f"Karar: {report.recommendation}")
    print("=" * 80)

if __name__ == "__main__":
    main()
