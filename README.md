# ProtFlow: Zero-Friction Bio-Pipeline Framework (Aşama 1)

ProtFlow, araç eklemeyi (local binary, container, remote web service) neredeyse sıfır sürtünmeyle gerçekleştiren ve tipli sözleşmelerle yöneten hesaplamalı biyoloji omurgasıdır.

## 1. Yedi Çekirdek Veri Tipi (`protflow.types`)

Tüm veri tipleri `protein_id` ortak anahtarı ile bağlanır ve `Provenance` (araç, sürüm, parametreler, tarih, uzak servis bilgisi) taşır:
1. **Sequence** (`FASTA`): Sekans verisi.
2. **Structure** (`mmCIF` / `PDB`): 3B koordinatlar (`gemmi` tabanlı model).
3. **ResidueScore** (`Parquet`): Kalıntı düzeyinde skorlar (`pLDDT`, çözünürlük, hidrofobisite).
4. **Pocket** (`Parquet`): Tahmin edilen bağlanma cepleri (koordinat, hacim, kalıntılar).
5. **Alignment** (`Tabular` / `Parquet`): Sekans ve yapı hizalama sonuçları.
6. **Annotation** (`GFF3` / `Parquet`): Gen ve bölge anotasyonları.
7. **Table** (`Parquet`): Genel metrikler ve birleşik skor tabloları.

## 2. Düğüm Manifesti (`protflow.manifest`)

Her araç bağımsız bir YAML manifestiyle tanımlanır:

```yaml
name: pocket_finder
version: "1.0.0"
inputs:
  structure: Structure
outputs:
  pockets: Pocket
execution:
  type: local # local | container | http
  command: "p2rank predict -f {inputs_structure} -o {output_dir}"
adapter: protflow.adapters.parsers:parse_dummy_pocket_finder
```

## 3. Güvenlik ve Gizlilik Katmanı

Bir düğüm `is_remote: true` veya `execution.type: http` olarak tanımlandığında:
* Konsol ve loglara belirgin `[PRIVACY WARNING]` uyarısı düşer.
* Yayımlanmamış sekansların dış ağlara sızması baştan görünür kılınır.

## 4. Kalıntı Kayması (Residue Drift) Doğrulaması (`protflow.adapters.residue`)

Araçlar arasındaki 0-based/1-based veya PDB/mmCIF `auth_seq_id` uyuşmazlıkları `validate_residue_table` ile doğrulanır:
* Kalıntının yapıda varlığı test edilir.
* Kalıntı amino asit harfi sekansla eşleşmiyorsa `ResidueDriftError` fırlatılır.

## 5. Aşama 2: Uçtan Uca Bakteriyel Referans Hat

Tek bir konfigürasyon dosyasıyla ([`bacterial_reference_pipeline.yaml`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/bacterial_reference_pipeline.yaml)) 8 adımlı referans hat çalıştırılır:
1. **Genom $\rightarrow$ Prodigal:** ORF tespiti ve CDS çevrimi (`Sequence`, `Annotation`).
2. **Subcellular Filter:** Sinyal peptidi / transmembran eleme (`Table`).
3. **AlphaFold DB:** Uzak REST API ile 3B yapı ve kalıntı bazlı pLDDT çekme (`Structure`, `ResidueScore`).
4. **pLDDT Filtresi:** Düşük kaliteli modelleri eleme.
5. **Cep Tespiti (P2Rank / fpocket):** Bağlanma cebi ve hacim tespiti (`Pocket`). Tek satırda P2Rank $\leftrightarrow$ fpocket değişimi!
6. **FreeSASA:** Yüzey alanı ve hidrofobisite oranı hesabı (`Table`).
7. **Foldseek:** Yapısal eşleme ve yenilik (`novelty`) skoru (`Alignment`, `Table`).
8. **Önceliklendirici (Prioritizer):** Çok kriterli aday öncelik skoru (`priority_table`).

### Örnek Çalıştırma
```bash
python run_reference_pipeline.py
```

## 6. Aşama 3: Optimizasyon, Kalibrasyon & Benchmark (`protflow.benchmark`)

1. **Kullanıcı Riski Çözümü (Yapı Tahmini Koruma Katmanı):**
   * AFDB'de bulunmayan yeni sekanslar için sessizce sahte model üretimi engellendi.
   * `allow_synthetic_fallback: false` durumunda hat açıkça [`StructureNotFoundError`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/protflow/adapters/pipeline_tools.py#L162-L165) fırlatır.
   * Yerel de novo modelleme için [`manifests/03_colabfold.yaml`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/manifests/03_colabfold.yaml) manifesti eklendi.
2. **Düğüm Düzeyinde Skor Normalizasyonu ([`normalize.py`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/protflow/normalize.py)):**
   * P2Rank, fpocket, pLDDT ve SASA değerleri ortak $[0.0, 1.0]$ ölçeğine dönüştürülür.
3. **Swiss-Prot Kıyas Seti ([`dataset.py`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/protflow/benchmark/dataset.py)):**
   * 6 EC enzim sınıfı (EC 1-6) ve 3 negatif kontrol (yapısal/düzensiz proteinler).
4. **Hat Varyantı Optimizasyonu ([`optimizer.py`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/protflow/benchmark/optimizer.py)):**
   * Farklı cep bulucuları ve pLDDT eşiklerini Pareto metriğiyle (F1, Precision, Latency) kıyaslar.
5. **Empirik Kalibrasyon Tablosu ([`calibration.py`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/protflow/benchmark/calibration.py)):**
   * Her protein sınıfı için ampirik doğruluk tablosu ve güvenilirlik derecesi (Grade A-D).
6. **Otomatik Düğüm Kıyaslayıcı ([`evaluator.py`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/protflow/benchmark/evaluator.py)):**
   * Sisteme yeni eklenen araçları mevcut standartla otomatik kıyaslayıp öneri raporu basar.

### Optimizasyon ve Kıyas Çalıştırma
```bash
python run_optimization_and_benchmark.py
```

## 7. Aşama 4: ComfyUI Tarzı Tipli Düğüm Editörü & Islak Lab (`protflow.ui`)

1. **Görsel Düğüm Editörü (Studio):**
   * Sıfır npm bağımlılığıyla doğrudan Python üzerinden çalışır (`python run_ui.py`).
   * **Tipli Port Doğrulaması:** Farklı tiplerdeki portlar bağlanamaz (örneğin `Pocket` $\rightarrow$ `Structure` bağlantısı kırmızı uyarı vererek reddedilir).
   * **Çalıştırma Modu Rozetleri:** Her düğümün üzerinde `LOCAL`, `CONTAINER` veya `REMOTE / GPU` açıkça görünür.
   * **Özel Düğüm Sihirbazı:** Arayüz üzerinden YAML düzenleme, örnek girdiyle canlı test etme (`/api/manifest/test`) ve kaydetme.
2. **Kiralık Sunucu / GPU Modülerliği (RFdiffusion, Chai, vb.):**
   * Ağır GPU modelleri motora gömülü değildir; [`manifests/09_rfdiffusion.yaml`](file:///C:/Users/efeca/.gemini/antigravity/scratch/protflow/manifests/09_rfdiffusion.yaml) gibi bağımsız manifestlerle uzak kiralık sunuculara (`HTTP / Container`) yönlendirilir.
3. **Islak Lab Alanin Taraması (Alanine Scanning):**
   * Hattın seçtiği ilk 10 aday enzimin aktif bölgesindeki kalıntılar için otomatik mutasyon paneli (`A:35 GLU->ALA`) üretir.
   * Islak laboratuvar deney sonuçları geri beslenerek kalibrasyon tablosu otomatik güncellenir.

### Web Editörünü Başlatma
```bash
python run_ui.py
# Tarayıcıda açın: http://127.0.0.1:8080
```

## 8. Snakemake & Nextflow İnce Katmanı (`protflow.export`)

Manifestler ve tüm hat doğrudan native pipeline motorlarına çevrilebilir:
* `pipeline.export_snakemake()` $\rightarrow$ `Snakefile`
* `pipeline.export_nextflow()` $\rightarrow$ `main.nf`

## 9. Testleri Çalıştırma

```bash
python -m pytest -v
```
*(26/26 test: Modeller, kalıntı kayması, referans hat, varyant optimizasyonu, kıyas, alanin taraması ve UI API testleri.)*
