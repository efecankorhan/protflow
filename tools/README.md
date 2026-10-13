# ProtFlow Biyoinformatik Araç Deposu (`tools/`)

Bu dizin, ProtFlow boru hatlarının kullandığı harici yerel (local) hesaplamalı biyoloji ve yapısal biyoinformatik araçlarının resmi Git alt modüllerini (`git submodules`) ve çalışma zamanı dosyalarını içerir.

ProtFlow, bu araçların kurulu/derlenmiş olup olmadığını otomatik olarak tespit eder (`auto-detection`). Sistemde veya bu dizinde bulunamazlarsa, hat kesintiye uğramaz; saf Python/Gemmi tabanlı geometrik algoritmalara veya uzak REST servislerine şeffaf bir şekilde düşüş (`fallback`) yapar.

---

## 1. Dahil Edilen Resmi Araç Repoları

| Araç | Konum / Alt Modül | Açıklama & Algoritma | Girdi Formatı | Çıktı Formatı |
| :--- | :--- | :--- | :--- | :--- |
| **fpocket** | `tools/fpocket` | Voronoi mozaiklemesi ve alfa küreleri (alpha spheres) ile protein kavitasyon ve bağlanma cebi tespiti | `.pdb`, `.cif` | `.pqr`, `*_info.txt`, `pockets/` |
| **P2Rank** | `tools/p2rank` | Rastgele Orman (Random Forest) makine öğrenmesi ve küresel yüzey noktaları tabanlı ligand bağlanma cebi tahmini | `.pdb`, `.cif` | `*_predictions.csv`, `*_residues.csv`, PyMOL `.pml` |
| **Foldseek** | `tools/foldseek` | 3B yapısal hizalama ve yapı arama (3Di yapısal alfabe ile TM-align / DALI'den 10.000x daha hızlı) | `.pdb`, `.cif` | `.m8` (BLAST benzeri tabular), `.tsv` |
| **ColabFold** | `tools/colabfold` | MMseqs2 çoklu sekans hizalama (MSA) + AlphaFold2 / AlphaFold-Multimer ile süper hızlı de novo 3B modelleme | `.fasta` | `.pdb`, `.cif`, pLDDT / PAE `.json` |
| **MMseqs2** | `tools/mmseqs2` | Milyarlarca sekans üzerinde ultra hızlı sekans benzerlik araması, profilleme ve kümeleme motoru | `.fasta` | `.m8`, `.tsv`, küme listeleri |
| **MAFFT** | `tools/mafft` | Kazutaka Katoh'un resmi BSD lisanslı çoklu sekans hizalama (MSA) motoru (FFT tabanlı) | `.fasta` | `.fasta`, `.aln` |
| **FAMSA v2** | `tools/famsa` | Deorowicz et al. (Bioinformatics 2021) ultra hızlı, büyük ölçekli çoklu sekans hizalama | `.fasta` | `.fasta`, `.aln` |

---

## 2. Araç Detayları ve Derleme / Çalıştırma Yönergeleri

### 1. `fpocket` (Discngine / Peter Schmidtke)
- **Ne İşe Yarar?** Protein yüzeyindeki kavitasyonları tespit eder. Her cebin hacmini, polarite/hidrofobisite oranını ve ilaca yatkınlık skorunu (*Druggability Score*, 0.0 - 1.0) hesaplar.
- **Kurulum / Derleme:**
  ```bash
  cd tools/fpocket
  # Linux / macOS / WSL:
  make
  # Derlenen çalıştırılabilir dosya: tools/fpocket/bin/fpocket
  ```
- **ProtFlow Entegrasyonu:**
  `protflow.adapters.pockets:parse_fpocket` fonksiyonu sırasıyla `tools/fpocket/bin/fpocket`, `tools/fpocket/fpocket` ve sistem `PATH` değişkenini arar. Bulamazsa Gemmi koordinat kavitasyon algoritmasıyla yaklaşık cep üretir.
- **Girdi:** 3B Koordinat Dosyası (`Structure`, `.pdb` veya `.cif`).
- **Çıktı:** `Pocket` nesnesi (`DataFrame`: `pocket_id`, `score`, `volume`, `center_x/y/z`, `residues`).

---

### 2. `P2Rank` (David Hoksza, Petr Skoda / RDK)
- **Ne İşe Yarar?** Protein yüzeyinde ligand bağlayabilecek aktif cepleri ve amino asit kalıntılarını tahmin eder. PDB yapıları için otomatik PyMOL görselleştirme betiği (`.pml`) üretir.
- **Gereksinimler:** Java JRE 8+ (`tools/jre` veya sistem `JAVA_HOME`).
- **Kurulum / Derleme:**
  ```bash
  cd tools/p2rank
  ./gradlew dist
  # veya hazır ikili paket: tools/p2rank_2.5.1/
  ```
- **ProtFlow Entegrasyonu:**
  `protflow.adapters.pockets:parse_p2rank` fonksiyonu `tools/p2rank_2.5.1` veya `tools/p2rank` altındaki `prank.bat` (Windows) ya da `prank` (Linux/Mac) betiğini çalıştırır.
- **Girdi:** `Structure` (`.pdb`).
- **Çıktı:** `Pocket` nesnesi + Provenance içinde PyMOL `.pml` betiği.

---

### 3. `Foldseek` (Martin Steinegger Lab)
- **Ne İşe Yarar?** 3B protein yapılarını yapısal alfabe (`3Di`) dizilerine dönüştürerek saniyeler içinde PDB ve AlphaFold DB üzerinde benzer yapıları tarar.
- **Kurulum / Derleme:**
  ```bash
  cd tools/foldseek
  mkdir build && cd build
  cmake -DCMAKE_BUILD_TYPE=RELEASE ..
  make -j
  ```
- **ProtFlow Entegrasyonu:**
  Yapısal yenilik (`novelty`) ve PDB homoloji eşleme adımlarında kullanılır (`manifests/07_foldseek.yaml`).
- **Girdi:** Sorgu yapısı (`Structure`, `.pdb` / `.cif`).
- **Çıktı:** `Alignment` ve `Table` (TM-score, RMSD, E-değeri, hizalama örtüşmesi).

---

### 4. `ColabFold` (Sergey Ovchinnikov & Martin Steinegger)
- **Ne İşe Yarar?** Homoloji modellemesi yetersiz kaldığında veya AlphaFold DB'de bulunmayan özgün proteinler için GPU destekli yerel veya bulut AlphaFold2 tahmini yapar.
- **Kurulum:**
  ```bash
  pip install colabfold[all]
  # veya Docker konteyneri: colabfold:latest
  ```
- **ProtFlow Entegrasyonu:**
  `manifests/03_colabfold.yaml` manifesti üzerinden `execution.type: local` veya `container` / `http` olarak çağrılır.
- **Girdi:** Amino asit sekansı (`Sequence`, `.fasta`).
- **Çıktı:** 3B Koordinat Yapısı (`Structure`, `.pdb`) ve kalıntı bazlı pLDDT tablosu (`ResidueScore`).

---

### 5. `MMseqs2` (Johannes Söding Lab)
- **Ne İşe Yarar?** Milyonlarca sekanslık veri tabanlarında (UniProt, Swiss-Prot, NCBI nr) BLAST'a göre 100-1000 kat daha hızlı homoloji ve ortolog arar.
- **Kurulum:**
  `tools/mmseqs/` klasörüne statik ikili dosya indirilebilir veya PATH üzerinden çalıştırılabilir.
- **ProtFlow Entegrasyonu:**
  `protflow.adapters.homology:parse_mmseqs` Swiss-Prot veritabanı eşleştirmesi (`manifests/11_mmseqs_homology.yaml`).
- **Girdi:** Sorgu sekansı (`Sequence`).
- **Çıktı:** `Table` / `Alignment` (Pfam alanları, E-değerleri, yüzdelik benzerlik).

---

## 3. Veri Formatları ve Standartlar

ProtFlow boru hatları boyunca tüm araçlar aşağıdaki tipli veri formatlarını kullanır:

```
[Genom / Nükleotid FASTA]
          │
          ▼ (01_prodigal: pyrodigal)
[Sequence (AA FASTA)] + [Annotation (GFF3 / JSON)]
          │
          ▼ (03_alphafold_db / 03_colabfold)
[Structure (PDB / mmCIF)] + [ResidueScore (pLDDT Parquet / CSV)]
          │
          ▼ (05_pocket_p2rank / 05_pocket_fpocket)
[Pocket (Parquet / CSV)] + [PyMOL PML Script]
          │
          ▼ (06_surface_sasa & 07_foldseek)
[Table / Alignment (Parquet / CSV / JSON)]
```

### Bölgesel Ayar & Excel Uyumluluğu (Türkçe Windows Desteği)
- **Biyoinformatik Standart:** Virgül (`,`) sütun ayırıcı, nokta (`.`) ondalık basamak (RFC-4180). Python, R, AWK, Cut araçları için uygundur.
- **Türkçe Excel Standart:** Noktalı virgül (`;`) sütun ayırıcı, virgül (`,`) ondalık basamak ve UTF-8 BOM (`\ufeff`). Microsoft Excel'de çift tıklandığında sütunlar otomatik bölünür.
- ProtFlow Web UI hem standart **CSV (,)** hem de **Excel (;)** formatında doğrudan indirme butonları sunar.
