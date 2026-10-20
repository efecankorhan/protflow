# ProtFlow: Zero-Friction Bio-Pipeline Framework
BİLİMSEL ÇIKTI ÜRETMEK İÇİN DOĞRU DEĞİLDİR ÖĞRENMEK İÇİN YAPILDI
ProtFlow; genom analizi, protein 3B yapı modellemesi, bağlanma cebi tespiti, biyofiziksel yüzey filtreleme, çoklu sekans hizalama (MSA) ve kodon optimizasyonunu tek bir çatı altında birleştiren, tipli sözleşmelerle çalışan hesaplamalı biyoloji ve protein mühendisliği boru hattı omurgasıdır.

Harici araçları (yerel binary, Docker konteyneri, uzak REST API) sıfır sürtünmeyle bağlar. C++ derleyiciniz veya GPU'nuz olmasa dahi **şeffaf düşüş (fallback)** katmanı sayesinde tüm hatları standart bir dizüstü bilgisayarda eksiksiz çalıştırır.

---

##  1. Sıfırdan Adım Adım Kurulum (Hiç Bilmeyen Biri İçin)

Her şeyi en baştan, temiz bir makinede adım adım kuralım.

### Adım 1: Önkoşullar
Bilgisayarınızda şunların kurulu olduğundan emin olun:
1. **Python 3.10 veya üzeri** ([python.org](https://www.python.org/downloads/) - Kurarken *"Add Python to PATH"* kutucuğunu işaretleyin).
2. **Git** ([git-scm.com](https://git-scm.com/)).

---

### Adım 2: Repoyu Bilgisayara İndirme
ProtFlow, ekosistemdeki popüler 7 aracı (`tools/` klasöründe) resmi Git alt modülleri olarak barındırır.

**Terminali (PowerShell, CMD veya Bash) açın ve şunu yazın:**
```bash
# Repoyu alt modülleriyle (submodules) birlikte tek seferde indirin:
git clone --recursive https://github.com/your-username/protflow.git
cd protflow
```

> **Eğer repoyu önceden düz `git clone` ile indirdiyseniz ve `tools/` klasörleri boş görünüyorsa şu komutu verin:**
> ```bash
> git submodule update --init --recursive
> ```

---

### Adım 3: Sanal Ortam (Virtual Environment) Oluşturma
Sisteminizdeki diğer Python paketlerinin karışmaması için temiz bir sanal ortam açıyoruz:

* **Windows (PowerShell):**
  ```powershell
  python -m venv .venv
  .venv\Scripts\Activate.ps1
  ```
  *(Eğer script çalıştırma yetki hatası alırsanız: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass` komutunu verip tekrar deneyin).*

* **Windows (CMD):**
  ```cmd
  python -m venv .venv
  .venv\Scripts\activate.bat
  ```

* **macOS / Linux:**
  ```bash
  python3 -m venv .venv
  source .venv/bin/activate
  ```

---

### Adım 4: Bağımlılıkları Yükleme ve Paketi Derleme
Sanal ortam aktifken şu iki komutu verin:

```bash
# 1. Temel kütüphaneleri yükleyin (pandas, pyarrow, gemmi, pyrodigal, pyyaml, requests)
pip install -r requirements.txt

# 2. ProtFlow'u geliştirici (editable) modda derleyip sisteme kurun
pip install -e .
```

Kurulumun başarılı olduğunu doğrulamak için:
```bash
python -c "import protflow; print('ProtFlow basariyla kuruldu!')"
```

---

##  2. Nasıl Çalıştırılır? (3 Farklı Kullanım Modu)

ProtFlow'u ihtiyacınıza göre 3 farklı şekilde kullanabilirsiniz:

###  Mod A: Görsel Web Stüdyosu (ComfyUI Tarzı Node Editor)
Hiç kod yazmadan, tarayıcı üzerinden blokları bağlayarak çalıştırmak için:

```bash
python run_ui.py
```
Konsolda `http://127.0.0.1:8080` adresini göreceksiniz. Tarayıcınızda (Chrome, Edge, Firefox) bu adresi açın.

**Web Arayüzünde Neler Yapabilirsiniz?**
1. **Canlı NCBI Genom İndirici:** Arama kutusuna bakteri adı (örn. *Thermus thermophilus*, *Bacillus subtilis*) veya accession no (örn. `NC_000913.3`) yazıp tek tıkla genomu doğrudan NCBI sunucularından çekin.
2. **Hazır Bakteri Profilleri:** *E. coli*, *B. subtilis*, *C. glutamicum* gibi endüstriyel konakçıları seçin.
3. **Tek Tıkla Boru Hattını Çalıştırma:** Genom $\rightarrow$ Gen Bulma $\rightarrow$ Filtreleme $\rightarrow$ AlphaFold 3B Yapı $\rightarrow$ Cep Bulma $\rightarrow$ Yüzey SASA $\rightarrow$ Homoloji $\rightarrow$ Önceliklendirme adımlarını saniyeler içinde baştan sona çalıştırın.
4. **Çıktıları Önizleme & İndirme:** Düğüm kartları üzerindeki butonlardan FASTA sekanslarını, GFF3 anotasyonlarını, 3B PDB dosyalarını ve CSV tablolarını indirin.
5. **Türkçe Excel Desteği:** Tek tıkla Türkçe Excel'de otomatik sütunlara ayrılan `;` ayrılmış CSV indirin.
6. **Düğüm Denetçisi (Inspector):** İstediğiniz düğümün YAML kodunu sağ panelde canlı düzenleyin ve `/api/manifest/test` ile test edin.

---

###  Mod B: Konsoldan Tek Komutla Referans Hat
Terminalden hızlı bir analiz koşturmak isterseniz:

```bash
python run_reference_pipeline.py
```
Bu betik `bacterial_reference_pipeline.yaml` dosyasını okur, temsili bir enzim lokusunu analiz eder ve konsola adayın **Öncelik Skoru**, **Tespit Edilen Cepleri** ve **AlphaFold Güvenilirlik Bilgilerini** tablo halinde basar.

---

###  Mod C: Optimizasyon, Kalibrasyon & Benchmark Süpürme
P2Rank mi yoksa fpocket mı daha iyi sonuç veriyor? pLDDT eşiğini 70 mi 80 mi yapmalıyız? Pareto kıyası için:

```bash
python run_optimization_and_benchmark.py
```
Swiss-Prot kıyas seti (EC 1-6 enzim sınıfları ve yapısal negatif kontroller) üzerinde tüm hat varyantlarını tarar, F1/Doğruluk/Gecikme metriklerine göre en iyi yapılandırmayı seçer ve kalibrasyon tablosunu üretir.

---

##  3. Yedi Çekirdek Biyo-Veri Tipi (`protflow.types`)

ProtFlow içinde veri, başıboş sözlükler veya kontrolsüz dosya yolları olarak değil, tip garantili sınıflar olarak akar. Her nesne bir `protein_id` taşır ve kökenini (`Provenance`: hangi araç, hangi sürüm, hangi parametreyle üretti) saklar:

| Veri Tipi | Açıklama | Standart Dosya Formatı |
| :--- | :--- | :--- |
| **`Sequence`** | Protein veya nükleotid dizilimi | `.fasta` |
| **`SequenceCollection`**| Çoklu sekans kütüphanesi / proteom | `.fasta` (Multi-FASTA) |
| **`Structure`** | 3B atomik koordinatlar (Gemmi destekli) | `.cif`, `.mmcif`, `.pdb` |
| **`ResidueScore`** | Kalıntı bazlı metrikler (pLDDT, esneklik) | `.parquet`, `.csv` |
| **`Pocket`** | Tahmin edilen aktif/bağlanma cepleri | `.parquet`, `.csv` |
| **`Alignment`** | Sekans/yapı hizalama ve benzerlik sonuçları | Tabular / Parquet |
| **`Annotation`** | Gen bölgeleri, CDS koordinatları, ekzonlar | `.gff3` |
| **`Table`** | Genel metrikler, CAI skorları, öncelik tabloları | `.parquet`, `.csv` |

---

## 4. Dahil Edilen Harici Araçlar (`tools/`) ve Akıllı Fallback

`tools/` klasöründe dünyanın en çok kullanılan açık kaynaklı biyoinformatik motorları Git submodule olarak yer alır:

| Araç | Konum | Görevi | Normal Çalışma | Sistemde Binary Yoksa Ne Olur? (Fallback) |
| :--- | :--- | :--- | :--- | :--- |
| **ColabFold** | `tools/colabfold` | De novo AlphaFold2 3B modelleme | GPU / Yerel Batch | AlphaFold DB REST API veya Geometrik Fallback |
| **P2Rank** | `tools/p2rank` | Makine öğrenmesiyle aktif cep bulma | Java / prank binary | Dahili Gemmi 3B Kavitasyon Motoru |
| **fpocket** | `tools/fpocket` | Voronoi geometrisiyle cep tespiti | C ikili dosyası | Dahili Gemmi 3B Kavitasyon Motoru |
| **Foldseek** | `tools/foldseek` | 3Di yapısal alfabe ile 3B yapı eşleme | C++ binary | Dahili Swiss-Prot Yapısal Eşleme Modeli |
| **MMseqs2** | `tools/mmseqs2` | Ultra hızlı homoloji ve sekans arama | C++ binary | Swiss-Prot Referans Arama Motoru |
| **MAFFT** | `tools/mafft` | Altın standart çoklu sekans hizalama | C binary | Resmi EMBL-EBI REST API $\rightarrow$ Saf Python BLOSUM62 Aligner |
| **FAMSA v2** | `tools/famsa` | Milyonlarca sekans için ultra hızlı MSA | C++ binary | MAFFT / EBI REST $\rightarrow$ Saf Python BLOSUM62 Aligner |

###  "Benim bilgisayarımda GPU yok, C++ derleyicisi (GCC) yok. Çalıştırabilir miyim?"
**EVET, %100 ÇALIŞTIRABİLİRSİNİZ.**  
ProtFlow'un mimarisi "zero-friction" ilkesine dayanır:
1. **Seviye 1 (Yerel İkili):** Bilgisayarınızda veya `tools/` altında derlenmiş binary varsa doğrudan onu çalıştırır.
2. **Seviye 2 (Uzak REST API):** Binary yoksa resmi araştırma enstitülerinin ücretsiz API'lerine bağlanır (örn. AlphaFold DB API, EMBL-EBI MAFFT API).
3. **Seviye 3 (Çevrimdışı Saf Python Fallback):** İnternetiniz dahi olmasa Python standard library + Gemmi + BLOSUM62 dinamik programlama motorları devreye girer. Boru hattı asla hata verip yarıda kesilmez.

Ayrıntılı C/C++ derleme yönergeleri için: [`tools/README.md`](tools/README.md)

---

## 🇹🇷 5. Excel & Bölgesel Ayar Uyumluluğu (Türkçe Windows Çözümü)

Biyoinformatik dünyasının evrensel standardı RFC-4180'dir (virgül `,` sütun ayırıcı, nokta `.` ondalık ayracı).  
Fakat **Türkçe Windows** kurulu bilgisayarlarda Microsoft Excel varsayılan olarak noktalı virgül (`;`) bekler; normal bir CSV'ye çift tıkladığınızda tüm veriyi tek bir sütuna yapıştırır ve sayıları bozar.

ProtFlow bu sorunu kökten çözer:
* Web arayüzünde hem **CSV (,)** (Python, R, AWK için) hem de **Excel (;)** butonları bulunur.
* **Excel (;)** seçildiğinde dosya otomatik olarak noktalı virgül (`;`) ayırıcı, virgül (`,`) ondalık basamak ve Türkçe karakterlerin bozulmasını önleyen **UTF-8 BOM (`\ufeff`)** başlığıyla indirilir. Excel'de çift tıkladığınızda sütunlar ve sayılar jilet gibi ayrılır.

---

##  6. Sunucu & HPC Dışa Aktarımı (Snakemake & Nextflow)

Yerel bilgisayarınızda hazırladığınız bir boru hattını üniversite veya kurum sunucusuna (HPC / SLURM) taşımak için tek satır kod yeterlidir:

```python
from protflow.pipeline import Pipeline

pipe = Pipeline.from_yaml("bacterial_reference_pipeline.yaml")

# Snakemake kural dosyası üretir:
with open("Snakefile", "w") as f:
    f.write(pipe.export_snakemake())

# Nextflow DSL2 scripti üretir:
with open("main.nf", "w") as f:
    f.write(pipe.export_nextflow())
```

---

##  7. Testleri Çalıştırma

Tüm sistemin doğruluğunu, tiplerini, güvenlik kontrollerini ve fallback mekanizmalarını test etmek için:

```bash
python -m pytest -v
```

**Sonuç:** `49 passed, 2 skipped` (Atlanan 2 test yalnızca sistemde harici Java JRE veya mmseqs binary'si kurulu olmadığında güvenle atlanan testlerdir; tüm çekirdek işlevler %100 başarıyla geçer).

---

##  8. Lisans & Telif Güvencesi

ProtFlow çekirdek motoru son derece esnek MIT lisansına sahiptir. `tools/` altındaki GPL lisanslı araçlar (MMseqs2, Foldseek, FAMSA) ProtFlow koduna doğrudan bağlanmaz (link edilmez); işletim sistemi düzeyinde bağımsız birer alt süreç (`subprocess CLI`) olarak çağrılır. Bu sayede boru hattınız GPL copyleft kısıtlamalarına maruz kalmaz, ticari ve akademik kullanımda tamamen özgür kalır.
