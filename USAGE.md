# ProtFlow Kullanım Kılavuzu (User Manual)

ProtFlow; hesaplamalı biyoloji, protein mühendisliği ve enzim tarama iş akışlarını modüler, tipli ve görsel bir düğüm (node) editörüyle tasarlamanızı ve yürütmenizi sağlayan platformdur.

---

## 1. Hızlı Başlangıç

### Gereksinimler
- Python 3.10+
- Temel bağımlılıklar: `gemmi`, `pandas`, `pyarrow`, `pyyaml`, `requests`, `numpy`

### Web Arayüzünü (Studio) Başlatma
Terminalinizde repo ana dizinindeyken:
```bash
python run_ui.py
```
Sunucu başladığında tarayıcınızda şu adresi açın:
```
http://127.0.0.1:8080
```

### Komut Satırından (CLI) Çalıştırma
Görsel arayüz kullanmadan doğrudan arka planda çalıştırmak için:
- **Referans Hattı Koşma:**
  ```bash
  python run_reference_pipeline.py
  ```
- **Kıyaslama ve Optimizasyon (Benchmark Sweep):**
  ```bash
  python run_optimization_and_benchmark.py
  ```
- **Tüm Testleri Koşma:**
  ```bash
  pytest
  ```

---

## 2. Görsel Editör (Canvas) Kullanımı

ProtFlow Studio, soldan sağa akan veri hatlarını (DAG) düğümler ve kablolarla yönetir.

### Düğüm Taşıma ve Gezinme
- **Kaydırma (Pan):** Canvas boşluğuna sol tıklayıp sürükleyin.
- **Yakınlaştırma (Zoom):** Fare tekerleği veya sol alttaki `+`, `-`, `Sığdır` düğmelerini kullanın.
- **Otomatik Düzenle:** Üst barda **"Otomatik Düzenle"** düğmesine tıklayarak tüm hattı katmanlı (soldan sağa) hizalayın.

### Portlar ve Tip Güvenliği
Her bağlantı noktası (port) bir veri tipini temsil eder ve özel rengiyle gösterilir:
| Veri Tipi | Renk | Anlamı |
|---|---|---|
| **Sequence** | Turkuaz (`#00cec9`) | FASTA formatında protein veya nükleotid sekansı |
| **Structure** | Turuncu (`#f39c12`) | 3B atom koordinatları (mmCIF / PDB) |
| **ResidueScore** | Sarı (`#f1c40f`) | Kalıntı düzeyinde skorlar (örn. pLDDT) |
| **Pocket** | Mor (`#9b59b6`) | Bağlanma cebi koordinatları ve hacimleri |
| **Alignment** | Mavi (`#0984e3`) | Yapısal veya sekans eşleşme kayıtları |
| **Annotation** | Pembe (`#e84393`) | Genomik ve özellik anotasyonları (GFF3) |
| **Table** | Zümrüt Yeşili (`#2ecc71`) | Genel metrikler, filtreleme ve skor tabloları |

> **Kural:** Yalnızca uyumlu tipler birbirine bağlanabilir. Uyumsuz bir bağlantı çekildiğinde arayüz işlemi iptal eder ve gerekçesini bildirimle açıklar.

### Bağlantı (Kablo) Kurma ve Silme
- **Bağlama:** Bir düğümün çıktı portuna (`out`) tıklayıp farenizi bırakmadan hedef düğümün girdi portuna (`in`) sürükleyip bırakın.
- **Silme (3 Yöntem):**
  1. **Çarpı Butonu (✕):** Fare imlecini kablonun ortasına getirdiğinizde beliren kırmızı `✕` butonuna doğrudan tıklayın.
  2. **Klavye:** Kabloya tıklayarak seçin (kırmızı renge döner), ardından klavyenizdeki `Delete` veya `Backspace` tuşuna basın.
  3. **Çift Tıklama / Sağ Tıklama:** Kablonun üzerine çift tıklayın veya sağ tıklayın.

---

## 3. Düğüm Denetim Paneli (Inspector)

Bir düğüme tıkladığınızda ekranın sağ tarafında denetim paneli açılır:
1. **Düğüm Başlığı & Çalışma Modu:** Düğümün adını değiştirebilir; `LOCAL`, `CONTAINER` veya `HTTP` modlarını seçebilirsiniz.
2. **Çalıştırma Komutu / URL:** Aracın yerel binary yolunu veya uzak REST API adresini düzenleyebilirsiniz.
3. **Parametreler:** Eşik değerlerini (örn. `min_plddt: 70.0`) anlık güncelleyebilirsiniz.
4. **Hızlı Araç Takası (1-Tık):** Aynı girdi ve çıktı tiplerine sahip alternatif araçlar (örneğin P2Rank ↔ fpocket) tek tıkla takas edilebilir; mevcut bağlantılar korunur.
5. **Örnek Girdiyle Test Et:** Düğümü hatta sokmadan önce bağımsız olarak test eder, tip uyumunu doğrular ve alt konsolda logları gösterir.

---

## 4. Özel Düğüm / Araç Ekleme

Yeni bir hesaplama aracı eklemek için Python veya UI koduna dokunmanıza gerek yoktur; yalnızca `manifests/` klasörüne bir YAML dosyası eklemeniz yeterlidir.

### Örnek Manifest Şablonu (`manifests/my_tool.yaml`)
```yaml
name: my_custom_tool
version: "1.0.0"
description: "Özel protein analiz aracı"
is_remote: false # Uzak API ise true (Gizlilik koruması devreye girer)
inputs:
  structure: Structure
outputs:
  analysis_table: Table
parameters:
  cutoff: 5.0
execution:
  type: local # local | container | http
  command: "my_binary --input {inputs.structure.file_path} --cutoff {parameters.cutoff} -o {output_dir}"
adapter: protflow.adapters.pipeline_tools:parse_my_tool
```

### Arayüzden Ekleme
1. Sol paneldeki **"+ Özel Düğüm Ekle"** butonuna tıklayın.
2. YAML tanımını yapıştırın veya sihirbaz alanlarını doldurun.
3. **"Doğrula & Test Et"** butonuna basarak şema ve tip kontrolü yapın.
4. **"Kaydet"** butonuna bastığınızda araç sol kütüphaneye ve canvas'a eklenir.

---

## 5. Referans Hat ve Sonuçları İnceleme

Üst bardaki **"Hattı Çalıştır"** butonuna tıkladığınızda:
1. Uzak servis içeren düğümler (`REMOTE`) varsa gizlilik onayı modalı açılır.
2. Onay verildiğinde düğümler sırayla çalışır ve durum rozetleri güncellenir (`running` → `PASSED`).
3. Filtre düğümlerinin altında elenen protein sayıları görünür (örn. `Girdi: 1 → Kalan: 1`).
4. Hat tamamlandığında **"Enzim Aday Öncelik Listesi"** açılır:
   - Nihai skor (`final_priority_score`)
   - Yapı kalitesi (`structure_quality`)
   - Cep skoru ve hacmi (`pocket_score`, `pocket_volume`)
   - Çözünürlük ve yenilik skorları
   - Aday durumu (`CANDIDATE` / `DEPURGED`)

---

## 6. Kıyaslama (Benchmark) & Islak Lab Alanin Taraması

- **Kıyaslama & Grafik:** Üst menüden **"Kıyaslama & Grafik"** modalını açarak farklı hat varyantlarının (örn. P2Rank vs fpocket) Swiss-Prot etiketli enzim veri seti üzerindeki F1 skoru, kesinlik ve gecikme sürelerini karşılaştırabilirsiniz.
- **Islak Lab / Alanin Taraması:** **"Islak Lab / Alanin"** butonuna tıklayarak ilk adayların aktif ceplerindeki kalıntılar için üretilen mutasyonları (`A:35 GLU->ALA`) inceleyebilir ve laboratuvar ölçümlerini kalibrasyon tablosuna geri besleyebilirsiniz.

---

## 7. Sorun Giderme ve İpuçları

1. **Bağlantı Kurulamıyor:** Port renklerini kontrol edin. Farklı renkli portlar bağlanamaz.
2. **Kablo Silinmiyor:** Kabloya tıklayıp kırmızı olmasını sağlayın ve klavyeden `Delete` tuşuna basın veya ortasındaki `✕` butonuna tıklayın.
3. **Uzak Sunucu / API Hatası:** Ortam değişkenlerinizi (`$API_KEY` vb.) tanımladığınızdan ve internet bağlantınızın açık olduğundan emin olun. Offline çalışmada yerel simülasyon korumaları devrededir.
