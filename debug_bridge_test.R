# ─────────────────────────────────────────────────────────────────────────────
# BRIDGE: Hızlı Hata Ayıklama (Debug) ve Test Konfigürasyonu
# ─────────────────────────────────────────────────────────────────────────────
# Bu script, bridge.R dosyasında yapılan tüm devasa mimari değişikliklerin
# (PCA hızlandırması, bellek optimizasyonu, p/n filtrelemesi, Jaccard vb.)
# başarıyla çalışıp çalışmadığını, "saatlerce donma" yaşatmadan saniyeler
# içinde test edebilmeniz için özel olarak hazırlanmıştır.
#
# Lütfen bu dosyayı RStudio'da açıp tümünü (Ctrl+A -> Ctrl+Enter) çalıştırın.

source("R/bridge.R")
source("R/utils.R")

# 1. Küçük / Sentetik Bir Veri Seti Oluşturulması
set.seed(42)
message("\n[1] Sentetik veri seti (100 Gen x 40 Hasta) oluşturuluyor...")
# 100 Gen, 40 Hasta (20'si A sınıfı, 20'si B sınıfı)
X <- matrix(rnorm(100 * 40), nrow = 100, 
            dimnames = list(paste0("G", 1:100), paste0("S", 1:40)))
y <- factor(rep(c("A", "B"), each = 20))

# 2. Sentetik Hastalık Sinyali (Differential Module) Eklenmesi
# G1'den G10'a kadar olan genler, B sınıfı hastalarda aktifleşsin.
X[1:10, 21:40] <- X[1:10, 21:40] + 1.5

# 3. Sentetik Biyolojik Ön Bilgi (Prior)
# G1'den G12'ye kadar olan genlerin literatürde etkileşimde olduğu bilinsin.
B_prior <- list(known_disease_pathway = paste0("G", 1:12))

message("\n[2] BRIDGE motoru (DEBUG MODE) başlatılıyor...\n")
# ─────────────────────────────────────────────────────────────────────────────
# DİKKAT: DEBUG MODU PARAMETRELERİ
# ─────────────────────────────────────────────────────────────────────────────
# - prefilter = "variance": Label-blind hızlandırma için en varyanslıları alır.
# - max_genes = 50: Veriyi 50 gene kısıtlar (Hızlı test için).
# - resolution = 1: Tek bir Leiden çözünürlüğü tarar (Kapsamlı grid yerine).
# - n_boot = 5: Stabilite testini sadece 5 kez yapar (Normalde 20 veya 100 olmalı).
# - tau = 0.5: Toleranslı stabilite skoru (Zor testlerde 0.6 veya 0.7 olmalı).
# ─────────────────────────────────────────────────────────────────────────────

res <- bridge_grouping(
  X = X, 
  y = y, 
  B = B_prior, 
  prefilter = "variance", 
  max_genes = 50,          # Debug: Düşük gen sayısı
  resolution = 1,          # Debug: Tek resolution
  n_boot = 5,              # Debug: Düşük bootstrap sayısı
  tau = 0.5,               # Debug: Esnek stabilite eşiği
  verbose = TRUE,
  details = TRUE           # Teşhis (diagnostics) tablolarını görmek için
)

message("\n[3] BRIDGE Analizi Tamamlandı! Sonuçlar:\n")

# A) Seçilen Modüller
print(str(res))

# B) Teşhis ve Çalışma Süresi Raporları (Diagnostics)
message("\n=== Modül Teşhis (Diagnostics) Raporu ===")
print(attr(res, "diagnostics"))

message("\n=== Test Başarıyla Sona Erdi ===")
