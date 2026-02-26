# FLAIR-HUB R — Occupation du sol multimodale depuis les données IGN

Code R pour exploiter le dataset **[FLAIR-HUB](https://huggingface.co/datasets/IGNF/FLAIR-HUB)** de l'IGN et produire des cartes d'occupation du sol et de cultures à partir des **ortho IGN** (RVB + IRC à 0.20m) en utilisant les modèles pré-entraînés FLAIR-HUB (Swin, ConvNeXTV2, U-TAE).

## Contexte

Le projet **FLAIR-HUB** de l'IGN est le plus grand dataset multi-capteurs d'occupation du sol avec des annotations très haute résolution (20 cm), couvrant **2 528 km²** de France avec **63 milliards de pixels annotés**.

### FLAIR-HUB en bref

| | FLAIR-HUB |
|---|---|
| **Couverture** | 2 528 km² (France) |
| **Annotations** | 63 milliards de pixels |
| **Résolution VHR** | 0.20 m/pixel |
| **Modalités** | 6 (aérien, SPOT, Sentinel-1/2, MNT, historique) |
| **Classes OCS** | 19 (CoSIA) |
| **Classes cultures** | 23 (LPIS/RPG) |
| **ROIs** | 2 822 |
| **Patches** | 512×512 px (102.4m × 102.4m) |
| **Format** | GeoTIFF |
| **Licence** | Etalab 2.0 |

### 6 modalités alignées

| Modalité | Résolution | Bandes | Type |
|---|---|---|---|
| **Aérien RGBI** | 0.20 m | R, G, B, NIR | UInt8 |
| **Aérien historique** | 0.40 m | Panchromatique (1950s) | UInt8 |
| **MNT (DSM/DTM)** | 0.20 m | DSM, DTM | Float32 |
| **SPOT RGBI** | 1.60 m | R, G, B, NIR | UInt16 BOA |
| **Sentinel-2 SITS** | 10 m | Séries temporelles | Multi-bandes |
| **Sentinel-1 SITS** | 10 m | Séries temporelles | Multi-bandes |

### 2 supervisions

- **CoSIA** : 19 classes d'occupation du sol (photo-interprétation experte)
- **LPIS** : 23 classes de cultures (données déclaratives RPG)

## Structure du projet

```
├── R/
│   ├── 01_download_flair_hub.R      # Téléchargement HF + chargement modalités
│   ├── 02_analyse_flair_hub.R       # Analyse, indices spectraux, visualisation
│   ├── 03_prediction_flair_hub.R    # Inférence modèles segmentation sémantique
│   └── 04_pipeline_aoi_to_landcover.R  # Pipeline complet AOI → carte OCS
├── data/
│   ├── flair_hub/                   # Dataset HF (patches FLAIR-HUB)
│   ├── ign/                         # Ortho IGN (RVB + IRC, 0.20m)
│   └── aoi.gpkg                     # Zone d'intérêt (polygone)
├── outputs/                         # Résultats et graphiques
├── DESCRIPTION                      # Métadonnées R package
├── NAMESPACE                        # Exports R
├── LICENSE                          # Licence Etalab 2.0
└── README.md
```

## Prérequis

### Packages R

```r
install.packages(c("terra", "sf", "httr2", "jsonlite", "curl", "fs"))
install.packages("reticulate")  # Interface Python
```

### Environnement Python (conda + PyTorch)

```bash
conda create -n FLAIRHUB python=3.10
conda activate FLAIRHUB
pip install torch==2.6.0 torchvision==0.21.0 torchaudio==2.6.0 \
  --extra-index-url https://download.pytorch.org/whl/cu126
pip install segmentation-models-pytorch timm \
  rasterio huggingface_hub numpy
```

### Connexion R → Python

```r
library(reticulate)
use_condaenv("FLAIRHUB", required = TRUE)
```

## Quelles données télécharger ?

L'aérien RGBI seul donne déjà **97% de la performance maximale** (64.1% vs 65.8% mIoU). Ajouter Sentinel, SPOT, etc. apporte un gain marginal.

| Config | Données nécessaires | mIoU | Recommandation |
|---|---|---|---|
| **LC-A** | Aérien RGBI seul | 64.1% | **Premier test** |
| **LC-B** | Aérien + MNT (DSM/DTM) | 65.1% | Meilleur rapport perf/complexité |
| **LC-D** | Aérien + MNT + Sentinel-2 | ~65% | Si cultures (labouré, etc.) |
| **LC-L** | Toutes modalités | 65.8% | Performance maximale |

## Utilisation

### 0. Tester rapidement avec le TOY DATASET

```r
source("R/01_download_flair_hub.R")

# Télécharger le petit jeu de données de test
toy_dir <- download_toy_dataset()

# Lister les fichiers disponibles par modalité
scan_flair_files(toy_dir)

# Voir toutes les configurations et leurs performances
show_configs()
```

### 1. Télécharger les données FLAIR-HUB

```r
source("R/01_download_flair_hub.R")

# --- Méthode rapide : télécharger par configuration ---
# Config LC-A (aérien seul, recommandé pour commencer)
download_config("LC-A", n_patches = 10)

# Config LC-B (aérien + MNT, meilleur rapport perf/complexité)
download_config("LC-B", n_patches = 10)

# Config LC-D (aérien + MNT + Sentinel-2, pour les cultures)
download_config("LC-D", n_patches = 10)

# Config LC-L (toutes modalités, performance max)
download_config("LC-L", n_patches = 5)

# --- Méthode manuelle : télécharger par modalité ---
download_flair_hub_subset(modality = "AERIAL_RGBI", n_patches = 10)
download_flair_hub_subset(modality = "DEM_ELEV", n_patches = 10)
download_flair_hub_subset(modality = "SENTINEL2_TS", n_patches = 10)

# Charger les différentes modalités
aerial <- load_aerial_rgbi("data/flair_hub/patch_aerial.tif")
spot   <- load_spot_rgbi("data/flair_hub/patch_spot.tif")
dem    <- load_dem_elev("data/flair_hub/patch_dem.tif")
s2     <- load_sentinel2_sits("data/flair_hub/patch_s2.tif")
labels <- load_label_cosia("data/flair_hub/patch_label.tif")
```

### 2. Analyser (indices spectraux, visualisation)

```r
source("R/02_analyse_flair_hub.R")

# NDVI depuis l'aérien RGBI
ndvi  <- compute_ndvi(aerial)
gndvi <- compute_gndvi(aerial)
savi  <- compute_savi(aerial)

# Masque de végétation
veg <- mask_vegetation(ndvi, threshold = 0.3)

# Visualisation multi-modalités
plot_multimodal_comparison(aerial, labels, spot = spot, dem = dem)

# Statistiques d'occupation du sol
lc_stats <- compute_landcover_stats(labels)

# Croisement OCS × hauteur (MNT)
cross <- cross_landcover_dem(labels, dem)
```

### 3. Prédire l'occupation du sol

```r
source("R/03_prediction_flair_hub.R")

# Configurer Python
setup_conda_env("FLAIRHUB")

# Télécharger un modèle pré-entraîné
# Option 1 : FLAIR-INC (simple, ResNet34-UNet, recommandé pour commencer)
model_path <- download_pretrained_model("FLAIR-INC_rgbi_15cl_resnet34-unet")

# Option 2 : FLAIR-HUB (multimodal, ConvNeXTV2-UPerNet)
# model_path <- download_pretrained_model("FLAIR-HUB_LC-A_IR_convnextv2tiny-upernet")

# Prédire depuis une image aérienne RGBI
landcover <- predict_landcover_from_aerial("data/ign/ortho_rgbi.tif", model_path)
```

### 4. Pipeline complet AOI → Carte d'occupation du sol

```r
source("R/04_pipeline_aoi_to_landcover.R")

# Config LC-A : RGBI seul (64.1% mIoU, rapide)
result <- pipeline_aoi_to_landcover("data/aoi.gpkg")

# Config LC-B : RGBI + MNT à 1m (65.1% mIoU, +1pt, recommandé)
result <- pipeline_aoi_to_landcover("data/aoi.gpkg",
  use_dem = TRUE, dem_res_m = 1)

# Le résultat contient :
# result$ortho_rvb   - Ortho RVB (0.20m)
# result$ortho_irc   - Ortho IRC (0.20m)
# result$ortho_rgbi  - Ortho RGBI 4 bandes (0.20m)
# result$ndvi        - NDVI
# result$dem         - MNT DSM+DTM (si use_dem=TRUE)
# result$landcover   - Carte d'occupation du sol
```

**Fichiers produits dans `outputs/` :**

| Fichier | Description |
|---|---|
| `ortho_rvb.tif` | Ortho RVB IGN (0.20m) |
| `ortho_irc.tif` | Ortho IRC IGN (0.20m) |
| `ortho_rgbi.tif` | Ortho RGBI 4 bandes (0.20m) |
| `dem_dsm_dtm.tif` | MNT DSM+DTM rééchantillonné à 0.2m (si `use_dem=TRUE`) |
| `ndvi.tif` | NDVI calculé depuis l'IRC |
| `landcover_predicted.tif` | Carte d'occupation du sol prédite |
| `resultats_aoi_flair_hub.pdf` | Visualisation récapitulative (4 panneaux) |

## Classes d'occupation du sol (CoSIA)

| ID | Classe | Couleur |
|---|---|---|
| 1 | Bâtiment | ![#db0e9a](https://via.placeholder.com/15/db0e9a/db0e9a.png) |
| 2 | Serre | ![#938e7b](https://via.placeholder.com/15/938e7b/938e7b.png) |
| 3 | Piscine | ![#f80c00](https://via.placeholder.com/15/f80c00/f80c00.png) |
| 4 | Surface imperméable | ![#a97101](https://via.placeholder.com/15/a97101/a97101.png) |
| 5 | Surface perméable | ![#1553ae](https://via.placeholder.com/15/1553ae/1553ae.png) |
| 6 | Sol nu | ![#194a26](https://via.placeholder.com/15/194a26/194a26.png) |
| 7 | Eau | ![#46e483](https://via.placeholder.com/15/46e483/46e483.png) |
| 8 | Neige | ![#f3a60d](https://via.placeholder.com/15/f3a60d/f3a60d.png) |
| 9 | Végétation herbacée | ![#660082](https://via.placeholder.com/15/660082/660082.png) |
| 10 | Terre agricole | ![#55ff00](https://via.placeholder.com/15/55ff00/55ff00.png) |
| 11 | Terre labourée | ![#fff30d](https://via.placeholder.com/15/fff30d/fff30d.png) |
| 12 | Vigne | ![#e4df7c](https://via.placeholder.com/15/e4df7c/e4df7c.png) |
| 13 | Verger | ![#3de6eb](https://via.placeholder.com/15/3de6eb/3de6eb.png) |
| 14 | Feuillu | ![#ffffff](https://via.placeholder.com/15/ffffff/ffffff.png) |
| 15 | Conifère | ![#8ab3a0](https://via.placeholder.com/15/8ab3a0/8ab3a0.png) |
| 16 | Lande | ![#6b714f](https://via.placeholder.com/15/6b714f/6b714f.png) |
| 17 | Ligneux mélangé | ![#c5dc42](https://via.placeholder.com/15/c5dc42/c5dc42.png) |
| 18 | Fleur / Garrigue | ![#9999ff](https://via.placeholder.com/15/9999ff/9999ff.png) |
| 19 | Non classé | ![#000000](https://via.placeholder.com/15/000000/000000.png) |

## Modèles pré-entraînés disponibles

### FLAIR-INC (simples, recommandés pour tester)

| Modèle HF | Entrée | Classes | Encodeur |
|---|---|---|---|
| `IGNF/FLAIR-INC_rgbi_15cl_resnet34-unet` | RGBI (4 bandes) | 15 | ResNet34-UNet |
| `IGNF/FLAIR-INC_rgbie_15cl_resnet34-unet` | RGBI+E (5 bandes) | 15 | ResNet34-UNet |
| `IGNF/FLAIR-INC_rgb_15cl_resnet34-unet` | RGB (3 bandes) | 15 | ResNet34-UNet |
| `IGNF/FLAIR-INC_rgb_12cl_resnet34-unet` | RGB (3 bandes) | 12 | ResNet34-UNet |

### FLAIR-HUB (multimodaux, plus performants)

| Modèle HF | Encodeur | Décodeur | Entrée |
|---|---|---|---|
| `IGNF/FLAIR-HUB_LC-G_utae` | U-TAE | UperFuse | Multimodal (S2 SITS) |
| `IGNF/FLAIR-HUB_LC-A_IR_convnextv2tiny-upernet` | ConvNeXTV2-T | UPerNet | Aérien RGBI |

Collection complète : [huggingface.co/collections/IGNF/flair-models](https://huggingface.co/collections/IGNF/flair-models-684035e78bd5bff99199ff87)

## Références

- **FLAIR-HUB** : Garioud, A., Giordano, S., David, N., & Gonthier, N. (2025). [arXiv:2506.07080](https://arxiv.org/abs/2506.07080)
- **Dataset HF** : [IGNF/FLAIR-HUB](https://huggingface.co/datasets/IGNF/FLAIR-HUB)
- **Code Python** : [IGNF/FLAIR-HUB](https://github.com/IGNF/FLAIR-HUB)
- **BD ORTHO® IGN** : [geoservices.ign.fr/bdortho](https://geoservices.ign.fr/bdortho)
- **Géoplateforme** : [data.geopf.fr](https://data.geopf.fr)
- **Page projet FLAIR** : [ignf.github.io/FLAIR](https://ignf.github.io/FLAIR/)
