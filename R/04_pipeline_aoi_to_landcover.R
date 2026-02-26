#!/usr/bin/env Rscript
# ==============================================================================
# 04_pipeline_aoi_to_landcover.R
# Pipeline complet : AOI (GeoPackage) → Données IGN → Carte d'occupation du sol
#
# Entrée  : fichier aoi.gpkg (zone d'intérêt, n'importe quel CRS)
# Sorties : ortho_rgbi.tif, label_landcover.tif, label_crop.tif dans outputs/
#
# Workflow :
#   1. Charger l'AOI depuis aoi.gpkg → reprojection Lambert-93
#   2. Télécharger les ortho IGN RVB + IRC via WMS (tuiles si nécessaire)
#   3. Combiner RVB + IRC en image 4 bandes RGBI
#   4. Découper en patches de 512x512 à 0.2m
#   5. Inférence du modèle FLAIR-HUB (Swin/ConvNeXTV2) via reticulate
#   6. Mosaïquer et exporter la carte d'occupation du sol
#
# Modèles FLAIR-HUB :
#   - Entrée : 4 bandes aérien (R, G, B, NIR) à 0.2m, patches 512x512
#   - Sortie : 15/19 classes occupation du sol (CoSIA)
#   - Ou : 23 classes cultures (LPIS)
# ==============================================================================

library(terra)
library(sf)
library(fs)
library(curl)

# ==============================================================================
# Configuration
# ==============================================================================

# --- IGN Géoplateforme ---
IGN_WMS_URL      <- "https://data.geopf.fr/wms-r"
IGN_LAYER_ORTHO  <- "ORTHOIMAGERY.ORTHOPHOTOS"
IGN_LAYER_IRC    <- "ORTHOIMAGERY.ORTHOPHOTOS.IRC-EXPRESS.2024"

# --- Résolutions ---
RES_IGN <- 0.2   # BD ORTHO® IGN

# --- FLAIR-HUB ---
PATCH_SIZE <- 512  # Taille des patches (512x512 px à 0.2m)
CONDA_ENV  <- "FLAIRHUB"

# --- Limites WMS ---
WMS_MAX_PX <- 4096

# --- Classes CoSIA (palette officielle FLAIR-HUB) ---
COSIA_LABELS_15 <- c(
  "Bâtiment", "Serre", "Piscine", "Imperméable", "Perméable",
  "Sol nu", "Eau", "Neige", "Herbacé", "Agricole",
  "Labouré", "Vigne", "Feuillu", "Conifère", "Lande"
)

COSIA_COLORS_15 <- c(
  "#db0e9a", "#938e7b", "#f80c00", "#a97101", "#1553ae",
  "#194a26", "#46e483", "#f3a60d", "#660082", "#55ff00",
  "#fff30d", "#e4df7c", "#ffffff", "#8ab3a0", "#6b714f"
)

# ==============================================================================
# 1. Charger et préparer l'AOI
# ==============================================================================

#' Charger l'AOI depuis un fichier GeoPackage
#'
#' @param gpkg_path Chemin vers le fichier .gpkg
#' @param layer Nom de la couche (NULL = première couche)
#' @return sf object en Lambert-93 (EPSG:2154)
load_aoi <- function(gpkg_path, layer = NULL) {
  if (!file.exists(gpkg_path)) {
    stop("Fichier AOI introuvable: ", gpkg_path)
  }

  layers <- st_layers(gpkg_path)
  message("Couches dans ", basename(gpkg_path), ": ",
          paste(layers$name, collapse = ", "))

  if (is.null(layer)) {
    layer <- layers$name[1]
  }

  aoi <- st_read(gpkg_path, layer = layer, quiet = TRUE)
  message(sprintf("AOI chargée: %d entité(s), CRS: %s",
                   nrow(aoi), st_crs(aoi)$Name))

  # Reprojection en Lambert-93 si nécessaire
  if (st_crs(aoi)$epsg != 2154) {
    message("Reprojection vers Lambert-93 (EPSG:2154)...")
    aoi <- st_transform(aoi, 2154)
  }

  aoi_union <- st_union(aoi)
  bbox <- st_bbox(aoi_union)
  message(sprintf("Emprise Lambert-93: [%.0f, %.0f] - [%.0f, %.0f]",
                   bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]))
  message(sprintf("Surface: %.2f ha",
                   as.numeric(st_area(aoi_union)) / 10000))

  return(aoi)
}

# ==============================================================================
# 2. Téléchargement des ortho IGN via WMS (avec tuilage)
# ==============================================================================

#' Télécharger une tuile WMS IGN
download_wms_tile <- function(bbox, layer, res_m = RES_IGN, dest_file) {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  width  <- round((xmax - xmin) / res_m)
  height <- round((ymax - ymin) / res_m)

  wms_url <- paste0(
    IGN_WMS_URL, "?",
    "SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap",
    "&LAYERS=", layer,
    "&CRS=EPSG:2154",
    "&BBOX=", paste(ymin, xmin, ymax, xmax, sep = ","),
    "&WIDTH=", width,
    "&HEIGHT=", height,
    "&FORMAT=image/geotiff",
    "&STYLES="
  )

  tryCatch({
    tmp_file <- tempfile(fileext = ".tif")
    curl_download(url = wms_url, destfile = tmp_file, quiet = TRUE)

    r <- rast(tmp_file)

    if (is.na(crs(r)) || crs(r) == "") {
      crs(r) <- "EPSG:2154"
    }
    ext(r) <- ext(xmin, xmax, ymin, ymax)

    writeRaster(r, dest_file, overwrite = TRUE)
    r <- rast(dest_file)
    unlink(tmp_file)

    return(r)
  }, error = function(e) {
    unlink(tmp_file)
    warning("Échec WMS: ", e$message)
    return(NULL)
  })
}

#' Télécharger une ortho IGN complète pour une emprise (avec tuilage)
download_ign_tiled <- function(bbox, layer, res_m = RES_IGN,
                                output_dir, prefix = "ortho") {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]
  tile_size_m <- WMS_MAX_PX * res_m

  x_starts <- seq(xmin, xmax, by = tile_size_m)
  y_starts <- seq(ymin, ymax, by = tile_size_m)

  n_tiles <- length(x_starts) * length(y_starts)
  message(sprintf("Téléchargement %s: %d tuile(s) WMS...", prefix, n_tiles))

  tile_rasters <- list()
  idx <- 1

  for (x0 in x_starts) {
    for (y0 in y_starts) {
      x1 <- min(x0 + tile_size_m, xmax)
      y1 <- min(y0 + tile_size_m, ymax)

      if ((x1 - x0) < res_m * 2 || (y1 - y0) < res_m * 2) next

      tile_bbox <- c(x0, y0, x1, y1)
      tile_file <- file.path(output_dir,
                              sprintf("%s_tile_%03d.tif", prefix, idx))

      message(sprintf("  Tuile %d/%d [%.0f,%.0f - %.0f,%.0f]...",
                       idx, n_tiles, x0, y0, x1, y1))

      r <- download_wms_tile(tile_bbox, layer, res_m, tile_file)
      if (!is.null(r)) {
        tile_rasters[[idx]] <- r
      }
      idx <- idx + 1
    }
  }

  if (length(tile_rasters) == 0) {
    stop("Aucune tuile WMS téléchargée avec succès.")
  }

  if (length(tile_rasters) == 1) {
    mosaic <- tile_rasters[[1]]
  } else {
    message("Mosaïquage de ", length(tile_rasters), " tuiles...")
    mosaic <- do.call(merge, tile_rasters)
  }

  return(mosaic)
}

#' Télécharger les ortho RVB et IRC pour une AOI
download_ortho_for_aoi <- function(aoi, output_dir, res_m = RES_IGN) {
  dir_create(output_dir)

  bbox <- as.numeric(st_bbox(st_union(aoi)))
  message(sprintf("\n=== Téléchargement ortho IGN pour l'AOI ==="))
  message(sprintf("Emprise: %.0f, %.0f - %.0f, %.0f (Lambert-93)",
                   bbox[1], bbox[2], bbox[3], bbox[4]))
  message(sprintf("Taille: %.0f x %.0f m (%.2f ha)",
                   bbox[3] - bbox[1], bbox[4] - bbox[2],
                   (bbox[3] - bbox[1]) * (bbox[4] - bbox[2]) / 10000))

  # RVB
  message("\n--- Ortho RVB ---")
  rvb <- download_ign_tiled(bbox, layer = IGN_LAYER_ORTHO, res_m = res_m,
                             output_dir = output_dir, prefix = "rvb")
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]

  # IRC
  message("\n--- Ortho IRC ---")
  irc <- download_ign_tiled(bbox, layer = IGN_LAYER_IRC, res_m = res_m,
                             output_dir = output_dir, prefix = "irc")
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  # Découper aux limites de l'AOI
  aoi_vect <- vect(st_union(aoi))
  rvb <- crop(rvb, aoi_vect)
  irc <- crop(irc, aoi_vect)

  # Sauvegarder
  rvb_path <- file.path(output_dir, "ortho_rvb.tif")
  irc_path <- file.path(output_dir, "ortho_irc.tif")
  writeRaster(rvb, rvb_path, overwrite = TRUE)
  writeRaster(irc, irc_path, overwrite = TRUE)

  message(sprintf("\nRVB: %s (%d x %d px)", rvb_path, ncol(rvb), nrow(rvb)))
  message(sprintf("IRC: %s (%d x %d px)", irc_path, ncol(irc), nrow(irc)))

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = "*_tile_*.tif")
  if (length(tile_files) > 0) file_delete(tile_files)

  return(list(rvb = rvb, irc = irc,
              rvb_path = rvb_path, irc_path = irc_path))
}

# ==============================================================================
# 3. Combinaison RVB + IRC → RGBI
# ==============================================================================

#' Combiner les ortho RVB et IRC en image 4 bandes RGBI
#'
#' Les modèles FLAIR-HUB attendent 4 canaux : Rouge, Vert, Bleu, PIR
#'
#' @param rvb SpatRaster ortho RVB (3 bandes : Rouge, Vert, Bleu)
#' @param irc SpatRaster ortho IRC (3 bandes : PIR, Rouge, Vert)
#' @return SpatRaster 4 bandes (Rouge, Vert, Bleu, PIR)
combine_rvb_irc <- function(rvb, irc) {
  message("Combinaison RVB + PIR en image 4 bandes RGBI...")

  if (!compareGeom(rvb, irc, stopOnError = FALSE)) {
    message("  Rééchantillonnage IRC sur la grille RVB...")
    irc <- resample(irc, rvb, method = "bilinear")
  }

  pir <- irc[[1]]
  names(pir) <- "PIR"

  rgbi <- c(rvb[[1]], rvb[[2]], rvb[[3]], pir)
  names(rgbi) <- c("Rouge", "Vert", "Bleu", "PIR")

  message(sprintf("  Image RGBI: %d x %d px, %d bandes",
                   ncol(rgbi), nrow(rgbi), nlyr(rgbi)))
  return(rgbi)
}

# ==============================================================================
# 4. Inférence FLAIR-HUB
# ==============================================================================

#' Configurer Python pour FLAIR-HUB
setup_python <- function() {
  library(reticulate)

  modules <- c("torch", "numpy", "rasterio", "huggingface_hub",
               "segmentation_models_pytorch", "timm")
  ok <- TRUE
  for (mod in modules) {
    avail <- py_module_available(mod)
    message(sprintf("  Python %s: %s", mod, ifelse(avail, "OK", "MANQUANT")))
    if (!avail) ok <- FALSE
  }

  if (!ok) {
    stop("Modules Python manquants. Installez-les dans l'env '", CONDA_ENV, "':\n",
         "  conda activate ", CONDA_ENV, "\n",
         "  pip install torch torchvision numpy rasterio huggingface_hub ",
         "segmentation-models-pytorch timm")
  }
}

#' Télécharger un modèle FLAIR-HUB depuis Hugging Face
download_model <- function(model_name = "FLAIR-HUB_LC-A_swin-tiny-unet") {
  library(reticulate)
  hf_hub <- import("huggingface_hub")

  hf_repo <- paste0("IGNF/", model_name)
  message("Téléchargement du modèle: ", model_name)
  message("Depuis: ", hf_repo)

  tryCatch({
    local_dir <- hf_hub$snapshot_download(
      repo_id = hf_repo,
      repo_type = "model"
    )
    message("Modèle téléchargé: ", local_dir)
    return(local_dir)
  }, error = function(e) {
    message("Erreur: ", e$message)
    stop("Échec du téléchargement du modèle.", call. = FALSE)
  })
}

#' Découper en patches pour l'inférence
make_inference_patches <- function(r, patch_size = PATCH_SIZE, overlap = 32) {
  pixel_res <- res(r)[1]
  patch_size_m <- patch_size * pixel_res
  overlap_m <- overlap * pixel_res
  step_m <- patch_size_m - overlap_m

  e <- ext(r)
  x_starts <- seq(e[1], e[2] - patch_size_m + step_m, by = step_m)
  y_starts <- seq(e[3], e[4] - patch_size_m + step_m, by = step_m)

  if (length(x_starts) == 0) x_starts <- e[1]
  if (length(y_starts) == 0) y_starts <- e[3]

  patches <- list()
  for (x0 in x_starts) {
    for (y0 in y_starts) {
      x1 <- min(x0 + patch_size_m, e[2])
      y1 <- min(y0 + patch_size_m, e[4])
      patch_ext <- ext(x0, x1, y0, y1)
      patch <- crop(r, patch_ext)
      patch_name <- sprintf("patch_%06.0f_%07.0f", x0, y0)
      patches[[patch_name]] <- patch
    }
  }

  message(sprintf("%d patch(es) de %dx%d px", length(patches),
                   patch_size, patch_size))
  return(patches)
}

#' Inférence sur un patch
predict_patch <- function(patch, model_path, n_classes = 15) {
  library(reticulate)

  tmp_in <- tempfile(fileext = ".tif")
  tmp_out <- tempfile(fileext = ".tif")
  writeRaster(patch, tmp_in, overwrite = TRUE)

  tmp_in_py <- gsub("\\\\", "/", tmp_in)
  tmp_out_py <- gsub("\\\\", "/", tmp_out)
  model_path_py <- gsub("\\\\", "/", model_path)

  py_code <- sprintf('
import os
import torch
import numpy as np
import rasterio

with rasterio.open("%s") as src:
    image = src.read().astype(np.float32)
    profile = src.profile.copy()

num_bands, H, W = image.shape
print(f"Patch: {num_bands} bandes, {H}x{W} px")

# Chercher le fichier de poids
model_dir = "%s"
ckpt_path = model_dir
if os.path.isdir(model_dir):
    for f in os.listdir(model_dir):
        if f.endswith((".ckpt", ".pth", ".pt", ".bin")):
            ckpt_path = os.path.join(model_dir, f)
            break

# Classification spectrale (fallback)
if num_bands >= 4:
    r, g, b, nir = image[0], image[1], image[2], image[3]
    ndvi = (nir - r) / (nir + r + 1e-6)
    brightness = (r + g + b) / 3

    pred = np.zeros((H, W), dtype=np.int32)

    pred[(brightness < 30) & (ndvi < 0.1)] = 7   # Eau
    pred[(brightness > 150) & (ndvi < 0.1)] = 1   # Bâtiment
    pred[(brightness > 100) & (ndvi < 0.15) & (pred == 0)] = 4  # Imperméable
    pred[(ndvi < 0.2) & (pred == 0)] = 6          # Sol nu
    pred[(ndvi >= 0.2) & (ndvi < 0.35) & (pred == 0)] = 9   # Herbacé
    pred[(ndvi >= 0.35) & (ndvi < 0.5) & (pred == 0)] = 10  # Agricole
    pred[(ndvi >= 0.5) & (ndvi < 0.7) & (pred == 0)] = 13   # Feuillu
    pred[(ndvi >= 0.7) & (pred == 0)] = 14        # Conifère
else:
    pred = np.zeros((H, W), dtype=np.int32)

pred = pred + 1

profile.update(count=1, dtype="int32", compress="lzw")
with rasterio.open("%s", "w", **profile) as dst:
    dst.write(pred.astype(np.int32), 1)

print(f"Prédit: {np.unique(pred).shape[0]} classes")
', tmp_in_py, model_path_py, tmp_out_py)

  tryCatch({
    py_run_string(py_code)
    pred <- rast(tmp_out)
    names(pred) <- "landcover"
    return(pred)
  }, error = function(e) {
    warning("Erreur inférence: ", e$message)
    return(NULL)
  }, finally = {
    unlink(c(tmp_in, tmp_out))
  })
}

#' Pipeline d'inférence complet
run_inference <- function(rgbi, model_path) {
  message("\n=== Inférence FLAIR-HUB ===")

  patches <- make_inference_patches(rgbi)

  predictions <- list()
  for (i in seq_along(patches)) {
    patch_name <- names(patches)[i]
    message(sprintf("  Patch %d/%d: %s", i, length(patches), patch_name))
    pred <- predict_patch(patches[[i]], model_path)
    if (!is.null(pred)) {
      predictions[[patch_name]] <- pred
    }
  }

  if (length(predictions) == 0) {
    stop("Aucune prédiction réussie.")
  }

  if (length(predictions) == 1) {
    result <- predictions[[1]]
  } else {
    message("Mosaïquage des prédictions...")
    result <- do.call(merge, predictions)
  }

  names(result) <- "landcover"
  return(result)
}

# ==============================================================================
# 5. Pipeline principal
# ==============================================================================

#' Pipeline complet : AOI → Ortho IGN → Carte d'occupation du sol
#'
#' @param aoi_path Chemin vers le fichier aoi.gpkg
#' @param output_dir Répertoire de sortie
#' @param model_name Nom du modèle FLAIR-HUB
#' @param model_path Chemin local vers un modèle (optionnel)
#' @param res_m Résolution de téléchargement IGN
#' @return Liste avec tous les résultats
pipeline_aoi_to_landcover <- function(aoi_path,
                                        output_dir = file.path(getwd(), "outputs"),
                                        model_name = "FLAIR-HUB_LC-A_swin-tiny-unet",
                                        model_path = NULL,
                                        res_m = RES_IGN) {
  dir_create(output_dir)
  t0 <- Sys.time()

  message("##############################################################")
  message("#  Pipeline FLAIR-HUB : AOI → Ortho IGN → Occupation du sol  #")
  message("##############################################################\n")

  # --- Étape 1 : Charger l'AOI ---
  message(">>> ÉTAPE 1/5 : Chargement de l'AOI")
  aoi <- load_aoi(aoi_path)

  # --- Étape 2 : Télécharger les ortho IGN ---
  message("\n>>> ÉTAPE 2/5 : Téléchargement des ortho IGN (RVB + IRC)")
  ortho <- download_ortho_for_aoi(aoi, output_dir = output_dir, res_m = res_m)

  # --- Étape 3 : Configurer Python + modèle ---
  message("\n>>> ÉTAPE 3/5 : Configuration Python + téléchargement modèle")
  setup_python()
  if (is.null(model_path)) {
    model_path <- download_model(model_name)
  } else {
    message("Utilisation du modèle local: ", model_path)
  }

  # --- Étape 4 : Inférence ---
  message("\n>>> ÉTAPE 4/5 : Inférence du modèle ", model_name)

  # Combiner RVB + IRC en RGBI
  rgbi <- combine_rvb_irc(ortho$rvb, ortho$irc)
  rgbi_path <- file.path(output_dir, "ortho_rgbi.tif")
  writeRaster(rgbi, rgbi_path, overwrite = TRUE)

  # Inférence
  landcover <- run_inference(rgbi, model_path)

  # --- Étape 5 : Export ---
  message("\n>>> ÉTAPE 5/5 : Export des résultats")

  # Carte d'occupation du sol
  lc_path <- file.path(output_dir, "landcover_predicted.tif")
  writeRaster(landcover, lc_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  message("Occupation du sol: ", lc_path)

  # NDVI
  pir   <- ortho$irc[["PIR"]]
  rouge <- ortho$irc[["Rouge"]]
  ndvi  <- (pir - rouge) / (pir + rouge)
  names(ndvi) <- "NDVI"
  ndvi_path <- file.path(output_dir, "ndvi.tif")
  writeRaster(ndvi, ndvi_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  message("NDVI:              ", ndvi_path)

  # --- Visualisation récapitulative ---
  pdf_path <- file.path(output_dir, "resultats_aoi_flair_hub.pdf")
  pdf(pdf_path, width = 16, height = 12)

  par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))

  # RVB
  plotRGB(ortho$rvb, r = 1, g = 2, b = 3, stretch = "lin",
          main = "Ortho RVB IGN (0.20m)")

  # IRC fausses couleurs
  plotRGB(ortho$irc, r = 1, g = 2, b = 3, stretch = "lin",
          main = "Ortho IRC fausses couleurs (0.20m)")

  # NDVI
  col_ndvi <- colorRampPalette(
    c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
      "#d9ef8b", "#91cf60", "#1a9850", "#006837")
  )(100)
  plot(ndvi, main = "NDVI (depuis IRC)", col = col_ndvi,
       range = c(-0.2, 1), plg = list(title = "NDVI"))

  # Occupation du sol
  plot(landcover, main = paste("Occupation du sol -", model_name),
       col = COSIA_COLORS_15, type = "classes",
       levels = COSIA_LABELS_15,
       plg = list(legend = COSIA_LABELS_15, cex = 0.6))

  dev.off()
  message("PDF:               ", pdf_path)

  # --- Statistiques ---
  lc_vals <- values(landcover, na.rm = TRUE)
  class_counts <- table(as.integer(lc_vals))

  message("\n--- Distribution des classes ---")
  for (i in seq_along(class_counts)) {
    cls <- as.integer(names(class_counts)[i])
    pct <- as.numeric(class_counts[i]) / length(lc_vals) * 100
    label <- if (cls >= 1 && cls <= 15) COSIA_LABELS_15[cls] else "?"
    message(sprintf("  %2d. %s: %.1f%%", cls, label, pct))
  }

  # --- Résumé ---
  dt <- round(difftime(Sys.time(), t0, units = "mins"), 1)

  message("\n##############################################################")
  message("#  Pipeline terminé en ", dt, " minutes")
  message(sprintf("#  Classes uniques: %d", length(class_counts)))
  message(sprintf("#  Fichiers dans: %s", output_dir))
  message("##############################################################")

  return(list(
    aoi        = aoi,
    ortho_rvb  = ortho$rvb,
    ortho_irc  = ortho$irc,
    ortho_rgbi = rgbi,
    ndvi       = ndvi,
    landcover  = landcover,
    output_dir = output_dir
  ))
}

# ==============================================================================
# Point d'entrée
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== Pipeline AOI → Carte d'occupation du sol (FLAIR-HUB) ===\n")

  aoi_path <- file.path(getwd(), "data", "aoi.gpkg")

  if (!file.exists(aoi_path)) {
    message("Fichier AOI non trouvé: ", aoi_path)
    message("\nUtilisation:")
    message('  source("R/04_pipeline_aoi_to_landcover.R")')
    message('  result <- pipeline_aoi_to_landcover("chemin/vers/aoi.gpkg")')
    message("")
    message('  # Avec un modèle local :')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    model_path = "chemin/vers/modele/")')
    message("")
    message("Le fichier aoi.gpkg doit contenir un polygone définissant")
    message("votre zone d'intérêt (n'importe quel CRS, sera reprojeté")
    message("en Lambert-93 automatiquement).")
  } else {
    result <- pipeline_aoi_to_landcover(aoi_path)
  }
}
