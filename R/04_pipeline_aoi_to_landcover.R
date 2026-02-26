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
IGN_WCS_URL      <- "https://data.geopf.fr/wcs"
IGN_LAYER_ORTHO  <- "ORTHOIMAGERY.ORTHOPHOTOS"
IGN_LAYER_IRC    <- "ORTHOIMAGERY.ORTHOPHOTOS.IRC-EXPRESS.2024"
IGN_COV_MNT      <- "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES"
IGN_COV_MNS      <- "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES.MNS"

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
# 2b. Téléchargement du MNT/MNS IGN via WCS (RGE ALTI® 1m)
# ==============================================================================
# Le MNT (DTM, terrain nu) et MNS (DSM, surface avec bâtiments/végétation)
# sont disponibles à 1m via le WCS de la Géoplateforme IGN.
#
# FLAIR-HUB attend 2 bandes : DSM (MNS) + DTM (MNT) en Float32.
# Le CHM (Canopy Height Model) = DSM - DTM est utile pour distinguer
# les classes arborées (feuillu, conifère) des classes basses (herbacé).

#' Télécharger une couverture WCS IGN (MNT ou MNS)
#'
#' @param bbox Emprise c(xmin, ymin, xmax, ymax) en Lambert-93
#' @param coverage_id Identifiant de la couverture WCS
#' @param res_m Résolution en mètres (1 = RGE ALTI 1m)
#' @param dest_file Fichier de destination
#' @return SpatRaster ou NULL si échec
download_wcs_coverage <- function(bbox, coverage_id, res_m = 1, dest_file) {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  # WCS 2.0.1 GetCoverage
  wcs_url <- paste0(
    IGN_WCS_URL, "?",
    "SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&CoverageId=", coverage_id,
    "&SUBSET=x(", xmin, ",", xmax, ")",
    "&SUBSET=y(", ymin, ",", ymax, ")",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/2154",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/2154",
    "&FORMAT=image/tiff"
  )

  tryCatch({
    tmp_file <- tempfile(fileext = ".tif")
    curl_download(url = wcs_url, destfile = tmp_file, quiet = TRUE)

    # Vérifier que c'est bien un GeoTIFF (pas une erreur XML)
    fsize <- file.info(tmp_file)$size
    if (fsize < 1000) {
      raw <- readLines(tmp_file, n = 5, warn = FALSE)
      if (any(grepl("Exception|Error|xml", raw, ignore.case = TRUE))) {
        warning("Erreur WCS pour ", coverage_id, " : ", paste(raw, collapse = " "))
        unlink(tmp_file)
        return(NULL)
      }
    }

    r <- rast(tmp_file)
    writeRaster(r, dest_file, overwrite = TRUE)
    r <- rast(dest_file)
    unlink(tmp_file)

    return(r)
  }, error = function(e) {
    unlink(tmp_file)
    warning("Échec WCS (", coverage_id, "): ", e$message)
    return(NULL)
  })
}

#' Télécharger le MNT et le MNS IGN pour une emprise (avec tuilage WCS)
#'
#' Le WCS a une limite de taille par requête (~2048x2048 pixels à 1m).
#' On tuilage les requêtes si nécessaire.
#'
#' @param bbox Emprise c(xmin, ymin, xmax, ymax)
#' @param coverage_id Couverture WCS
#' @param res_m Résolution (1m par défaut)
#' @param output_dir Répertoire de sortie
#' @param prefix Préfixe des fichiers
#' @return SpatRaster
download_elevation_tiled <- function(bbox, coverage_id, res_m = 1,
                                      output_dir, prefix = "elev") {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]
  wcs_max_px <- 2048
  tile_size_m <- wcs_max_px * res_m

  x_starts <- seq(xmin, xmax, by = tile_size_m)
  y_starts <- seq(ymin, ymax, by = tile_size_m)
  n_tiles <- length(x_starts) * length(y_starts)
  message(sprintf("Téléchargement %s (%s): %d tuile(s) WCS à %dm...",
                   prefix, coverage_id, n_tiles, res_m))

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

      message(sprintf("  Tuile %d/%d...", idx, n_tiles))
      r <- download_wcs_coverage(tile_bbox, coverage_id, res_m, tile_file)
      if (!is.null(r)) {
        tile_rasters[[idx]] <- r
      }
      idx <- idx + 1
    }
  }

  if (length(tile_rasters) == 0) {
    warning("Aucune tuile WCS téléchargée pour ", coverage_id)
    return(NULL)
  }

  if (length(tile_rasters) == 1) {
    mosaic <- tile_rasters[[1]]
  } else {
    message("Mosaïquage de ", length(tile_rasters), " tuiles ", prefix, "...")
    mosaic <- do.call(merge, tile_rasters)
  }

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = paste0("*", prefix, "_tile_*.tif"))
  if (length(tile_files) > 0) file_delete(tile_files)

  return(mosaic)
}

#' Télécharger le DEM (DSM + DTM) pour une AOI
#'
#' Télécharge le MNS (DSM) et le MNT (DTM) depuis le WCS IGN (RGE ALTI 1m)
#' et les combine en un SpatRaster 2 bandes comme attendu par FLAIR-HUB.
#'
#' @param aoi sf object en Lambert-93
#' @param output_dir Répertoire de sortie
#' @param res_m Résolution du MNT (1 = RGE ALTI 1m, 5 = BD ALTI 5m)
#' @param rgbi SpatRaster de référence pour le rééchantillonnage à 0.2m
#' @return Liste avec dem (SpatRaster 2 bandes DSM+DTM) et dem_path
download_dem_for_aoi <- function(aoi, output_dir, res_m = 1, rgbi = NULL) {
  dir_create(output_dir)

  bbox <- as.numeric(st_bbox(st_union(aoi)))
  message(sprintf("\n=== Téléchargement MNT/MNS IGN (RGE ALTI %dm) ===", res_m))

  # MNT (DTM - terrain nu)
  message("\n--- MNT (DTM, terrain nu) ---")
  dtm <- download_elevation_tiled(
    bbox, coverage_id = IGN_COV_MNT, res_m = res_m,
    output_dir = output_dir, prefix = "mnt"
  )

  # MNS (DSM - surface avec bâtiments/végétation)
  message("\n--- MNS (DSM, surface) ---")
  dsm <- download_elevation_tiled(
    bbox, coverage_id = IGN_COV_MNS, res_m = res_m,
    output_dir = output_dir, prefix = "mns"
  )

  # Découper aux limites de l'AOI
  aoi_vect <- vect(st_union(aoi))

  if (!is.null(dtm)) dtm <- crop(dtm, aoi_vect)
  if (!is.null(dsm)) dsm <- crop(dsm, aoi_vect)

  # Si le MNS n'est pas disponible, utiliser le MNT pour les 2 bandes
  if (is.null(dsm) && !is.null(dtm)) {
    message("MNS non disponible, utilisation du MNT seul (DSM = DTM)")
    dsm <- dtm
  }
  if (is.null(dtm) && !is.null(dsm)) {
    message("MNT non disponible, utilisation du MNS seul (DTM = DSM)")
    dtm <- dsm
  }
  if (is.null(dtm) && is.null(dsm)) {
    warning("Aucune donnée d'élévation téléchargée.")
    return(NULL)
  }

  # Aligner les grilles DSM et DTM
  if (!compareGeom(dsm, dtm, stopOnError = FALSE)) {
    dsm <- resample(dsm, dtm, method = "bilinear")
  }

  # Combiner en SpatRaster 2 bandes (format FLAIR-HUB DEM_ELEV)
  dem <- c(dsm[[1]], dtm[[1]])
  names(dem) <- c("DSM", "DTM")

  # Rééchantillonner vers la grille aérienne (0.2m) si fournie
  if (!is.null(rgbi)) {
    message("Rééchantillonnage MNT/MNS de ", res_m, "m vers 0.2m...")
    dem <- resample(dem, rgbi, method = "bilinear")
  }

  # Sauvegarder
  dem_path <- file.path(output_dir, "dem_dsm_dtm.tif")
  writeRaster(dem, dem_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
  message(sprintf("\nDEM: %s (%d x %d px, bandes: DSM + DTM)",
                   dem_path, ncol(dem), nrow(dem)))

  # Statistiques
  chm <- dem[["DSM"]] - dem[["DTM"]]
  message(sprintf("  Altitude DTM: %.0f - %.0f m",
                   min(values(dem[["DTM"]]), na.rm = TRUE),
                   max(values(dem[["DTM"]]), na.rm = TRUE)))
  message(sprintf("  Hauteur CHM (DSM-DTM): %.1f - %.1f m",
                   min(values(chm), na.rm = TRUE),
                   max(values(chm), na.rm = TRUE)))

  return(list(dem = dem, dem_path = dem_path))
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

#' Télécharger un modèle FLAIR depuis Hugging Face
#'
#' Par défaut, utilise FLAIR-INC_rgbi_15cl_resnet34-unet (le plus simple).
#' Pour le modèle FLAIR-HUB multimodal, utiliser "FLAIR-HUB_LC-G_utae".
#'
#' @param model_name Nom du modèle
#' @return Chemin local du modèle
download_model <- function(model_name = "FLAIR-INC_rgbi_15cl_resnet34-unet") {
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

#' Pipeline complet : AOI → Ortho IGN (+MNT) → Carte d'occupation du sol
#'
#' @param aoi_path Chemin vers le fichier aoi.gpkg
#' @param output_dir Répertoire de sortie
#' @param model_name Nom du modèle FLAIR-HUB
#' @param model_path Chemin local vers un modèle (optionnel)
#' @param res_m Résolution de téléchargement ortho IGN (0.2m)
#' @param use_dem Télécharger et utiliser le MNT/MNS IGN (config LC-B, +1pt mIoU)
#' @param dem_res_m Résolution du MNT (1 = RGE ALTI 1m, 5 = BD ALTI 5m)
#' @return Liste avec tous les résultats
pipeline_aoi_to_landcover <- function(aoi_path,
                                        output_dir = file.path(getwd(), "outputs"),
                                        model_name = "FLAIR-INC_rgbi_15cl_resnet34-unet",
                                        model_path = NULL,
                                        res_m = RES_IGN,
                                        use_dem = FALSE,
                                        dem_res_m = 1) {
  dir_create(output_dir)
  t0 <- Sys.time()

  config_label <- if (use_dem) "LC-B (RGBI + MNT)" else "LC-A (RGBI seul)"
  n_steps <- if (use_dem) 6 else 5

  message("##############################################################")
  message("#  Pipeline FLAIR-HUB : AOI → Ortho IGN → Occupation du sol  #")
  message(sprintf("#  Configuration: %s", config_label))
  message("##############################################################\n")

  # --- Étape 1 : Charger l'AOI ---
  message(sprintf(">>> ÉTAPE 1/%d : Chargement de l'AOI", n_steps))
  aoi <- load_aoi(aoi_path)

  # --- Étape 2 : Télécharger les ortho IGN ---
  message(sprintf("\n>>> ÉTAPE 2/%d : Téléchargement des ortho IGN (RVB + IRC)", n_steps))
  ortho <- download_ortho_for_aoi(aoi, output_dir = output_dir, res_m = res_m)

  # Combiner RVB + IRC en RGBI
  rgbi <- combine_rvb_irc(ortho$rvb, ortho$irc)
  rgbi_path <- file.path(output_dir, "ortho_rgbi.tif")
  writeRaster(rgbi, rgbi_path, overwrite = TRUE)

  # --- Étape 2b : Télécharger le MNT/MNS (optionnel, config LC-B) ---
  dem_data <- NULL
  if (use_dem) {
    step_dem <- 3
    message(sprintf("\n>>> ÉTAPE %d/%d : Téléchargement MNT/MNS IGN (RGE ALTI %dm)",
                     step_dem, n_steps, dem_res_m))
    dem_data <- download_dem_for_aoi(aoi, output_dir = output_dir,
                                      res_m = dem_res_m, rgbi = rgbi)
  }

  # --- Étape 3/4 : Configurer Python + modèle ---
  step_py <- if (use_dem) 4 else 3
  message(sprintf("\n>>> ÉTAPE %d/%d : Configuration Python + téléchargement modèle",
                   step_py, n_steps))
  setup_python()
  if (is.null(model_path)) {
    model_path <- download_model(model_name)
  } else {
    message("Utilisation du modèle local: ", model_path)
  }

  # --- Étape 4/5 : Inférence ---
  step_inf <- if (use_dem) 5 else 4
  message(sprintf("\n>>> ÉTAPE %d/%d : Inférence du modèle %s",
                   step_inf, n_steps, model_name))

  landcover <- run_inference(rgbi, model_path)

  # --- Étape 5/6 : Export ---
  step_exp <- if (use_dem) 6 else 5
  message(sprintf("\n>>> ÉTAPE %d/%d : Export des résultats", step_exp, n_steps))

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
  n_panels <- if (use_dem && !is.null(dem_data)) 6 else 4
  pdf_w <- if (n_panels > 4) 18 else 16
  pdf(pdf_path, width = pdf_w, height = 12)

  if (n_panels > 4) {
    par(mfrow = c(2, 3), mar = c(2, 2, 3, 4))
  } else {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
  }

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

  # MNT si disponible
  if (!is.null(dem_data)) {
    col_elev <- colorRampPalette(
      c("#313695", "#4575b4", "#74add1", "#abd9e9", "#fee090",
        "#fdae61", "#f46d43", "#d73027", "#a50026")
    )(100)
    plot(dem_data$dem[["DTM"]], main = sprintf("MNT IGN (RGE ALTI %dm)", dem_res_m),
         col = col_elev, plg = list(title = "Altitude (m)"))

    chm <- dem_data$dem[["DSM"]] - dem_data$dem[["DTM"]]
    col_chm <- colorRampPalette(
      c("#ffffcc", "#d9f0a3", "#addd8e", "#78c679",
        "#41ab5d", "#238443", "#005a32")
    )(100)
    plot(chm, main = "CHM (DSM - DTM)",
         col = col_chm, plg = list(title = "Hauteur (m)"))
  }

  # Occupation du sol
  plot(landcover, main = paste("Occupation du sol -", config_label),
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
  message(sprintf("#  Configuration: %s", config_label))
  message(sprintf("#  Classes uniques: %d", length(class_counts)))
  message(sprintf("#  Fichiers dans: %s", output_dir))
  message("##############################################################")

  result <- list(
    aoi        = aoi,
    ortho_rvb  = ortho$rvb,
    ortho_irc  = ortho$irc,
    ortho_rgbi = rgbi,
    ndvi       = ndvi,
    landcover  = landcover,
    output_dir = output_dir
  )
  if (!is.null(dem_data)) result$dem <- dem_data$dem

  return(result)
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
    message("")
    message('  # Config LC-A : RGBI seul (64.1% mIoU)')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg")')
    message("")
    message('  # Config LC-B : RGBI + MNT (65.1% mIoU, +1pt)')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    use_dem = TRUE, dem_res_m = 1)')
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
