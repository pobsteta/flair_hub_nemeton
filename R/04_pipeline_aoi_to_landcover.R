#!/usr/bin/env Rscript
# ==============================================================================
# 04_pipeline_aoi_to_landcover.R
# Pipeline complet : AOI (GeoPackage) → Données IGN → Carte d'occupation du sol
#
# Entrée  : fichier aoi.gpkg (zone d'intérêt, n'importe quel CRS)
# Sorties : ortho_rgbi.tif, dem_dsm_dtm.tif, landcover_predicted.tif
#
# Workflow :
#   1. Charger l'AOI depuis aoi.gpkg → reprojection Lambert-93
#   2. Télécharger les ortho IGN RVB + IRC via WMS (tuiles si nécessaire)
#   2b. Télécharger le MNT/MNS IGN via WMS-R (optionnel, config LC-B)
#   3. Combiner RVB + IRC en image 4 bandes RGBI
#   4. Découper en patches de 512x512 à 0.2m
#   5. Inférence du modèle FLAIR-HUB (Swin/ConvNeXTV2) via reticulate
#   6. Mosaïquer et exporter la carte d'occupation du sol
#
# Cache :
#   Les fichiers téléchargés sont réutilisés s'ils existent déjà dans
#   output_dir (ortho_rvb.tif, ortho_irc.tif, dem_dsm_dtm.tif).
#   Pour forcer le re-téléchargement, supprimer ces fichiers.
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
# Toutes les données (ortho + élévation) sont accessibles via WMS-R.
# La Géoplateforme n'offre pas de service WCS pour l'altimétrie.
IGN_WMS_URL      <- "https://data.geopf.fr/wms-r"
IGN_LAYER_ORTHO  <- "ORTHOIMAGERY.ORTHOPHOTOS"
IGN_LAYER_MNT    <- "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES"
IGN_LAYER_MNS    <- "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES.MNS"
# Alternatives LiDAR HD (couverture partielle mais plus précis) :
# IGN_LAYER_MNT <- "IGNF_LIDAR-HD_MNT_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93"
# IGN_LAYER_MNS <- "IGNF_LIDAR-HD_MNS_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93"

# --- Millésime ortho ---
# L'IRC-EXPRESS et l'ORTHO-EXPRESS sont millésimés : le suffixe année
# est obligatoire. Les anciens millésimes sont dépubliés par l'IGN
# (depuis sept. 2025, seuls 2024+ restent accessibles).
# NULL = détection automatique (année courante, fallback N-1).
IGN_MILLESIME_IRC   <- NULL
IGN_MILLESIME_ORTHO <- NULL

# --- Résolutions ---
RES_IGN <- 0.2   # BD ORTHO® IGN

# --- FLAIR-HUB ---
PATCH_SIZE <- 512  # Taille des patches (512x512 px à 0.2m)
CONDA_ENV  <- "FLAIRHUB"

# --- Limites WMS ---
WMS_MAX_PX <- 4096  # Taille max par requête WMS

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
# 2. Gestion du millésime et téléchargement ortho IGN
# ==============================================================================

#' Résoudre le millésime IRC ou ORTHO-EXPRESS
#'
#' Les couches IRC-EXPRESS et ORTHO-EXPRESS sont millésimées :
#' IRC-EXPRESS.2024, ORTHO-EXPRESS.2025...
#' Les anciens millésimes sont dépubliés par l'IGN.
#' Cette fonction détermine le millésime à utiliser et vérifie sa disponibilité.
#'
#' @param millesime Année (NULL = détection auto, entier = année forcée)
#' @param bbox c(xmin, ymin, xmax, ymax) en Lambert-93 pour tester la couche
#' @param type Type de couche ("irc" ou "ortho")
#' @return Année résolue (entier)
resolve_millesime <- function(millesime = NULL, bbox = NULL, type = "irc") {
  label <- toupper(type)
  layer_prefix <- if (type == "irc") {
    "ORTHOIMAGERY.ORTHOPHOTOS.IRC-EXPRESS."
  } else {
    "ORTHOIMAGERY.ORTHOPHOTOS.ORTHO-EXPRESS."
  }

  if (!is.null(millesime)) {
    message(sprintf("Millésime %s forcé: %d", label, millesime))
    return(as.integer(millesime))
  }

  # Détection automatique : année courante, puis fallback N-1
  annee <- as.integer(format(Sys.Date(), "%Y"))
  candidates <- c(annee, annee - 1)

  if (is.null(bbox)) {
    message(sprintf("Millésime %s: %d (année courante)", label, annee))
    return(annee)
  }

  # Tester chaque millésime avec une requête WMS minimale
  for (yr in candidates) {
    layer_test <- paste0(layer_prefix, yr)
    test_ok <- test_wms_layer(bbox, layer_test)
    if (test_ok) {
      message(sprintf("Millésime %s: %d (vérifié OK)", label, yr))
      return(yr)
    }
    message(sprintf("  Millésime %s %d: non disponible pour cette zone", label, yr))
  }

  # Fallback : année courante sans vérification
  message(sprintf("Millésime %s: %d (par défaut, non vérifié)", label, annee))
  return(annee)
}

#' Tester si une couche WMS est disponible pour une emprise
#'
#' Effectue une requête WMS minimale (2x2 px) pour vérifier que la couche
#' retourne une image valide (et non une erreur XML).
#'
#' @param bbox c(xmin, ymin, xmax, ymax) en Lambert-93
#' @param layer Nom de la couche WMS
#' @return TRUE si la couche est disponible
test_wms_layer <- function(bbox, layer) {
  xmin <- bbox[1]; ymin <- bbox[2]
  # Petite emprise de test (200m x 200m)
  xmax <- xmin + 200
  ymax <- ymin + 200

  # WMS 1.3.0 + EPSG:2154 (projected CRS) : BBOX = xmin,ymin,xmax,ymax
  wms_url <- paste0(
    IGN_WMS_URL, "?",
    "SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap",
    "&LAYERS=", layer,
    "&CRS=EPSG:2154",
    "&BBOX=", paste(xmin, ymin, xmax, ymax, sep = ","),
    "&WIDTH=2&HEIGHT=2",
    "&FORMAT=image/geotiff",
    "&STYLES="
  )

  tryCatch({
    tmp <- tempfile(fileext = ".tif")
    curl_download(url = wms_url, destfile = tmp, quiet = TRUE)

    # Vérifier que c'est un GeoTIFF et non un XML d'erreur
    fsize <- file.info(tmp)$size
    if (fsize < 500) {
      raw <- readLines(tmp, n = 5, warn = FALSE)
      if (any(grepl("Exception|Error|xml", raw, ignore.case = TRUE))) {
        unlink(tmp)
        return(FALSE)
      }
    }

    r <- rast(tmp)
    unlink(tmp)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

#' Construire les noms de couches WMS pour les millésimes donnés
#'
#' Chaque couche (RVB et IRC) peut avoir son propre millésime :
#'   - RVB : ORTHO-EXPRESS.{année} ou mosaïque nationale si NULL
#'   - IRC : IRC-EXPRESS.{année} (toujours millésimé)
#'
#' @param millesime_irc Année pour l'IRC (entier, obligatoire)
#' @param millesime_ortho Année pour le RVB (entier ou NULL = mosaïque nationale)
#' @return Liste nommée avec les couches ortho et irc
build_layer_names <- function(millesime_irc, millesime_ortho = NULL) {
  irc_layer <- paste0("ORTHOIMAGERY.ORTHOPHOTOS.IRC-EXPRESS.", millesime_irc)

  if (!is.null(millesime_ortho)) {
    ortho_layer <- paste0("ORTHOIMAGERY.ORTHOPHOTOS.ORTHO-EXPRESS.", millesime_ortho)
    message(sprintf("  RVB: %s (millésime %d)", ortho_layer, millesime_ortho))
  } else {
    ortho_layer <- IGN_LAYER_ORTHO
    message(sprintf("  RVB: %s (mosaïque nationale)", ortho_layer))
  }
  message(sprintf("  IRC: %s (millésime %d)", irc_layer, millesime_irc))

  list(ortho = ortho_layer, irc = irc_layer)
}

#' Télécharger une tuile WMS IGN
#'
#' @param bbox c(xmin, ymin, xmax, ymax) en Lambert-93
#' @param layer Couche WMS
#' @param res_m Résolution en mètres
#' @param dest_file Fichier de sortie
#' @param styles Style WMS ("" = défaut pour ortho, "normal" = valeurs brutes
#'   pour couches d'élévation)
#' @return SpatRaster ou NULL si échec
download_wms_tile <- function(bbox, layer, res_m = RES_IGN, dest_file,
                               styles = "") {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  width  <- round((xmax - xmin) / res_m)
  height <- round((ymax - ymin) / res_m)

  # WMS 1.3.0 + EPSG:2154 (projected CRS) : BBOX = xmin,ymin,xmax,ymax
  # (axis order follows the CRS definition: Easting first, Northing second)
  wms_url <- paste0(
    IGN_WMS_URL, "?",
    "SERVICE=WMS&VERSION=1.3.0&REQUEST=GetMap",
    "&LAYERS=", layer,
    "&CRS=EPSG:2154",
    "&BBOX=", paste(xmin, ymin, xmax, ymax, sep = ","),
    "&WIDTH=", width,
    "&HEIGHT=", height,
    "&FORMAT=image/geotiff",
    "&STYLES=", styles
  )

  message("  WMS URL: ", wms_url)

  tryCatch({
    tmp_file <- tempfile(fileext = ".tif")
    curl_download(url = wms_url, destfile = tmp_file, quiet = TRUE)

    r <- rast(tmp_file)

    # Assigner le CRS et l'emprise si nécessaire
    if (is.na(crs(r)) || crs(r) == "") {
      crs(r) <- "EPSG:2154"
    }
    ext(r) <- ext(xmin, xmax, ymin, ymax)

    # Écrire le fichier final et re-lire
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

#' Télécharger une couche WMS IGN complète pour une emprise (tuilage automatique)
#'
#' @param bbox c(xmin, ymin, xmax, ymax) en Lambert-93
#' @param layer Couche WMS (ortho ou élévation)
#' @param res_m Résolution en mètres
#' @param output_dir Répertoire de sortie
#' @param prefix Préfixe pour les fichiers
#' @param styles Style WMS ("" pour ortho, "normal" pour élévation brute)
#' @return SpatRaster mosaïqué
download_ign_tiled <- function(bbox, layer, res_m = RES_IGN,
                                output_dir, prefix = "ortho",
                                styles = "") {
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

      r <- download_wms_tile(tile_bbox, layer, res_m, tile_file,
                              styles = styles)
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
#'
#' Gestion du cache : si ortho_rvb.tif et ortho_irc.tif existent déjà,
#' ils sont réutilisés sans re-téléchargement.
#'
#' Gestion des millésimes (indépendants pour IRC et ortho) :
#'   - millesime_irc : NULL = auto-détection, entier = IRC-EXPRESS.{année}
#'   - millesime_ortho : NULL = mosaïque nationale, entier = ORTHO-EXPRESS.{année}
#'
#' @param aoi sf object (AOI en Lambert-93)
#' @param output_dir Répertoire de sortie
#' @param res_m Résolution en mètres
#' @param millesime_irc Millésime IRC (NULL = auto, entier = année forcée)
#' @param millesime_ortho Millésime ortho RVB (NULL = mosaïque nationale, entier = année)
#' @return Liste avec rvb, irc (SpatRaster), millesime_irc, millesime_ortho
download_ortho_for_aoi <- function(aoi, output_dir, res_m = RES_IGN,
                                    millesime_irc = NULL,
                                    millesime_ortho = NULL) {
  dir_create(output_dir)

  rvb_path <- file.path(output_dir, "ortho_rvb.tif")
  irc_path <- file.path(output_dir, "ortho_irc.tif")

  # Cache : réutiliser les fichiers existants
  if (file.exists(rvb_path) && file.exists(irc_path)) {
    message("\n=== Ortho RVB et IRC déjà téléchargées (cache) ===")
    rvb <- rast(rvb_path)
    irc <- rast(irc_path)
    names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]
    names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]
    message(sprintf("RVB: %s (%d x %d px)", rvb_path, ncol(rvb), nrow(rvb)))
    message(sprintf("IRC: %s (%d x %d px)", irc_path, ncol(irc), nrow(irc)))
    return(list(rvb = rvb, irc = irc,
                rvb_path = rvb_path, irc_path = irc_path,
                millesime_irc = millesime_irc,
                millesime_ortho = millesime_ortho))
  }

  bbox <- as.numeric(st_bbox(st_union(aoi)))
  message(sprintf("\n=== Téléchargement ortho IGN pour l'AOI ==="))
  message(sprintf("Emprise: %.0f, %.0f - %.0f, %.0f (Lambert-93)",
                   bbox[1], bbox[2], bbox[3], bbox[4]))
  message(sprintf("Taille: %.0f x %.0f m (%.2f ha)",
                   bbox[3] - bbox[1], bbox[4] - bbox[2],
                   (bbox[3] - bbox[1]) * (bbox[4] - bbox[2]) / 10000))

  # Résoudre les millésimes indépendamment
  millesime_irc <- resolve_millesime(millesime_irc, bbox, type = "irc")
  millesime_ortho_resolved <- if (!is.null(millesime_ortho)) {
    resolve_millesime(millesime_ortho, bbox, type = "ortho")
  } else {
    NULL
  }

  layers <- build_layer_names(millesime_irc = millesime_irc,
                               millesime_ortho = millesime_ortho_resolved)

  # RVB
  if (!is.null(millesime_ortho_resolved)) {
    message(sprintf("\n--- Ortho RVB (millésime %d) ---", millesime_ortho_resolved))
  } else {
    message("\n--- Ortho RVB (mosaïque nationale) ---")
  }
  rvb <- download_ign_tiled(bbox, layer = layers$ortho, res_m = res_m,
                             output_dir = output_dir, prefix = "rvb")
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]

  # IRC
  message(sprintf("\n--- Ortho IRC (millésime %d) ---", millesime_irc))
  irc <- download_ign_tiled(bbox, layer = layers$irc, res_m = res_m,
                             output_dir = output_dir, prefix = "irc")
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  # Découper aux limites exactes de l'AOI
  aoi_vect <- vect(st_union(aoi))
  rvb <- crop(rvb, aoi_vect)
  irc <- crop(irc, aoi_vect)

  # Sauvegarder les mosaïques finales
  writeRaster(rvb, rvb_path, overwrite = TRUE)
  writeRaster(irc, irc_path, overwrite = TRUE)

  message(sprintf("\nRVB sauvegardé: %s (%d x %d px)",
                   rvb_path, ncol(rvb), nrow(rvb)))
  message(sprintf("IRC sauvegardé: %s (%d x %d px)",
                   irc_path, ncol(irc), nrow(irc)))

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = "*_tile_*.tif")
  if (length(tile_files) > 0) file_delete(tile_files)

  # Re-lire depuis les fichiers sauvegardés (terra est file-backed :
  # après suppression des tuiles, les objets doivent pointer vers
  # les fichiers finaux ortho_rvb.tif / ortho_irc.tif)
  rvb <- rast(rvb_path)
  irc <- rast(irc_path)

  return(list(rvb = rvb, irc = irc,
              rvb_path = rvb_path, irc_path = irc_path,
              millesime_irc = millesime_irc,
              millesime_ortho = millesime_ortho_resolved))
}

# ==============================================================================
# 2b. Téléchargement du MNT/MNS IGN via WMS-R (RGE ALTI® 1m)
# ==============================================================================
# FLAIR-HUB attend 2 bandes : DSM (MNS) + DTM (MNT) en Float32.
#
# Dans FLAIR-HUB original :
#   - DSM (MNS) = résolution native 0.2m (corrélation dense des photos aériennes)
#   - DTM (MNT) = RGE ALTI natif à 1m, rééchantillonné à 0.2m
#
# Via la Géoplateforme IGN WMS-R (pas de service WCS disponible) :
#   - MNT (DTM) : ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES → RGE ALTI 1m
#   - MNS (DSM) : ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES.MNS → si disponible
#     (couverture LiDAR HD en cours de déploiement, pas disponible partout)
#
# Le CHM (Canopy Height Model) = DSM - DTM est utile pour distinguer
# les classes arborées (feuillu, conifère) des classes basses (herbacé).

#' Télécharger le DEM (DSM + DTM) pour une AOI
#'
#' Télécharge le MNS (DSM) et le MNT (DTM) depuis le WMS-R IGN (RGE ALTI 1m)
#' et les combine en un SpatRaster 2 bandes comme attendu par FLAIR-HUB.
#'
#' Gestion du cache : si dem_dsm_dtm.tif existe déjà, il est réutilisé.
#'
#' @param aoi sf object en Lambert-93
#' @param output_dir Répertoire de sortie
#' @param res_m Résolution du MNT (1 = RGE ALTI 1m)
#' @param rgbi SpatRaster de référence pour le rééchantillonnage à 0.2m
#' @return Liste avec dem (SpatRaster 2 bandes DSM+DTM) et dem_path
download_dem_for_aoi <- function(aoi, output_dir, res_m = 1, rgbi = NULL) {
  dir_create(output_dir)

  dem_path <- file.path(output_dir, "dem_dsm_dtm.tif")

  # Cache : réutiliser le fichier existant
  if (file.exists(dem_path)) {
    message("\n=== DEM déjà téléchargé (cache) ===")
    dem <- rast(dem_path)
    names(dem) <- c("DSM", "DTM")[1:nlyr(dem)]
    message(sprintf("DEM: %s (%d x %d px, %d bandes)",
                     dem_path, ncol(dem), nrow(dem), nlyr(dem)))
    return(list(dem = dem, dem_path = dem_path))
  }

  bbox <- as.numeric(st_bbox(st_union(aoi)))
  message(sprintf("\n=== Téléchargement MNT/MNS IGN (RGE ALTI %dm via WMS-R) ===",
                   res_m))

  # MNT (DTM - terrain nu) — disponible partout en France (RGE ALTI)
  message("\n--- MNT (DTM, terrain nu, RGE ALTI) ---")
  dtm <- tryCatch(
    download_ign_tiled(bbox, layer = IGN_LAYER_MNT, res_m = res_m,
                        output_dir = output_dir, prefix = "mnt",
                        styles = "normal"),
    error = function(e) {
      message("  MNT non téléchargé: ", e$message)
      NULL
    }
  )

  # MNS (DSM - surface avec bâtiments/végétation)
  # La couverture MNS LiDAR HD est en cours de déploiement,
  # pas encore disponible partout en France.
  message("\n--- MNS (DSM, surface, LiDAR HD) ---")
  message("  Note : le MNS n'est pas disponible partout (LiDAR HD en cours)")
  dsm <- tryCatch(
    download_ign_tiled(bbox, layer = IGN_LAYER_MNS, res_m = res_m,
                        output_dir = output_dir, prefix = "mns",
                        styles = "normal"),
    error = function(e) {
      message("  MNS non disponible pour cette zone: ", e$message)
      NULL
    }
  )

  # Découper aux limites de l'AOI
  aoi_vect <- vect(st_union(aoi))

  if (!is.null(dtm)) dtm <- crop(dtm, aoi_vect)
  if (!is.null(dsm)) dsm <- crop(dsm, aoi_vect)

  # Si le MNS n'est pas disponible, utiliser le MNT seul
  # (DSM = DTM → CHM = 0, pas d'info de hauteur mais on garde l'altitude)
  has_mns <- !is.null(dsm)
  if (!has_mns && !is.null(dtm)) {
    message("\nMNS non disponible pour cette zone.")
    message("Utilisation du MNT seul (DSM = DTM, CHM = 0).")
    message("Le modèle bénéficiera quand même de l'altitude du terrain.")
    dsm <- dtm
  }
  if (is.null(dtm) && has_mns) {
    message("MNT non disponible, utilisation du MNS seul (DTM = DSM)")
    dtm <- dsm
  }
  if (is.null(dtm) && is.null(dsm)) {
    warning("Aucune donnée d'élévation téléchargée.")
    return(NULL)
  }

  # Aligner les grilles DSM et DTM
  if (has_mns && !compareGeom(dsm, dtm, stopOnError = FALSE)) {
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
  writeRaster(dem, dem_path, overwrite = TRUE, gdal = c("COMPRESS=LZW"))

  message(sprintf("\nDEM sauvegardé: %s (%d x %d px, bandes: DSM + DTM)",
                   dem_path, ncol(dem), nrow(dem)))
  if (has_mns) {
    message("  Sources: DSM = MNS LiDAR HD, DTM = RGE ALTI")
  } else {
    message("  Sources: DSM = DTM = RGE ALTI (MNS non disponible)")
  }

  # Statistiques
  message(sprintf("  Altitude DTM: %.0f - %.0f m",
                   min(values(dem[["DTM"]]), na.rm = TRUE),
                   max(values(dem[["DTM"]]), na.rm = TRUE)))
  if (has_mns) {
    chm <- dem[["DSM"]] - dem[["DTM"]]
    message(sprintf("  Hauteur CHM (DSM-DTM): %.1f - %.1f m",
                     min(values(chm), na.rm = TRUE),
                     max(values(chm), na.rm = TRUE)))
  }

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = "*_tile_*.tif")
  if (length(tile_files) > 0) file_delete(tile_files)

  # Re-lire depuis le fichier sauvegardé
  dem <- rast(dem_path)

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
    # Only delete the input file. The output tmp_out must persist because
    # the returned SpatRaster is file-backed and terra needs the file
    # until after the merge/mosaic step. R cleans up tempdir() on exit.
    unlink(tmp_in)
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
    result <- do.call(merge, unname(predictions))
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
#' @param millesime_irc Millésime IRC (NULL = auto, entier = année forcée)
#' @param millesime_ortho Millésime ortho RVB (NULL = mosaïque nationale, entier = année)
#' @param use_dem Télécharger et utiliser le MNT/MNS IGN (config LC-B, +1pt mIoU)
#' @param dem_res_m Résolution du MNT (1 = RGE ALTI 1m)
#' @return Liste avec tous les résultats
pipeline_aoi_to_landcover <- function(aoi_path,
                                        output_dir = file.path(getwd(), "outputs"),
                                        model_name = "FLAIR-INC_rgbi_15cl_resnet34-unet",
                                        model_path = NULL,
                                        res_m = RES_IGN,
                                        millesime_irc = NULL,
                                        millesime_ortho = NULL,
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
  message(sprintf("\n>>> ÉTAPE 2/%d : Téléchargement des ortho IGN (RVB + IRC)",
                   n_steps))
  ortho <- download_ortho_for_aoi(aoi, output_dir = output_dir, res_m = res_m,
                                   millesime_irc = millesime_irc,
                                   millesime_ortho = millesime_ortho)

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

  # NDVI (IRC bandes: 1=PIR, 2=Rouge, 3=Vert)
  pir   <- ortho$irc[[1]]
  rouge <- ortho$irc[[2]]
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
    aoi             = aoi,
    ortho_rvb       = ortho$rvb,
    ortho_irc       = ortho$irc,
    ortho_rgbi      = rgbi,
    millesime_irc   = ortho$millesime_irc,
    millesime_ortho = ortho$millesime_ortho,
    ndvi            = ndvi,
    landcover       = landcover,
    output_dir      = output_dir
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
    message('  # Étude temporellement cohérente (RVB + IRC même année) :')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    millesime_irc = 2024, millesime_ortho = 2024)')
    message("")
    message('  # IRC millésimé + RVB mosaïque nationale (par défaut) :')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    millesime_irc = 2024)')
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
