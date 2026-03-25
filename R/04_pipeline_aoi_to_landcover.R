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
#   4. Découper en patches de 512x512 à 0.2m (overlap ajusté au buffer)
#   5. Inférence du modèle FLAIR-HUB (Swin/ConvNeXTV2) via reticulate
#      → buffer : rognage des bords de chaque patch (partie centrale conservée)
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
IGN_LAYER_IRC    <- "ORTHOIMAGERY.ORTHOPHOTOS.IRC"
IGN_LAYER_MNT    <- "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES"
IGN_LAYER_MNS    <- "ELEVATION.ELEVATIONGRIDCOVERAGE.HIGHRES.MNS"
# Alternatives LiDAR HD (couverture partielle mais plus précis) :
# IGN_LAYER_MNT <- "IGNF_LIDAR-HD_MNT_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93"
# IGN_LAYER_MNS <- "IGNF_LIDAR-HD_MNS_ELEVATION.ELEVATIONGRIDCOVERAGE.LAMB93"

# --- Millésimes ortho (NULL = mosaïque la plus récente) ---
MILLESIME_ORTHO  <- NULL
MILLESIME_IRC    <- NULL

# --- Résolutions ---
RES_IGN <- 0.2   # BD ORTHO® IGN

# --- FLAIR-HUB ---
PATCH_SIZE <- 512  # Taille des patches (512x512 px à 0.2m)
CONDA_ENV  <- "FLAIRHUB"

# --- Limites WMS ---
WMS_MAX_PX <- 4096  # Taille max par requête WMS

# --- Buffer / overlap pour le blending (inférence) ---
# Contrôle l'overlap entre patches : overlap = max(64, 2 × BUFFER_PX).
# Le mosaïquage utilise une fenêtre de Hann (cosinus 2D) : chaque pixel
# reçoit la prédiction du patch dont le centre est le plus proche.
# Cela élimine naturellement les artefacts de bord du réseau de neurones
# et produit une mosaïque sans couture visible.
#   - 0   : overlap minimal de 64 px (blending léger)
#   - 32  : overlap = 64 px  (identique à 0)
#   - 64  : recommandé (overlap = 128 px, bon compromis qualité/vitesse)
#   - 128 : agressif (overlap = 256 px, meilleur blending mais plus lent)
BUFFER_PX <- 128

# --- Registre des modèles FLAIR supportés ---
# Chaque entrée décrit l'architecture du modèle pour l'instanciation automatique.
# Le champ 'decoder' correspond à la classe smp (segmentation_models_pytorch).
# Le champ 'encoder' correspond au nom timm / smp de l'encodeur.
FLAIR_MODELS <- list(
  "FLAIR-INC_rgbi_15cl_resnet34-unet" = list(
    hf_repo = "IGNF/FLAIR-INC_rgbi_15cl_resnet34-unet",
    encoder = "resnet34", decoder = "Unet",
    in_channels = 4, n_out = 19, n_classes = 15
  ),
  "FLAIR-INC_rgb_15cl_resnet34-unet" = list(
    hf_repo = "IGNF/FLAIR-INC_rgb_15cl_resnet34-unet",
    encoder = "resnet34", decoder = "Unet",
    in_channels = 3, n_out = 19, n_classes = 15
  ),
  "FLAIR-INC_rgbie_15cl_resnet34-unet" = list(
    hf_repo = "IGNF/FLAIR-INC_rgbie_15cl_resnet34-unet",
    encoder = "resnet34", decoder = "Unet",
    in_channels = 5, n_out = 19, n_classes = 15
  ),
  "FLAIR-INC_rgb_12cl_resnet34-unet" = list(
    hf_repo = "IGNF/FLAIR-INC_rgb_12cl_resnet34-unet",
    encoder = "resnet34", decoder = "Unet",
    in_channels = 3, n_out = 19, n_classes = 12
  )
)

# --- Classes CoSIA (palette officielle FLAIR-HUB) ---
COSIA_LABELS_15 <- c(
  "Bâtiment", "Serre", "Piscine", "Imperméable", "Perméable",
  "Sol nu", "Eau", "Neige", "Herbacé", "Agricole",
  "Labouré", "Vigne", "Feuillu", "Conifère", "Lande"
)

# Couleurs FLAIR-1 remappées vers l'ordre CoSIA thématique
COSIA_COLORS_15 <- c(
  "#db0e9a",  #  1 Bâtiment    (FLAIR-1: building)
  "#9999ff",  #  2 Serre       (FLAIR-1: greenhouse)
  "#3de6eb",  #  3 Piscine     (FLAIR-1: swimming pool)
  "#f80c00",  #  4 Imperméable (FLAIR-1: impervious)
  "#938e7b",  #  5 Perméable   (FLAIR-1: pervious)
  "#a97101",  #  6 Sol nu      (FLAIR-1: bare soil)
  "#1553ae",  #  7 Eau         (FLAIR-1: water)
  "#ffffff",  #  8 Neige       (FLAIR-1: snow)
  "#55ff00",  #  9 Herbacé     (FLAIR-1: herbaceous)
  "#fff30d",  # 10 Agricole    (FLAIR-1: agricultural)
  "#e4df7c",  # 11 Labouré     (FLAIR-1: plowed)
  "#660082",  # 12 Vigne       (FLAIR-1: vineyard)
  "#46e483",  # 13 Feuillu     (FLAIR-1: deciduous)
  "#194a26",  # 14 Conifère    (FLAIR-1: coniferous)
  "#f3a60d"   # 15 Lande       (FLAIR-1: brushwood)
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

#' Construire le nom de couche WMS IGN selon le millésime
#'
#' @param type "ortho" ou "irc"
#' @param millesime NULL (mosaïque la plus récente) ou entier (ex: 2024)
#' @return Nom de couche WMS
ign_layer_name <- function(type = c("ortho", "irc"), millesime = NULL) {
  type <- match.arg(type)
  if (is.null(millesime)) {
    if (type == "ortho") return(IGN_LAYER_ORTHO)
    else                 return(IGN_LAYER_IRC)
  }
  millesime <- as.character(millesime)
  if (type == "ortho") {
    return(paste0("ORTHOIMAGERY.ORTHOPHOTOS", millesime))
  } else {
    return(paste0("ORTHOIMAGERY.ORTHOPHOTOS.IRC.", millesime))
  }
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
                               styles = "", max_retries = 3) {
  xmin <- bbox[1]; ymin <- bbox[2]; xmax <- bbox[3]; ymax <- bbox[4]

  width  <- round((xmax - xmin) / res_m)
  height <- round((ymax - ymin) / res_m)

  # WMS 1.3.0 avec CRS EPSG:2154 : BBOX = xmin,ymin,xmax,ymax
  # (EPSG:2154 est un CRS projeté avec axes Easting, Northing)
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

  # Handle curl avec HTTP/1.1 forcé (évite les erreurs HTTP/2 du serveur IGN)
  h <- curl::new_handle()
  curl::handle_setopt(h, http_version = 2L)  # CURL_HTTP_VERSION_1_1

  tmp_file <- tempfile(fileext = ".tif")

  for (attempt in seq_len(max_retries)) {
    result <- tryCatch({
      curl::curl_download(url = wms_url, destfile = tmp_file,
                          quiet = TRUE, handle = h)

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
      e
    })

    # Si succès (SpatRaster retourné), on sort
    if (!inherits(result, "error")) return(result)

    # Sinon retry avec backoff exponentiel
    if (attempt < max_retries) {
      wait_s <- 2^attempt
      message(sprintf("  Retry %d/%d dans %ds (%s)",
                       attempt, max_retries, wait_s, result$message))
      Sys.sleep(wait_s)
    } else {
      unlink(tmp_file)
      warning("Échec WMS après ", max_retries, " tentatives: ", result$message)
      return(NULL)
    }
  }
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

#' Vérifier qu'un raster WMS contient des données réelles
#'
#' Certaines couches millésimées ne couvrent pas toutes les zones.
#' Le WMS retourne alors un raster valide mais vide (pixels à 0 ou NA).
#' Cette fonction détecte ce cas pour déclencher un fallback.
#'
#' @param r SpatRaster à valider
#' @param min_pct Pourcentage minimum de pixels non-vides requis (défaut: 5%)
#' @return TRUE si le raster contient suffisamment de données
validate_wms_data <- function(r, min_pct = 5) {
  if (is.null(r)) return(FALSE)
  vals <- values(r[[1]])
  n_valid <- sum(!is.na(vals) & vals > 0)
  pct <- n_valid / length(vals) * 100
  if (pct < min_pct) {
    message(sprintf("  Données insuffisantes : %.1f%% de pixels valides (%d/%d)",
                    pct, n_valid, length(vals)))
  }
  return(pct >= min_pct)
}

#' Télécharger les ortho RVB et IRC pour une AOI
#'
#' Gestion du cache : si ortho_rvb.tif et ortho_irc.tif existent déjà,
#' ils sont réutilisés sans re-téléchargement.
#'
#' Gestion des millésimes (indépendants pour IRC et ortho) :
#'   - millesime_ortho : NULL = mosaïque nationale, entier = année spécifique
#'   - millesime_irc : NULL = mosaïque la plus récente, entier = année spécifique
#'   Si le millésime demandé n'est pas disponible (erreur WMS ou données vides),
#'   fallback automatique sur la mosaïque courante.
#'
#' @param aoi sf object (AOI en Lambert-93)
#' @param output_dir Répertoire de sortie
#' @param res_m Résolution en mètres
#' @param millesime_ortho NULL ou entier (année de l'ortho RVB)
#' @param millesime_irc NULL ou entier (année de l'ortho IRC)
#' @return Liste avec rvb, irc (SpatRaster), chemins et info millésime
download_ortho_for_aoi <- function(aoi, output_dir, res_m = RES_IGN,
                                    millesime_ortho = MILLESIME_ORTHO,
                                    millesime_irc = MILLESIME_IRC) {
  dir_create(output_dir)

  layer_ortho <- ign_layer_name("ortho", millesime_ortho)
  layer_irc   <- ign_layer_name("irc",   millesime_irc)
  label_ortho <- if (is.null(millesime_ortho)) "plus récent" else millesime_ortho
  label_irc   <- if (is.null(millesime_irc))   "plus récent" else millesime_irc

  rvb_path <- file.path(output_dir, "ortho_rvb.tif")
  irc_path <- file.path(output_dir, "ortho_irc.tif")

  # Cache : réutiliser les fichiers existants
  if (file.exists(rvb_path) && file.exists(irc_path)) {
    message("\n=== Ortho IGN déjà présentes (cache) ===")
    message(sprintf("  RVB: %s", rvb_path))
    message(sprintf("  IRC: %s", irc_path))
    message("Réutilisation des fichiers existants.")

    rvb <- rast(rvb_path)
    irc <- rast(irc_path)
    names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]
    names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

    return(list(rvb = rvb, irc = irc,
                rvb_path = rvb_path, irc_path = irc_path,
                millesime_ortho = millesime_ortho,
                millesime_irc = millesime_irc,
                layer_ortho = layer_ortho,
                layer_irc = layer_irc))
  }

  bbox <- as.numeric(st_bbox(st_union(aoi)))

  message(sprintf("\n=== Téléchargement ortho IGN pour l'AOI ==="))
  message(sprintf("Emprise: %.0f, %.0f - %.0f, %.0f (Lambert-93)",
                   bbox[1], bbox[2], bbox[3], bbox[4]))
  message(sprintf("Taille: %.0f x %.0f m (%.2f ha)",
                   bbox[3] - bbox[1], bbox[4] - bbox[2],
                   (bbox[3] - bbox[1]) * (bbox[4] - bbox[2]) / 10000))
  message(sprintf("Millésime RVB: %s (couche: %s)", label_ortho, layer_ortho))
  message(sprintf("Millésime IRC: %s (couche: %s)", label_irc, layer_irc))

  # --- RVB (avec fallback si millésime indisponible ou données vides) ---
  message("\n--- Ortho RVB ---")
  rvb <- tryCatch(
    download_ign_tiled(bbox, layer = layer_ortho, res_m = res_m,
                       output_dir = output_dir, prefix = "rvb"),
    error = function(e) {
      message("  Erreur téléchargement RVB: ", e$message)
      NULL
    }
  )

  # Fallback : si erreur ou données vides (millésime non couvert pour cette zone)
  if (!is.null(millesime_ortho) &&
      (is.null(rvb) || !validate_wms_data(rvb))) {
    message(sprintf("  Millésime %s indisponible pour cette zone, fallback sur %s",
                    millesime_ortho, IGN_LAYER_ORTHO))
    layer_ortho <- IGN_LAYER_ORTHO
    label_ortho <- "plus récent (fallback)"
    # Nettoyer les tuiles du premier essai
    tile_files <- dir_ls(output_dir, glob = "rvb_tile_*.tif")
    if (length(tile_files) > 0) file_delete(tile_files)
    rvb <- download_ign_tiled(bbox, layer = IGN_LAYER_ORTHO, res_m = res_m,
                               output_dir = output_dir, prefix = "rvb")
  }
  if (is.null(rvb)) stop("Impossible de télécharger l'ortho RVB")
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]

  # --- IRC (avec fallback si millésime indisponible ou données vides) ---
  message("\n--- Ortho IRC ---")
  irc <- tryCatch(
    download_ign_tiled(bbox, layer = layer_irc, res_m = res_m,
                       output_dir = output_dir, prefix = "irc"),
    error = function(e) {
      message("  Erreur téléchargement IRC: ", e$message)
      NULL
    }
  )

  # Fallback : si erreur ou données vides (millésime non couvert pour cette zone)
  if (!is.null(millesime_irc) &&
      (is.null(irc) || !validate_wms_data(irc))) {
    message(sprintf("  Millésime %s indisponible pour cette zone, fallback sur %s",
                    millesime_irc, IGN_LAYER_IRC))
    layer_irc <- IGN_LAYER_IRC
    label_irc <- "plus récent (fallback)"
    # Nettoyer les tuiles du premier essai
    tile_files <- dir_ls(output_dir, glob = "irc_tile_*.tif")
    if (length(tile_files) > 0) file_delete(tile_files)
    irc <- download_ign_tiled(bbox, layer = IGN_LAYER_IRC, res_m = res_m,
                               output_dir = output_dir, prefix = "irc")
  }
  if (is.null(irc)) stop("Impossible de télécharger l'ortho IRC")
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  # Découper aux limites exactes de l'AOI
  aoi_vect <- vect(st_union(aoi))
  rvb <- crop(rvb, aoi_vect)
  irc <- crop(irc, aoi_vect)

  # Sauvegarder les mosaïques finales
  writeRaster(rvb, rvb_path, overwrite = TRUE)
  writeRaster(irc, irc_path, overwrite = TRUE)

  # Re-lire depuis les fichiers sauvegardés (terra est file-backed)
  rvb <- rast(rvb_path)
  irc <- rast(irc_path)
  names(rvb)[1:min(3, nlyr(rvb))] <- c("Rouge", "Vert", "Bleu")[1:min(3, nlyr(rvb))]
  names(irc)[1:min(3, nlyr(irc))] <- c("PIR", "Rouge", "Vert")[1:min(3, nlyr(irc))]

  message(sprintf("\nRVB sauvegardé: %s (%d x %d px)", rvb_path, ncol(rvb), nrow(rvb)))
  message(sprintf("IRC sauvegardé: %s (%d x %d px)", irc_path, ncol(irc), nrow(irc)))

  # Nettoyer les tuiles temporaires
  tile_files <- dir_ls(output_dir, glob = "*_tile_*.tif")
  if (length(tile_files) > 0) file_delete(tile_files)

  return(list(rvb = rvb, irc = irc,
              rvb_path = rvb_path, irc_path = irc_path,
              millesime_ortho = millesime_ortho,
              millesime_irc = millesime_irc,
              layer_ortho = layer_ortho,
              layer_irc = layer_irc))
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

  # Calculer le CHM (hauteur au-dessus du sol = DSM - DTM)
  if (has_mns) {
    chm <- dsm - dtm
    chm[chm < 0] <- 0
    names(chm) <- "CHM"
    chm_vals <- values(chm, na.rm = TRUE)
    message(sprintf("  CHM (DSM-DTM): range [%.1f, %.1f]m",
                     min(chm_vals), max(chm_vals)))
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
  dtm_vals <- values(dem[["DTM"]], na.rm = TRUE)
  if (length(dtm_vals) > 0 && any(is.finite(dtm_vals))) {
    message(sprintf("  Altitude DTM: %.0f - %.0f m",
                     min(dtm_vals, na.rm = TRUE),
                     max(dtm_vals, na.rm = TRUE)))
    if (has_mns) {
      chm <- dem[["DSM"]] - dem[["DTM"]]
      chm_vals <- values(chm, na.rm = TRUE)
      if (length(chm_vals) > 0 && any(is.finite(chm_vals))) {
        message(sprintf("  Hauteur CHM (DSM-DTM): %.1f - %.1f m",
                         min(chm_vals, na.rm = TRUE),
                         max(chm_vals, na.rm = TRUE)))
      }
    }
  } else {
    message("  ATTENTION: DEM contient uniquement des NA")
    message("  Le WMS d'élévation a peut-être retourné des données vides.")
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

#' Trouver le nom du fichier checkpoint dans un dépôt HF via l'API
#'
#' @param hf_repo Identifiant du dépôt HF (ex: "IGNF/FLAIR-INC_rgbi_15cl_resnet34-unet")
#' @return Nom du fichier checkpoint (.ckpt, .pth, .pt, .bin, .safetensors)
find_checkpoint_name <- function(hf_repo) {
  url <- paste0("https://huggingface.co/api/models/", hf_repo)
  resp <- tryCatch(
    httr2::request(url) |> httr2::req_perform(),
    error = function(e) NULL
  )
  if (is.null(resp)) return(NULL)

  info <- jsonlite::fromJSON(httr2::resp_body_string(resp))
  files <- info$siblings$rfilename
  ckpt_files <- files[grepl("\\.(ckpt|pth|pt|bin|safetensors)$", files)]
  if (length(ckpt_files) == 0) return(NULL)
  # Prioriser : .ckpt > .pth > .pt > .bin > .safetensors
  ext_priority <- c("\\.ckpt$", "\\.pth$", "\\.pt$", "\\.bin$", "\\.safetensors$")
  for (pat in ext_priority) {
    matches <- ckpt_files[grepl(pat, ckpt_files)]
    if (length(matches) > 0) return(matches[1])
  }
  return(ckpt_files[1])
}

#' Télécharger un modèle FLAIR depuis Hugging Face
#'
#' Utilise le package R hfhub (natif, sans Python) pour télécharger
#' le fichier checkpoint. Fallback sur Python huggingface_hub si hfhub
#' n'est pas installé.
#'
#' @param model_name Nom du modèle (voir FLAIR_MODELS pour la liste)
#' @return Chemin local du modèle (fichier ou répertoire)
#' Résoudre les symlinks/junctions HuggingFace sur Windows
#'
#' Sur Windows sans Developer Mode, le cache HF utilise des junctions NTFS
#' dans snapshots/ qui provoquent NotADirectoryError en Python.
#' Cette fonction trouve le vrai fichier blob correspondant.
#'
#' @param path Chemin retourné par hfhub::hub_download()
#' @return Chemin résolu vers le blob réel
resolve_hf_symlink <- function(path) {
  if (!grepl("snapshots", path)) return(path)

  # Utiliser suppressWarnings pour éviter les avertissements sur les junctions NTFS
  # (file.info et readLines échouent bruyamment sur les junctions Windows)

  # 1. Sys.readlink()
  real <- tryCatch(Sys.readlink(path), error = function(e) "")
  if (nchar(real) > 0) {
    real_resolved <- suppressWarnings(normalizePath(real, mustWork = FALSE))
    if (file.exists(real_resolved)) {
      message("  Symlink résolu: ", basename(path), " -> ", basename(real_resolved))
      return(real_resolved)
    }
  }

  # 2. Pointer file (HF stocke le hash du blob dans un petit fichier texte)
  sz <- suppressWarnings(tryCatch(file.info(path)$size, error = function(e) NA))
  if (!is.na(sz) && sz < 1000) {
    content <- suppressWarnings(
      tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) ""))
    if (length(content) > 0 && nchar(content[1]) > 10) {
      blob_path <- file.path(dirname(dirname(dirname(path))),
                             "blobs", trimws(content[1]))
      if (file.exists(blob_path)) {
        message("  Pointer file résolu -> blobs/", substr(trimws(content[1]), 1, 12), "...")
        return(normalizePath(blob_path, mustWork = FALSE))
      }
    }
  }

  # 3. Plus gros blob (les poids sont le plus gros fichier)
  blobs_dir <- file.path(dirname(dirname(dirname(path))), "blobs")
  if (dir.exists(blobs_dir)) {
    blobs <- list.files(blobs_dir, full.names = TRUE)
    blobs <- blobs[!grepl("\\.(lock|incomplete)$", blobs)]
    if (length(blobs) > 0) {
      sizes <- file.info(blobs)$size
      biggest <- blobs[which.max(sizes)]
      if (max(sizes, na.rm = TRUE) > 1e6) {
        message(sprintf("  Blob résolu: %s (%.1f Mo)",
                        basename(biggest), max(sizes) / 1e6))
        return(normalizePath(biggest, mustWork = FALSE))
      }
    }
  }

  normalizePath(path, mustWork = FALSE)
}

download_model <- function(model_name = "FLAIR-INC_rgbie_15cl_resnet34-unet") {
  config <- FLAIR_MODELS[[model_name]]
  hf_repo <- if (!is.null(config)) config$hf_repo else paste0("IGNF/", model_name)

  message("Téléchargement du modèle: ", model_name)
  message("Depuis: ", hf_repo)

  # --- Méthode 1 : hfhub R natif (préféré) ---
  if (requireNamespace("hfhub", quietly = TRUE)) {
    ckpt_name <- find_checkpoint_name(hf_repo)
    if (!is.null(ckpt_name)) {
      message("  Téléchargement via hfhub (R natif): ", ckpt_name)
      tryCatch({
        local_path <- hfhub::hub_download(hf_repo, ckpt_name)
        # Vérifier que le fichier existe réellement sur le disque
        if (file.exists(local_path)) {
          # Résoudre les symlinks HF sur Windows (junctions → blobs/)
          resolved <- resolve_hf_symlink(local_path)
          message("  Modèle téléchargé: ", resolved)
          return(resolved)
        } else {
          message("  ATTENTION: hfhub a retourné un chemin mais le fichier n'existe pas:")
          message("    ", local_path)
          message("  Le dépôt HuggingFace est peut-être privé/gated ou le téléchargement a échoué.")
          message("  → fallback Python")
        }
      }, error = function(e) {
        message("  hfhub échoué: ", e$message, " → fallback Python")
      })
    } else {
      message("  Aucun fichier checkpoint trouvé via l'API HuggingFace pour: ", hf_repo)
      message("  Vérifiez que le dépôt existe et contient un fichier .ckpt/.pth/.pt/.bin/.safetensors")
    }
  }

  # --- Méthode 2 : Python huggingface_hub (fallback) ---
  library(reticulate)
  hf_hub <- import("huggingface_hub")
  tryCatch({
    local_dir <- hf_hub$snapshot_download(repo_id = hf_repo, repo_type = "model")
    # Vérifier que le répertoire contient bien des fichiers de poids
    if (dir.exists(local_dir)) {
      ckpt_files <- list.files(local_dir,
                               pattern = "\\.(ckpt|pth|pt|bin|safetensors)$",
                               recursive = TRUE)
      if (length(ckpt_files) == 0) {
        stop("Le dépôt '", hf_repo, "' a été téléchargé mais ne contient aucun ",
             "fichier de poids (.ckpt/.pth/.pt/.bin/.safetensors).\n",
             "Contenu du répertoire: ",
             paste(list.files(local_dir, recursive = TRUE), collapse = ", "),
             call. = FALSE)
      }
      message("  Modèle téléchargé: ", local_dir)
      message("  Fichiers de poids: ", paste(ckpt_files, collapse = ", "))
    } else if (file.exists(local_dir)) {
      message("  Modèle téléchargé (fichier): ", local_dir)
    } else {
      stop("Le téléchargement a retourné un chemin inexistant: ", local_dir, call. = FALSE)
    }
    return(local_dir)
  }, error = function(e) {
    stop("Échec du téléchargement du modèle '", model_name, "' depuis '", hf_repo, "'.\n",
         "Erreur: ", e$message, "\n",
         "Vérifiez:\n",
         "  1. Le dépôt existe sur HuggingFace\n",
         "  2. Vous êtes authentifié (huggingface-cli login) si le dépôt est privé/gated\n",
         "  3. Votre connexion internet fonctionne",
         call. = FALSE)
  })
}

#' Découper en patches pour l'inférence
#'
#' Garantit que tous les patches font exactement patch_size × patch_size pixels.
#' Le dernier patch de chaque axe est calé contre le bord du raster (avec un
#' overlap supplémentaire si nécessaire) plutôt que d'être tronqué, ce qui
#' évite les prédictions aberrantes causées par un padding excessif.
make_inference_patches <- function(r, patch_size = PATCH_SIZE, overlap = 32) {
  pixel_res <- res(r)[1]
  patch_size_m <- patch_size * pixel_res
  overlap_m <- overlap * pixel_res
  step_m <- patch_size_m - overlap_m

  e <- ext(r)
  raster_w <- e[2] - e[1]
  raster_h <- e[4] - e[3]

  # Grille régulière
  x_starts <- seq(e[1], e[2] - patch_size_m + step_m, by = step_m)
  y_starts <- seq(e[3], e[4] - patch_size_m + step_m, by = step_m)

  if (length(x_starts) == 0) x_starts <- e[1]
  if (length(y_starts) == 0) y_starts <- e[3]

  # Garantir un dernier patch plein (512×512) calé contre le bord du raster.
  # Sans cela, le dernier patch peut être très petit (ex: 50 px de large)
  # et le padding reflect produit des prédictions aberrantes.
  if (raster_w > patch_size_m) {
    x_last <- e[2] - patch_size_m
    if (abs(tail(x_starts, 1) - x_last) > pixel_res * 0.5) {
      x_starts <- c(x_starts, x_last)
    }
  }
  if (raster_h > patch_size_m) {
    y_last <- e[4] - patch_size_m
    if (abs(tail(y_starts, 1) - y_last) > pixel_res * 0.5) {
      y_starts <- c(y_starts, y_last)
    }
  }

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

#' Inférence sur un patch avec un modèle FLAIR (multi-architecture)
#'
#' Instancie le modèle smp selon la config (encoder/decoder dynamiques),
#' charge les poids, normalise l'image, exécute l'inférence PyTorch,
#' et remappe les classes FLAIR-1 vers la nomenclature CoSIA.
#'
#' @param patch SpatRaster (1 patch)
#' @param model_path Chemin vers le checkpoint (.ckpt) ou le répertoire du modèle
#' @param model_config Liste avec encoder, decoder, in_channels, n_out
#'   (NULL = config par défaut ResNet34 + UNet)
#' @return SpatRaster 1 bande (landcover) ou NULL en cas d'erreur
#' Charger le modèle FLAIR en mémoire Python (une seule fois)
#'
#' @param model_path Chemin vers le checkpoint
#' @param model_config Liste avec encoder, decoder, in_channels, n_out
#' @return TRUE si le modèle a été chargé, FALSE sinon
load_flair_model <- function(model_path, model_config) {
  library(reticulate)

  model_path_py <- gsub("\\\\", "/", model_path)

  py_code <- '
import os
import torch
import numpy as np
import segmentation_models_pytorch as smp

# ======================================================================
# Charger le modèle FLAIR (une seule fois)
# ======================================================================
model_dir = "__MODEL_PATH__"
ckpt_path = None
_flair_model = None
_flair_model_loaded = False

if os.path.isdir(model_dir):
    ckpt_exts = [".ckpt", ".pth", ".pt", ".bin", ".safetensors"]
    candidates = []
    for f in os.listdir(model_dir):
        for i, ext in enumerate(ckpt_exts):
            if f.endswith(ext):
                candidates.append((i, f))
                break
    if candidates:
        candidates.sort()
        ckpt_path = os.path.join(model_dir, candidates[0][1])
elif os.path.isfile(model_dir):
    ckpt_path = model_dir
else:
    raise FileNotFoundError(f"Aucun fichier de poids trouvé: {model_dir}")

print(f"Fichier modèle: {os.path.basename(ckpt_path)}")

in_ch = __IN_CHANNELS__
n_out = __N_OUT__

decoder_class = getattr(smp, "__DECODER__")
_flair_model = decoder_class(
    encoder_name="__ENCODER__",
    encoder_weights=None,
    in_channels=in_ch,
    classes=n_out,
)

# Charger le checkpoint
if ckpt_path.endswith(".safetensors"):
    from safetensors.torch import load_file
    state_dict = load_file(ckpt_path)
else:
    checkpoint = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    if isinstance(checkpoint, dict):
        if "state_dict" in checkpoint:
            state_dict = checkpoint["state_dict"]
        elif "model_state_dict" in checkpoint:
            state_dict = checkpoint["model_state_dict"]
        else:
            state_dict = checkpoint
    else:
        state_dict = checkpoint.state_dict() if hasattr(checkpoint, "state_dict") else {}

# Auto-détection du préfixe
model_keys = set(_flair_model.state_dict().keys())
ckpt_keys = list(state_dict.keys())

candidate_prefixes = ["", "model.", "net.", "module.", "backbone.",
                      "model.model.", "network.", "seg_model."]
for k in ckpt_keys[:20]:
    parts = k.split(".")
    for i in range(1, min(4, len(parts))):
        p = ".".join(parts[:i]) + "."
        if p not in candidate_prefixes:
            candidate_prefixes.append(p)

best_prefix = ""
best_match = 0
for try_prefix in candidate_prefixes:
    matches = sum(1 for k in ckpt_keys
                  if k.startswith(try_prefix) and
                  k[len(try_prefix):] in model_keys)
    if matches > best_match:
        best_match = matches
        best_prefix = try_prefix

print(f"  Préfixe: \'{best_prefix}\' ({best_match}/{len(model_keys)} clés)")

cleaned = {}
for k, v in state_dict.items():
    if best_prefix and k.startswith(best_prefix):
        cleaned[k[len(best_prefix):]] = v
    elif not best_prefix:
        cleaned[k] = v

missing, unexpected = _flair_model.load_state_dict(cleaned, strict=False)
n_missing = len(missing) if missing else 0
if n_missing > len(model_keys) * 0.5:
    raise RuntimeError(f"Trop de clés manquantes ({n_missing}/{len(model_keys)})")

_flair_model.eval()
_flair_model_loaded = True
print(f"Modèle chargé: __DECODER__(__ENCODER__), {in_ch} bandes, {n_out} classes")
'
  py_code <- gsub("__MODEL_PATH__", model_path_py, py_code, fixed = TRUE)
  py_code <- gsub("__DECODER__", model_config$decoder, py_code, fixed = TRUE)
  py_code <- gsub("__ENCODER__", model_config$encoder, py_code, fixed = TRUE)
  py_code <- gsub("__IN_CHANNELS__", as.character(model_config$in_channels), py_code, fixed = TRUE)
  py_code <- gsub("__N_OUT__", as.character(model_config$n_out), py_code, fixed = TRUE)

  tryCatch({
    py_run_string(py_code)
    return(TRUE)
  }, error = function(e) {
    message("ERREUR chargement modèle: ", e$message)
    return(FALSE)
  })
}

predict_patch <- function(patch, model_path, model_config = NULL, ...) {
  if (is.null(model_config)) {
    model_config <- list(encoder = "resnet34", decoder = "Unet",
                         in_channels = 4, n_out = 19)
  }
  library(reticulate)

  tmp_in <- tempfile(fileext = ".tif")
  tmp_out <- tempfile(fileext = ".tif")
  writeRaster(patch, tmp_in, overwrite = TRUE)

  tmp_in_py <- gsub("\\\\", "/", tmp_in)
  tmp_out_py <- gsub("\\\\", "/", tmp_out)

  py_code <- '
import numpy as np
import torch
import rasterio

# Charger l image
with rasterio.open("__INPUT_PATH__") as src:
    image = src.read().astype(np.float32)
    profile = src.profile.copy()

num_bands, H, W = image.shape

# Inférence avec le modèle pré-chargé
if _flair_model_loaded:
    in_ch = __IN_CHANNELS__
    # Normalisation FLAIR (centre-réduit, statistiques TRAIN+VAL)
    norm_means = np.array([105.08, 110.87, 101.82, 106.38, 53.26], dtype=np.float32)
    norm_stds  = np.array([52.17, 45.38, 44.0, 39.69, 79.3], dtype=np.float32)

    img = image[:in_ch]
    for c in range(img.shape[0]):
        if c < len(norm_means):
            img[c] = (img[c] - norm_means[c]) / norm_stds[c]

    pad_h = max(0, 512 - H)
    pad_w = max(0, 512 - W)
    if pad_h > 0 or pad_w > 0:
        img = np.pad(img, ((0, 0), (0, pad_h), (0, pad_w)), mode="reflect")

    tensor = torch.from_numpy(img).unsqueeze(0)

    with torch.no_grad():
        logits = _flair_model(tensor)

    pred_flair = logits.squeeze(0).cpu().numpy().argmax(axis=0)

    if pad_h > 0 or pad_w > 0:
        pred_flair = pred_flair[:H, :W]

    # Remap FLAIR-1 → CoSIA
    remap = np.array([
        1, 5, 4, 6, 7, 14, 13, 15, 12, 9, 10, 11, 3, 8,
        0, 0, 0, 2, 0], dtype=np.int32)
    pred = remap[pred_flair]
else:
    # Fallback NDVI + CHM (classification par seuillage spectral et hauteur)
    if num_bands >= 4:
        r, g, b, nir = image[0], image[1], image[2], image[3]
        ndvi = (nir - r) / (nir + r + 1e-6)
        brightness = (r + g + b) / 3.0
        has_chm = num_bands >= 5
        chm = image[4] if has_chm else np.zeros((H, W), dtype=np.float32)

        pred = np.zeros((H, W), dtype=np.int32)

        # --- Classes avec CHM (hauteur au-dessus du sol) ---
        if has_chm:
            # Eau : sombre, pas de végétation, au sol
            pred[(brightness < 30) & (ndvi < 0.1) & (chm < 1)] = 7
            # Piscine : eau claire (bleu)
            pred[(b > r) & (b > nir) & (ndvi < 0) & (chm < 1) & (pred == 0)] = 3
            # Bâtiment : élevé, pas de végétation
            pred[(chm > 3) & (ndvi < 0.15) & (pred == 0)] = 1
            # Serre : élevé modéré, très brillant, pas de végétation
            pred[(chm > 2) & (chm < 6) & (brightness > 140) & (ndvi < 0.2) & (pred == 0)] = 2
            # Conifère : très haut, très vert
            pred[(chm > 10) & (ndvi >= 0.5) & (pred == 0)] = 14
            # Feuillu : haut, vert
            pred[(chm > 5) & (ndvi >= 0.3) & (pred == 0)] = 13
            # Lande/broussaille : hauteur moyenne, végétation
            pred[(chm > 1) & (chm <= 5) & (ndvi >= 0.2) & (pred == 0)] = 15
            # Imperméable : au sol, clair, pas de végétation
            pred[(chm < 1) & (brightness > 100) & (ndvi < 0.15) & (pred == 0)] = 4
            # Sol nu : au sol, peu de végétation
            pred[(chm < 1) & (ndvi < 0.2) & (pred == 0)] = 6
            # Herbacé : au sol, végétation faible à moyenne
            pred[(chm < 2) & (ndvi >= 0.2) & (ndvi < 0.4) & (pred == 0)] = 9
            # Agricole : au sol, végétation moyenne
            pred[(chm < 2) & (ndvi >= 0.4) & (pred == 0)] = 10
            # Reste non classé
            pred[pred == 0] = 9  # par défaut herbacé
        else:
            # --- Fallback NDVI seul (sans CHM) ---
            pred[(brightness < 30) & (ndvi < 0.1)] = 7
            pred[(brightness > 150) & (ndvi < 0.1) & (pred == 0)] = 1
            pred[(brightness > 100) & (ndvi < 0.15) & (pred == 0)] = 4
            pred[(ndvi < 0.2) & (pred == 0)] = 6
            pred[(ndvi >= 0.2) & (ndvi < 0.35) & (pred == 0)] = 9
            pred[(ndvi >= 0.35) & (ndvi < 0.5) & (pred == 0)] = 10
            pred[(ndvi >= 0.5) & (ndvi < 0.7) & (pred == 0)] = 13
            pred[(ndvi >= 0.7) & (pred == 0)] = 14
            pred[pred == 0] = 15

        n_unique = len(np.unique(pred[pred > 0]))
        mode = "NDVI+CHM" if has_chm else "NDVI seul"
        print(f"Fallback ({mode}): {n_unique} classes")
    else:
        pred = np.full((H, W), 15, dtype=np.int32)

# Construire un profil propre (ne pas hériter les paramètres du fichier source
# qui sont incompatibles avec une sortie 1 bande int32)
out_profile = {
    "driver": "GTiff",
    "dtype": "int32",
    "width": W,
    "height": H,
    "count": 1,
    "crs": profile.get("crs"),
    "transform": profile.get("transform"),
    "compress": "lzw",
}
with rasterio.open("__OUTPUT_PATH__", "w", **out_profile) as dst:
    dst.write(pred.astype(np.int32), 1)

print(f"Prédit: {np.unique(pred).shape[0]} classes uniques")
'
  py_code <- gsub("__INPUT_PATH__", tmp_in_py, py_code, fixed = TRUE)
  py_code <- gsub("__OUTPUT_PATH__", tmp_out_py, py_code, fixed = TRUE)
  py_code <- gsub("__IN_CHANNELS__", as.character(model_config$in_channels), py_code, fixed = TRUE)

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

#' Pipeline d'inférence complet avec blending par fenêtre de Hann
#'
#' Utilise une fenêtre de Hann (cosinus) 2D pour pondérer chaque patch
#' par la distance à son centre. Pour chaque pixel, la prédiction retenue
#' est celle du patch dont le centre est le plus proche (poids maximal).
#' Cela élimine naturellement les artefacts de bord du réseau de neurones
#' et produit une mosaïque sans coutures visibles.
#'
#' @param rgbi SpatRaster 4 bandes (Rouge, Vert, Bleu, PIR)
#' @param model_path Chemin vers le modèle FLAIR-HUB
#' @param buffer_px Contrôle l'overlap entre patches : overlap = max(64, 2 × buffer_px).
#'   Plus la valeur est élevée, plus la zone de recouvrement est large et le
#'   blending efficace (mais l'inférence est plus lente car plus de patches).
#'   Valeur recommandée : 64 (overlap = 128 px).
#' @param model_config Liste avec encoder, decoder, in_channels, n_out
#'   (NULL = config par défaut ResNet34 + UNet)
#' @return SpatRaster 1 bande (landcover)
run_inference <- function(rgbi, model_path, buffer_px = 0, model_config = NULL,
                          fallback = FALSE) {
  message("\n=== Inférence FLAIR-HUB ===")

  # Charger le modèle une seule fois en mémoire Python
  if (is.null(model_config)) {
    model_config <- list(encoder = "resnet34", decoder = "Unet",
                         in_channels = 4, n_out = 19, n_classes = 15)
  }
  if (fallback) {
    message("  Mode FALLBACK forcé : classification NDVI",
            if (nlyr(rgbi) >= 5) "+CHM" else "", " (pas de réseau de neurones)")
    # Forcer _flair_model_loaded = False en Python
    library(reticulate)
    py_run_string("_flair_model_loaded = False")
  } else {
    model_ok <- load_flair_model(model_path, model_config)
    if (!model_ok) {
      stop("Impossible de charger le modèle. Vérifiez le fichier de poids.", call. = FALSE)
    }
  }

  # Overlap : au moins 64 px, ou 2× buffer_px pour un bon recouvrement
  effective_overlap <- max(64, if (buffer_px > 0) 2 * buffer_px else 64)
  message(sprintf("  Overlap: %d px, vote pondéré par classe (Hann)", effective_overlap))

  patches <- make_inference_patches(rgbi, overlap = effective_overlap)
  n_total <- length(patches)

  # Cas trivial : un seul patch, pas de blending nécessaire
  if (n_total == 1) {
    message(sprintf("  Patch 1/1: %s", names(patches)[1]))
    result <- predict_patch(patches[[1]], model_path, model_config = model_config)
    if (is.null(result)) stop("Échec de l'inférence sur le seul patch.")
    names(result) <- "landcover"
    return(result)
  }

  # --- Vote pondéré par classe avec fenêtre de Hann ---
  # Principe : pour chaque pixel couvert par plusieurs patches, on
  # accumule le poids Hann pour chaque classe prédite. La classe avec
  # le poids cumulé le plus élevé l'emporte (soft majority voting).
  # Cela élimine les artefacts de bord car les prédictions du centre
  # (poids élevé) dominent naturellement celles des bords (poids faible).

  out_ext  <- ext(rgbi)
  out_res  <- res(rgbi)[1]
  n_rows   <- nrow(rgbi)
  n_cols   <- ncol(rgbi)

  # Nombre de classes CoSIA (1..15)
  n_classes <- if (!is.null(model_config)) {
    max(model_config$n_classes, 15)
  } else 15

  # Matrice de votes : accumule le poids Hann par classe pour chaque pixel
  # Dimension : n_classes × (n_rows * n_cols) — stockage en colonnes pour R
  vote_mat <- matrix(0, nrow = n_classes, ncol = n_rows * n_cols)

  n_ok <- 0
  for (i in seq_along(patches)) {
    patch_name <- names(patches)[i]
    message(sprintf("  Patch %d/%d: %s", i, n_total, patch_name))

    pred <- predict_patch(patches[[i]], model_path, model_config = model_config)
    if (is.null(pred)) next
    n_ok <- n_ok + 1

    nr <- nrow(pred)
    nc <- ncol(pred)

    # Fenêtre de Hann 2D : poids(r,c) = wy[r] * wx[c]
    wy <- 0.5 - 0.5 * cos(2 * pi * seq(0.5, nr - 0.5) / nr)
    wx <- 0.5 - 0.5 * cos(2 * pi * seq(0.5, nc - 0.5) / nc)
    w_mat <- outer(wy, wx)

    # Prédictions sous forme de matrice (lignes = rangées de pixels)
    pred_vec <- values(pred)[, 1]
    pred_mat <- matrix(pred_vec, nrow = nr, ncol = nc, byrow = TRUE)

    # Offsets (0-based) dans le raster de sortie
    pe <- ext(pred)
    col_off <- round((pe[1] - out_ext[1]) / out_res)
    row_off <- round((out_ext[4] - pe[4]) / out_res)

    # Indices 1-based dans les matrices de sortie
    r1 <- row_off + 1L;  r2 <- row_off + nr
    c1 <- col_off + 1L;  c2 <- col_off + nc

    # Sécurité : borner aux dimensions de sortie
    pr1 <- 1L; pr2 <- nr; pc1 <- 1L; pc2 <- nc
    if (r1 < 1L)      { pr1 <- 2L - r1;  r1 <- 1L }
    if (r2 > n_rows)   { pr2 <- nr - (r2 - n_rows);  r2 <- n_rows }
    if (c1 < 1L)      { pc1 <- 2L - c1;  c1 <- 1L }
    if (c2 > n_cols)   { pc2 <- nc - (c2 - n_cols);  c2 <- n_cols }

    # Sous-matrices du patch alignées sur la zone de sortie
    w_sub    <- w_mat[pr1:pr2, pc1:pc2]
    pred_sub <- pred_mat[pr1:pr2, pc1:pc2]

    # Accumuler les votes pondérés par classe (vectorisé)
    # Indices linéaires row-major (= ordre des cellules terra) :
    # cell_index = (row - 1) * n_cols + col
    out_rows <- r1:r2
    out_cols <- c1:c2
    idx_mat <- outer((out_rows - 1L) * n_cols, out_cols, "+")  # (rows, cols)

    pred_flat <- as.integer(pred_sub)
    w_flat <- as.numeric(w_sub)
    idx_flat <- as.integer(idx_mat)

    # Ajouter les votes par classe (15 itérations max, vectorisé par classe)
    for (k in seq_len(n_classes)) {
      mask <- pred_flat == k
      if (any(mask)) {
        idx_k <- idx_flat[mask]
        # Accumuler les poids avec tapply pour gérer les doublons éventuels
        vote_mat[k, idx_k] <- vote_mat[k, idx_k] + w_flat[mask]
      }
    }
  }

  if (n_ok == 0) stop("Aucune prédiction réussie.")

  # Classe finale = classe avec le vote pondéré le plus élevé
  class_vec <- apply(vote_mat, 2, which.max)
  # Pixels sans vote → 0
  no_vote <- colSums(vote_mat) == 0
  class_vec[no_vote] <- 0L

  message(sprintf("  Vote pondéré terminé (%d patches, %d classes détectées)",
                   n_ok, length(unique(class_vec[class_vec > 0]))))

  # Convertir en SpatRaster (column-major = ordre naturel de R)
  result <- rast(rgbi[[1]])
  values(result) <- as.integer(class_vec)
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
#' @param dem_res_m Résolution du MNT (1 = RGE ALTI 1m)
#' @param millesime_ortho NULL (mosaïque la plus récente) ou entier (ex: 2024)
#' @param millesime_irc NULL (mosaïque la plus récente) ou entier (ex: 2024)
#' @param buffer_px Contrôle l'overlap entre patches pour le blending
#'   (défaut: BUFFER_PX). overlap = max(64, 2 × buffer_px). Le mosaïquage
#'   utilise une fenêtre de Hann 2D : chaque pixel reçoit la prédiction du
#'   patch dont le centre est le plus proche, éliminant les artefacts de bord.
#'   Valeur recommandée : 64 (overlap = 128 px).
#' @return Liste avec tous les résultats
pipeline_aoi_to_landcover <- function(aoi_path,
                                        output_dir = file.path(getwd(), "outputs"),
                                        model_name = "FLAIR-INC_rgbie_15cl_resnet34-unet",
                                        model_path = NULL,
                                        res_m = RES_IGN,
                                        use_dem = TRUE,
                                        dem_res_m = 1,
                                        millesime_ortho = MILLESIME_ORTHO,
                                        millesime_irc = MILLESIME_IRC,
                                        buffer_px = BUFFER_PX,
                                        fallback = FALSE) {
  dir_create(output_dir)
  t0 <- Sys.time()

  config_label <- if (use_dem) "LC-B (RGBI + MNT)" else "LC-A (RGBI seul)"
  n_steps <- if (use_dem) 6 else 5

  message("##############################################################")
  message("#  Pipeline FLAIR-HUB : AOI → Ortho IGN → Occupation du sol  #")
  message(sprintf("#  Configuration: %s", config_label))
  if (buffer_px > 0) {
    message(sprintf("#  Buffer inférence: %d px par bord", buffer_px))
  }
  message("##############################################################\n")

  # --- Étape 1 : Charger l'AOI ---
  message(sprintf(">>> ÉTAPE 1/%d : Chargement de l'AOI", n_steps))
  aoi <- load_aoi(aoi_path)

  # --- Étape 2 : Télécharger les ortho IGN ---
  message(sprintf("\n>>> ÉTAPE 2/%d : Téléchargement des ortho IGN (RVB + IRC)",
                   n_steps))
  ortho <- download_ortho_for_aoi(aoi, output_dir = output_dir, res_m = res_m,
                                    millesime_ortho = millesime_ortho,
                                    millesime_irc = millesime_irc)

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

  # Récupérer la config architecture du modèle
  model_config <- FLAIR_MODELS[[model_name]]

  # --- Combiner RGBI + Élévation si le modèle le requiert (RGBIE, 5 canaux) ---
  inference_input <- rgbi  # par défaut : 4 bandes RGBI

  if (!is.null(model_config) && model_config$in_channels >= 5) {
    if (!is.null(dem_data)) {
      # 5ème canal : CHM (hauteur au-dessus du sol = DSM - DTM)
      # Plus informatif que le DTM brut pour la classification :
      # - Arbres : CHM 5-30m, Bâtiments : CHM 3-15m, Sol : CHM ~0m
      # - Pas d'artefact de dalle car le DTM (RGE ALTI) est continu
      dsm <- dem_data$dem[["DSM"]]
      dtm <- dem_data$dem[["DTM"]]
      chm <- dsm - dtm
      chm[chm < 0] <- 0   # Borner les valeurs négatives
      chm[chm > 60] <- 60  # Borner les valeurs aberrantes
      names(chm) <- "CHM"

      # Nettoyer les éventuels NA
      if (any(is.na(values(chm)))) {
        chm <- focal(chm, w = 5, fun = "mean", na.policy = "only", na.rm = TRUE)
        chm[is.na(chm)] <- 0
      }

      chm_range <- range(values(chm), na.rm = TRUE)
      message(sprintf("  CHM (DSM-DTM): range [%.1f, %.1f]m", chm_range[1], chm_range[2]))

      inference_input <- c(rgbi, chm)
      names(inference_input) <- c("Rouge", "Vert", "Bleu", "PIR", "CHM")

      rgbie_path <- file.path(output_dir, "ortho_rgbie.tif")
      writeRaster(inference_input, rgbie_path, overwrite = TRUE)

      message(sprintf("  Image RGBIE: %d x %d px, %d bandes (RGBI + Élévation DSM)",
                       ncol(inference_input), nrow(inference_input),
                       nlyr(inference_input)))
    } else {
      warning("Le modèle ", model_name, " attend ", model_config$in_channels,
              " canaux (RGBIE) mais le MNT n'est pas disponible.\n",
              "  Inférence avec 4 canaux (RGBI) seulement.\n",
              "  Activez use_dem = TRUE pour de meilleurs résultats.")
    }
  }

  # --- Étape 4/5 : Inférence ---
  step_inf <- if (use_dem) 5 else 4
  message(sprintf("\n>>> ÉTAPE %d/%d : Inférence du modèle %s",
                   step_inf, n_steps, model_name))

  landcover <- run_inference(inference_input, model_path, buffer_px = buffer_px,
                             model_config = model_config, fallback = fallback)

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
  pdf_path <- file.path(output_dir, paste0("resultats_", model_name, ".pdf"))
  has_dem <- !is.null(dem_data)
  model_uses_elev <- !is.null(model_config) && model_config$in_channels >= 5
  show_dem <- has_dem && (use_dem || model_uses_elev)
  n_panels <- if (show_dem) 6 else 4
  pdf_w <- if (n_panels > 4) 18 else 16
  pdf(pdf_path, width = pdf_w, height = 12)

  if (n_panels > 4) {
    par(mfrow = c(2, 3), mar = c(2, 2, 3, 4), oma = c(0, 0, 3, 0))
  } else {
    par(mfrow = c(2, 2), mar = c(2, 2, 3, 4), oma = c(0, 0, 3, 0))
  }

  # RVB
  tryCatch(
    plotRGB(ortho$rvb, r = 1, g = 2, b = 3, stretch = "lin",
            main = "Ortho RVB IGN (0.20m)"),
    error = function(e) {
      plot.new(); title(main = paste("RVB - erreur:", e$message))
    }
  )

  # IRC fausses couleurs
  # L'IRC IGN a 3 bandes : PIR, Rouge, Vert.
  # plotRGB(r=1,g=2,b=3) → fausses couleurs (PIR en rouge, végétation en rouge vif).
  # Si le WMS a retourné un canal alpha (4 bandes), on garde seulement les 3 premières.
  irc_for_plot <- ortho$irc
  if (nlyr(irc_for_plot) > 3) {
    irc_for_plot <- irc_for_plot[[1:3]]
  }
  tryCatch(
    plotRGB(irc_for_plot, r = 1, g = 2, b = 3, stretch = "lin",
            main = "Ortho IRC fausses couleurs (0.20m)"),
    error = function(e) {
      # Fallback : afficher la bande PIR seule si plotRGB échoue
      tryCatch({
        col_pir <- colorRampPalette(c("black", "red", "yellow", "white"))(100)
        plot(ortho$irc[[1]], main = "PIR (bande 1 IRC)",
             col = col_pir, plg = list(title = "PIR"))
      }, error = function(e2) {
        plot.new()
        title(main = paste("IRC - erreur:", e$message))
      })
    }
  )

  # NDVI
  col_ndvi <- colorRampPalette(
    c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
      "#d9ef8b", "#91cf60", "#1a9850", "#006837")
  )(100)
  plot(ndvi, main = "NDVI (depuis IRC)", col = col_ndvi,
       range = c(-0.2, 1), plg = list(title = "NDVI"))

  # MNT / DSM si disponible
  if (show_dem) {
    col_elev <- colorRampPalette(
      c("#313695", "#4575b4", "#74add1", "#abd9e9", "#fee090",
        "#fdae61", "#f46d43", "#d73027", "#a50026")
    )(100)

    if (model_uses_elev) {
      # Modèle RGBIE : montrer le DSM (5ème bande d'entrée) et le MNT
      plot(dem_data$dem[["DSM"]],
           main = "DSM - Élévation (5ème bande modèle)",
           col = col_elev, plg = list(title = "Altitude (m)"))
      plot(dem_data$dem[["DTM"]],
           main = sprintf("MNT IGN (RGE ALTI %dm)", dem_res_m),
           col = col_elev, plg = list(title = "Altitude (m)"))
    } else {
      # Modèle sans élévation : MNT + CHM
      plot(dem_data$dem[["DTM"]],
           main = sprintf("MNT IGN (RGE ALTI %dm)", dem_res_m),
           col = col_elev, plg = list(title = "Altitude (m)"))

      chm <- dem_data$dem[["DSM"]] - dem_data$dem[["DTM"]]
      col_chm <- colorRampPalette(
        c("#ffffcc", "#d9f0a3", "#addd8e", "#78c679",
          "#41ab5d", "#238443", "#005a32")
      )(100)
      plot(chm, main = "CHM (DSM - DTM)",
           col = col_chm, plg = list(title = "Hauteur (m)"))
    }
  }

  # Occupation du sol — raster catégoriel (uniquement les classes présentes)
  lc_plot <- landcover
  lc_present <- sort(unique(na.omit(as.integer(values(landcover)))))
  lc_ids <- lc_present
  lc_labels <- vapply(lc_present, function(cls) {
    if (cls == 0) "Non classifié" else COSIA_LABELS_15[cls]
  }, character(1))
  lc_colors <- vapply(lc_present, function(cls) {
    if (cls == 0) "#808080" else COSIA_COLORS_15[cls]
  }, character(1))
  levels(lc_plot) <- data.frame(id = lc_ids, label = lc_labels)
  plot(lc_plot, main = paste0("Occupation du sol - ", model_name, "\n", config_label),
       col = lc_colors, type = "classes",
       plg = list(legend = lc_labels, cex = 0.6, border = NA))

  # Titre global du PDF avec le nom du modèle
  mtext(paste0("FLAIR-HUB : ", model_name, "  —  ", config_label),
        outer = TRUE, cex = 1.4, font = 2, line = 1)

  dev.off()
  message("PDF:               ", pdf_path)

  # --- Statistiques ---
  lc_vals <- values(landcover, na.rm = TRUE)
  class_counts <- table(as.integer(lc_vals))

  message("\n--- Distribution des classes ---")
  for (i in seq_along(class_counts)) {
    cls <- as.integer(names(class_counts)[i])
    pct <- as.numeric(class_counts[i]) / length(lc_vals) * 100
    label <- if (cls == 0) "Non classifié" else if (cls >= 1 && cls <= 15) COSIA_LABELS_15[cls] else "?"
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
    inference_input = inference_input,
    ndvi            = ndvi,
    landcover       = landcover,
    output_dir      = output_dir,
    model_name      = model_name,
    config_label    = config_label
  )
  if (!is.null(dem_data)) result$dem <- dem_data$dem

  # --- Export patchwork en PDF (ggplot2) ---
  tryCatch({
    p <- plot_results(result, maxcell = 100000)
    if (!is.null(p)) {
      gg_pdf <- file.path(output_dir, paste0("resultats_", model_name, "_ggplot.pdf"))
      gg_w <- if (!is.null(dem_data)) 18 else 14
      ggplot2::ggsave(gg_pdf, plot = p, width = gg_w, height = 10, device = "pdf")
      message("PDF (ggplot):      ", gg_pdf)
    }
  }, error = function(e) {
    message("Export patchwork ignoré: ", e$message)
  })

  return(result)
}

# ==============================================================================
# Visualisation interactive (ggplot2 + tidyterra + patchwork)
# ==============================================================================

#' Afficher les résultats du pipeline dans RStudio
#'
#' Crée un assemblage patchwork identique au PDF exporté.
#' Nécessite ggplot2, tidyterra et patchwork.
#'
#' @param result Liste retournée par \code{pipeline_aoi_to_landcover()}
#' @return Un objet patchwork (affiché automatiquement dans RStudio)
#' @export
plot_results <- function(result, maxcell = 100000) {

  # --- Vérification des packages ---
  pkgs <- c("ggplot2", "tidyterra", "patchwork")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Packages manquants pour l'affichage RStudio: ",
            paste(missing, collapse = ", "))
    message('  install.packages(c("ggplot2", "tidyterra", "patchwork"))')
    return(invisible(NULL))
  }

  library(ggplot2)
  library(tidyterra)
  library(patchwork)

  # Palette NDVI
  col_ndvi <- c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
                "#d9ef8b", "#91cf60", "#1a9850", "#006837")

  # Limiter le nombre de cellules pour éviter le timeout RStudio
  # 100K = rapide, 250K = moyen, 500K+ = lent/timeout possible
  mc <- maxcell

  # --- Panel 1 : Ortho RVB ---
  p_rvb <- ggplot() +
    geom_spatraster_rgb(data = result$ortho_rvb, r = 1, g = 2, b = 3,
                        max_col_value = 255, maxcell = mc) +
    ggtitle("Ortho RVB IGN (0.20m)") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))

  # --- Panel 2 : IRC fausses couleurs ---
  # Si le WMS a retourné un canal alpha (4 bandes), garder les 3 premières
  irc_data <- result$ortho_irc
  if (nlyr(irc_data) > 3) irc_data <- irc_data[[1:3]]

  p_irc <- tryCatch({
    ggplot() +
      geom_spatraster_rgb(data = irc_data, r = 1, g = 2, b = 3,
                          max_col_value = 255, maxcell = mc) +
      ggtitle("Ortho IRC fausses couleurs (0.20m)") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
  }, error = function(e) {
    # Fallback : afficher la bande PIR seule
    ggplot() +
      geom_spatraster(data = result$ortho_irc[[1]], maxcell = mc) +
      scale_fill_gradient(low = "black", high = "red", name = "PIR",
                          na.value = "transparent") +
      ggtitle("PIR (bande 1 IRC)") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
  })

  # --- Panel 3 : NDVI ---
  p_ndvi <- ggplot() +
    geom_spatraster(data = result$ndvi, maxcell = mc) +
    scale_fill_gradientn(colours = col_ndvi, na.value = "transparent",
                         limits = c(-0.2, 1), name = "NDVI") +
    ggtitle("NDVI (depuis IRC)") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "right")

  # --- Panel 4 (optionnel) : MNT ---
  p_dtm <- NULL
  p_chm <- NULL
  if (!is.null(result$dem)) {
    col_elev <- c("#313695", "#4575b4", "#74add1", "#abd9e9", "#fee090",
                  "#fdae61", "#f46d43", "#d73027", "#a50026")

    dtm <- result$dem[["DTM"]]
    p_dtm <- ggplot() +
      geom_spatraster(data = dtm, maxcell = mc) +
      scale_fill_gradientn(colours = col_elev, na.value = "transparent",
                           name = "Altitude (m)") +
      ggtitle("MNT IGN (RGE ALTI)") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
            legend.position = "right")

    chm <- result$dem[["DSM"]] - result$dem[["DTM"]]
    col_chm <- c("#ffffcc", "#d9f0a3", "#addd8e", "#78c679",
                 "#41ab5d", "#238443", "#005a32")
    p_chm <- ggplot() +
      geom_spatraster(data = chm, maxcell = mc) +
      scale_fill_gradientn(colours = col_chm, na.value = "transparent",
                           name = "Hauteur (m)") +
      ggtitle("CHM (DSM - DTM)") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
            legend.position = "right")
  }

  # --- Panel Occupation du sol ---
  lc <- result$landcover

  # Identifier les classes présentes (uniquement celles dans les données)
  lc_present <- sort(unique(na.omit(as.integer(values(lc)))))
  lc_present <- lc_present[lc_present >= 1 & lc_present <= 15]

  # Reclasser en facteur avec uniquement les classes présentes
  lc_factor <- as.factor(lc)
  levels(lc_factor) <- data.frame(
    id    = lc_present,
    label = COSIA_LABELS_15[lc_present]
  )

  cls_colors <- setNames(COSIA_COLORS_15[lc_present],
                         COSIA_LABELS_15[lc_present])

  p_lc <- ggplot() +
    geom_spatraster(data = lc_factor, maxcell = mc) +
    scale_fill_manual(values = cls_colors,
                      na.value = "transparent", name = "Classe",
                      drop = TRUE) +
    guides(fill = guide_legend(override.aes = list(colour = NA))) +
    ggtitle(paste0("Occupation du sol - ", result$model_name %||% "FLAIR-HUB")) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          legend.position = "right",
          legend.text = element_text(size = 7),
          legend.key = element_rect(colour = NA))

  # --- Assemblage patchwork ---
  ann_title <- paste0("FLAIR-HUB : ", result$model_name %||% "Résultats du pipeline")
  ann_subtitle <- paste0(
    "Occupation du sol par segmentation sémantique (IGN)",
    if (!is.null(result$config_label)) paste0(" - ", result$config_label) else ""
  )

  if (!is.null(p_dtm)) {
    # Layout 2x3 (avec DEM)
    combined <- (p_rvb | p_irc | p_ndvi) /
                (p_dtm | p_chm | p_lc) +
      plot_annotation(
        title    = ann_title,
        subtitle = ann_subtitle,
        theme    = theme(
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40")
        )
      )
  } else {
    # Layout 2x2 (sans DEM)
    combined <- (p_rvb | p_irc) /
                (p_ndvi | p_lc) +
      plot_annotation(
        title    = ann_title,
        subtitle = ann_subtitle,
        theme    = theme(
          plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40")
        )
      )
  }

  return(combined)
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
    message('  # Modèle par défaut (RGBI+Élévation, 15 classes, ResNet34+UNet)')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg")')
    message("")
    message('  # Modèle RGBI seul (sans MNT)')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    model_name = "FLAIR-INC_rgbi_15cl_resnet34-unet", use_dem = FALSE)')
    message("")
    message('  # Modèle RGB 15 classes (sans IRC ni MNT)')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    model_name = "FLAIR-INC_rgb_15cl_resnet34-unet", use_dem = FALSE)')
    message("")
    message('  # Avec un modèle local :')
    message('  result <- pipeline_aoi_to_landcover("data/aoi.gpkg",')
    message('    model_path = "chemin/vers/modele/")')
    message("")
    message("Modèles disponibles :")
    for (nm in names(FLAIR_MODELS)) {
      cfg <- FLAIR_MODELS[[nm]]
      message(sprintf("  - %s  [%s(%s), %dch]",
                        nm, cfg$decoder, cfg$encoder, cfg$in_channels))
    }
    message("")
    message("Le fichier aoi.gpkg doit contenir un polygone définissant")
    message("votre zone d'intérêt (n'importe quel CRS, sera reprojeté")
    message("en Lambert-93 automatiquement).")
  } else {
    result <- pipeline_aoi_to_landcover(aoi_path)
  }
}
