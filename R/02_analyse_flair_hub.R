#!/usr/bin/env Rscript
# ==============================================================================
# 02_analyse_flair_hub.R
# Analyse et visualisation des données FLAIR-HUB :
#   - Images aériennes RGBI (0.2m)
#   - Images SPOT RGBI (1.6m)
#   - Sentinel-2 / Sentinel-1 séries temporelles
#   - MNT (DSM, DTM)
#   - Photos aériennes historiques 1950s
#   - Labels CoSIA (19 classes occupation du sol)
#   - Labels LPIS (23 classes cultures)
# ==============================================================================

library(terra)
library(sf)
library(fs)

# ==============================================================================
# Configuration
# ==============================================================================

DATA_DIR    <- file.path(getwd(), "data")
DATA_DIR_HF <- file.path(DATA_DIR, "flair_hub")
OUTPUT_DIR  <- file.path(getwd(), "outputs")
dir_create(OUTPUT_DIR)

# Palettes de couleurs pour les classes CoSIA (19 classes)
COSIA_COLORS <- c(
  "#db0e9a",  # 1  Bâtiment

"#938e7b",  # 2  Serre
  "#f80c00",  # 3  Piscine
  "#a97101",  # 4  Surface imperméable
  "#1553ae",  # 5  Surface perméable
  "#194a26",  # 6  Sol nu
  "#46e483",  # 7  Eau
  "#f3a60d",  # 8  Neige
  "#660082",  # 9  Végétation herbacée
  "#55ff00",  # 10 Terre agricole
  "#fff30d",  # 11 Terre labourée
  "#e4df7c",  # 12 Vigne
  "#3de6eb",  # 13 Verger
  "#ffffff",  # 14 Feuillu
  "#8ab3a0",  # 15 Conifère
  "#6b714f",  # 16 Lande
  "#c5dc42",  # 17 Ligneux mélangé
  "#9999ff",  # 18 Fleur / Garrigue
  "#000000"   # 19 Non classé
)

COSIA_LABELS <- c(
  "Bâtiment", "Serre", "Piscine", "Imperméable", "Perméable",
  "Sol nu", "Eau", "Neige", "Herbacé", "Agricole",
  "Labouré", "Vigne", "Verger", "Feuillu", "Conifère",
  "Lande", "Ligneux mélangé", "Fleur/Garrigue", "Non classé"
)

# Palettes LPIS (23 classes)
LPIS_COLORS <- c(
  "#808080",  # 0  Fond
  "#FFD700",  # 1  Blé tendre
  "#FF8C00",  # 2  Maïs
  "#DAA520",  # 3  Orge
  "#FFFF00",  # 4  Colza
  "#FFA500",  # 5  Tournesol
  "#90EE90",  # 6  Prairie temp.
  "#228B22",  # 7  Prairie perm.
  "#D2B48C",  # 8  Jachère
  "#8FBC8F",  # 9  Soja
  "#006400",  # 10 Légumineuse fourr.
  "#FF69B4",  # 11 Betterave
  "#CD853F",  # 12 Pomme de terre
  "#8B4513",  # 13 Autre perm.
  "#BDB76B",  # 14 Culture mélangée
  "#F0E68C",  # 15 Autre oléagineux
  "#32CD32",  # 16 Protéagineux
  "#00CED1",  # 17 Riz
  "#FF1493",  # 18 Légume/Fleur
  "#FF4500",  # 19 Sorgho
  "#B8860B",  # 20 Blé dur
  "#DEB887",  # 21 Avoine
  "#D2691E"   # 22 Triticale
)

LPIS_LABELS <- c(
  "Fond", "Blé tendre", "Maïs", "Orge", "Colza", "Tournesol",
  "Prairie temp.", "Prairie perm.", "Jachère", "Soja",
  "Légumineuse fourr.", "Betterave", "Pomme de terre",
  "Autre perm.", "Culture mélangée", "Autre oléagineux",
  "Protéagineux", "Riz", "Légume/Fleur", "Sorgho",
  "Blé dur", "Avoine", "Triticale"
)

# ==============================================================================
# 1. Indices spectraux (depuis aérien RGBI ou SPOT RGBI)
# ==============================================================================

#' Calculer le NDVI depuis une image RGBI
#'
#' NDVI = (PIR - Rouge) / (PIR + Rouge)
#'
#' @param rgbi_raster SpatRaster avec bandes Rouge et PIR
#' @return SpatRaster du NDVI (valeurs entre -1 et 1)
compute_ndvi <- function(rgbi_raster) {
  pir <- rgbi_raster[["PIR"]]
  rouge <- rgbi_raster[["Rouge"]]

  ndvi <- (pir - rouge) / (pir + rouge)
  names(ndvi) <- "NDVI"

  vals <- values(ndvi, na.rm = TRUE)
  message(sprintf("NDVI calculé: min=%.3f, max=%.3f, moy=%.3f",
                   min(vals), max(vals), mean(vals)))
  return(ndvi)
}

#' Calculer le GNDVI (Green NDVI)
#'
#' GNDVI = (PIR - Vert) / (PIR + Vert)
#'
#' @param rgbi_raster SpatRaster avec bandes Vert et PIR
#' @return SpatRaster du GNDVI
compute_gndvi <- function(rgbi_raster) {
  pir <- rgbi_raster[["PIR"]]
  vert <- rgbi_raster[["Vert"]]

  gndvi <- (pir - vert) / (pir + vert)
  names(gndvi) <- "GNDVI"
  return(gndvi)
}

#' Calculer le SAVI (Soil-Adjusted Vegetation Index)
#'
#' SAVI = ((PIR - R) / (PIR + R + L)) * (1 + L)
#'
#' @param rgbi_raster SpatRaster avec bandes Rouge et PIR
#' @param L Facteur d'ajustement sol (0.5 par défaut)
#' @return SpatRaster du SAVI
compute_savi <- function(rgbi_raster, L = 0.5) {
  pir <- rgbi_raster[["PIR"]]
  rouge <- rgbi_raster[["Rouge"]]

  savi <- ((pir - rouge) / (pir + rouge + L)) * (1 + L)
  names(savi) <- "SAVI"
  return(savi)
}

#' Créer un masque de végétation à partir du NDVI
#'
#' @param ndvi_raster SpatRaster du NDVI
#' @param threshold Seuil NDVI pour considérer de la végétation
#' @return SpatRaster binaire (1 = végétation)
mask_vegetation <- function(ndvi_raster, threshold = 0.3) {
  veg_mask <- ndvi_raster >= threshold
  names(veg_mask) <- "vegetation"

  pct <- sum(values(veg_mask, na.rm = TRUE)) / sum(!is.na(values(veg_mask))) * 100
  message(sprintf("Végétation détectée (NDVI >= %.2f): %.1f%%", threshold, pct))
  return(veg_mask)
}

# ==============================================================================
# 2. Visualisation
# ==============================================================================

#' Visualiser une image en couleurs naturelles (RGB)
#'
#' @param raster_rgb SpatRaster avec au moins 3 bandes
#' @param title Titre du graphique
#' @param bands Indices des bandes RGB
plot_rgb <- function(raster_rgb, title = "Image RGB", bands = c(1, 2, 3)) {
  if (nlyr(raster_rgb) >= 3) {
    plotRGB(raster_rgb, r = bands[1], g = bands[2], b = bands[3],
            stretch = "lin", main = title)
  } else {
    plot(raster_rgb, main = title)
  }
}

#' Visualiser les labels CoSIA (19 classes d'occupation du sol)
#'
#' @param label_raster SpatRaster des labels CoSIA
#' @param title Titre
plot_label_cosia <- function(label_raster, title = "Occupation du sol (CoSIA)") {
  vals <- values(label_raster, na.rm = TRUE)
  present_classes <- sort(unique(as.integer(vals)))

  # Filtrer les couleurs et labels pour les classes présentes
  colors <- COSIA_COLORS[present_classes]
  labels <- COSIA_LABELS[present_classes]

  plot(label_raster, main = title, col = COSIA_COLORS,
       type = "classes", levels = COSIA_LABELS,
       plg = list(legend = labels, cex = 0.7))
}

#' Visualiser les labels LPIS (23 classes de cultures)
#'
#' @param label_raster SpatRaster des labels LPIS
#' @param title Titre
plot_label_lpis <- function(label_raster, title = "Types de cultures (LPIS)") {
  vals <- values(label_raster, na.rm = TRUE)
  present_classes <- sort(unique(as.integer(vals)))

  # Classes LPIS indexées à partir de 0
  colors <- LPIS_COLORS[present_classes + 1]
  labels <- LPIS_LABELS[present_classes + 1]

  plot(label_raster, main = title, col = LPIS_COLORS,
       type = "classes", levels = LPIS_LABELS,
       plg = list(legend = labels, cex = 0.6))
}

#' Visualiser le NDVI
#'
#' @param ndvi_raster SpatRaster du NDVI
#' @param title Titre
plot_ndvi <- function(ndvi_raster, title = "NDVI") {
  col_ndvi <- colorRampPalette(
    c("#d73027", "#fc8d59", "#fee08b", "#ffffbf",
      "#d9ef8b", "#91cf60", "#1a9850", "#006837")
  )(100)

  plot(ndvi_raster, main = title, col = col_ndvi, range = c(-0.2, 1),
       plg = list(title = "NDVI"))
}

#' Visualiser le MNT (DSM/DTM)
#'
#' @param dem_raster SpatRaster avec bandes DSM et DTM
#' @param title Titre
plot_dem <- function(dem_raster, title = "MNT") {
  col_elev <- colorRampPalette(
    c("#313695", "#4575b4", "#74add1", "#abd9e9",
      "#e0f3f8", "#ffffbf", "#fee090", "#fdae61",
      "#f46d43", "#d73027", "#a50026")
  )(100)

  if (nlyr(dem_raster) >= 2) {
    par(mfrow = c(1, 2), mar = c(2, 2, 3, 4))
    plot(dem_raster[[1]], main = paste(title, "- DSM"),
         col = col_elev, plg = list(title = "Altitude (m)"))
    plot(dem_raster[[2]], main = paste(title, "- DTM"),
         col = col_elev, plg = list(title = "Altitude (m)"))
    par(mfrow = c(1, 1))
  } else {
    plot(dem_raster, main = title, col = col_elev,
         plg = list(title = "Altitude (m)"))
  }
}

#' Comparaison multi-modalités FLAIR-HUB
#'
#' @param aerial SpatRaster aérien RGBI (0.2m)
#' @param label_cosia SpatRaster labels CoSIA
#' @param spot SpatRaster SPOT RGBI (1.6m, optionnel)
#' @param dem SpatRaster MNT (optionnel)
#' @param tile_name Nom de la tuile
plot_multimodal_comparison <- function(aerial, label_cosia,
                                        spot = NULL, dem = NULL,
                                        tile_name = "") {
  n_panels <- 2 + !is.null(spot) + !is.null(dem)
  ncols <- min(n_panels, 3)
  nrows <- ceiling(n_panels / ncols)

  par(mfrow = c(nrows, ncols), mar = c(2, 2, 3, 4))

  # Aérien RGBI
  plot_rgb(aerial, title = paste("Aérien RGBI 0.2m", tile_name))

  # Labels CoSIA
  plot_label_cosia(label_cosia, title = paste("CoSIA", tile_name))

  # SPOT si disponible
  if (!is.null(spot)) {
    plot_rgb(spot, title = paste("SPOT RGBI 1.6m", tile_name))
  }

  # MNT si disponible
  if (!is.null(dem) && nlyr(dem) >= 2) {
    chm <- dem[["DSM"]] - dem[["DTM"]]
    names(chm) <- "CHM"
    col_chm <- colorRampPalette(
      c("#f7fcb9", "#addd8e", "#41ab5d", "#006837", "#004529")
    )(100)
    plot(chm, main = paste("CHM (DSM-DTM)", tile_name),
         col = col_chm, plg = list(title = "Hauteur (m)"))
  }

  par(mfrow = c(1, 1))
}

#' Visualiser une série temporelle Sentinel-2
#'
#' @param s2_raster SpatRaster multi-bandes (séries temporelles)
#' @param n_dates Nombre de dates à afficher
#' @param title Titre
plot_sentinel2_timeseries <- function(s2_raster, n_dates = 6,
                                       title = "Sentinel-2 SITS") {
  n_bands <- nlyr(s2_raster)
  n_show <- min(n_dates, n_bands)

  # Sélectionner des dates espacées uniformément
  indices <- round(seq(1, n_bands, length.out = n_show))

  ncols <- min(n_show, 3)
  nrows <- ceiling(n_show / ncols)

  par(mfrow = c(nrows, ncols), mar = c(2, 2, 3, 2))
  for (i in indices) {
    plot(s2_raster[[i]], main = paste(title, "- Bande", i),
         col = terrain.colors(100))
  }
  par(mfrow = c(1, 1))
}

# ==============================================================================
# 3. Statistiques
# ==============================================================================

#' Statistiques d'occupation du sol (CoSIA)
#'
#' @param label_raster SpatRaster des labels CoSIA
#' @return data.frame avec distribution des classes
compute_landcover_stats <- function(label_raster) {
  vals <- values(label_raster, na.rm = TRUE)
  total <- length(vals)

  class_counts <- table(as.integer(vals))
  class_ids <- as.integer(names(class_counts))

  stats <- data.frame(
    class_id = class_ids,
    label = COSIA_LABELS[class_ids],
    n_pixels = as.integer(class_counts),
    pct = round(as.numeric(class_counts) / total * 100, 2),
    stringsAsFactors = FALSE
  )

  stats <- stats[order(-stats$pct), ]

  message("=== Distribution occupation du sol (CoSIA) ===")
  message(sprintf("  Total: %d pixels", total))
  for (i in seq_len(min(nrow(stats), 10))) {
    message(sprintf("  %2d. %s: %.1f%%",
                     stats$class_id[i], stats$label[i], stats$pct[i]))
  }

  return(stats)
}

#' Statistiques des types de cultures (LPIS)
#'
#' @param label_raster SpatRaster des labels LPIS
#' @return data.frame avec distribution des classes
compute_crop_stats <- function(label_raster) {
  vals <- values(label_raster, na.rm = TRUE)
  total <- length(vals)

  class_counts <- table(as.integer(vals))
  class_ids <- as.integer(names(class_counts))

  stats <- data.frame(
    class_id = class_ids,
    label = LPIS_LABELS[class_ids + 1],
    n_pixels = as.integer(class_counts),
    pct = round(as.numeric(class_counts) / total * 100, 2),
    stringsAsFactors = FALSE
  )

  stats <- stats[order(-stats$pct), ]

  message("=== Distribution des cultures (LPIS) ===")
  message(sprintf("  Total: %d pixels", total))
  for (i in seq_len(min(nrow(stats), 10))) {
    message(sprintf("  %2d. %s: %.1f%%",
                     stats$class_id[i], stats$label[i], stats$pct[i]))
  }

  return(stats)
}

#' Croiser occupation du sol et MNT
#'
#' @param label_raster SpatRaster des labels CoSIA
#' @param dem_raster SpatRaster MNT (DSM, DTM)
#' @return data.frame avec hauteur moyenne par classe
cross_landcover_dem <- function(label_raster, dem_raster) {
  # Calculer le CHM (hauteur de canopée)
  if (nlyr(dem_raster) >= 2) {
    chm <- dem_raster[["DSM"]] - dem_raster[["DTM"]]
  } else {
    chm <- dem_raster
  }

  # Aligner si nécessaire
  if (!compareGeom(label_raster, chm, stopOnError = FALSE)) {
    message("Rééchantillonnage du MNT vers la résolution des labels...")
    chm <- resample(chm, label_raster, method = "bilinear")
  }

  # Statistiques par classe
  label_vals <- values(label_raster, na.rm = TRUE)
  chm_vals <- values(chm, na.rm = TRUE)

  classes <- sort(unique(as.integer(label_vals)))
  result <- data.frame(
    class_id = integer(0),
    label = character(0),
    n_pixels = integer(0),
    height_mean = numeric(0),
    height_median = numeric(0),
    height_sd = numeric(0),
    height_max = numeric(0),
    stringsAsFactors = FALSE
  )

  for (cls in classes) {
    mask <- as.integer(label_vals) == cls
    h <- chm_vals[mask]
    h <- h[!is.na(h)]

    if (length(h) > 0) {
      result <- rbind(result, data.frame(
        class_id = cls,
        label = ifelse(cls >= 1 && cls <= 19, COSIA_LABELS[cls], "Unknown"),
        n_pixels = length(h),
        height_mean = round(mean(h), 2),
        height_median = round(median(h), 2),
        height_sd = round(sd(h), 2),
        height_max = round(max(h), 2),
        stringsAsFactors = FALSE
      ))
    }
  }

  message("=== Croisement occupation du sol x hauteur (CHM) ===")
  print(result)
  return(result)
}

# ==============================================================================
# 4. Export
# ==============================================================================

export_raster <- function(raster_obj, filename, output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)
  writeRaster(raster_obj, out_path, overwrite = TRUE)
  message("Raster exporté: ", out_path)
  return(out_path)
}

export_stats <- function(stats_df, filename = "statistics.csv",
                          output_dir = OUTPUT_DIR) {
  out_path <- file.path(output_dir, filename)
  write.csv(stats_df, out_path, row.names = FALSE)
  message("Statistiques exportées: ", out_path)
  return(out_path)
}

# ==============================================================================
# Exécution principale
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== FLAIR-HUB : Analyse et visualisation ===\n")

  # Rechercher les fichiers disponibles
  all_tifs <- unlist(lapply("*.tif", function(ext) {
    dir_ls(DATA_DIR, recurse = TRUE, glob = ext)
  }))

  if (length(all_tifs) == 0) {
    message("Aucun fichier trouvé dans ", DATA_DIR)
    message("Exécutez d'abord: Rscript R/01_download_flair_hub.R")
    message("\nDémonstration avec des données simulées...\n")

    # --- Simulation ---
    set.seed(42)

    # Simuler un patch aérien RGBI 512x512 à 0.2m
    demo_aerial <- rast(nrows = 512, ncols = 512, nlyrs = 4,
                         xmin = 843000, xmax = 843102.4,
                         ymin = 6518000, ymax = 6518102.4,
                         crs = "EPSG:2154")
    names(demo_aerial) <- c("Rouge", "Vert", "Bleu", "PIR")
    values(demo_aerial) <- cbind(
      Rouge = pmax(0, pmin(255, rnorm(ncell(demo_aerial), 100, 30))),
      Vert = pmax(0, pmin(255, rnorm(ncell(demo_aerial), 110, 30))),
      Bleu = pmax(0, pmin(255, rnorm(ncell(demo_aerial), 90, 25))),
      PIR = pmax(0, pmin(255, rnorm(ncell(demo_aerial), 160, 40)))
    )

    # Simuler les labels CoSIA
    demo_cosia <- rast(nrows = 512, ncols = 512,
                        xmin = 843000, xmax = 843102.4,
                        ymin = 6518000, ymax = 6518102.4,
                        crs = "EPSG:2154")
    values(demo_cosia) <- sample(c(1, 4, 6, 9, 10, 14, 15),
                                   ncell(demo_cosia), replace = TRUE,
                                   prob = c(0.1, 0.1, 0.05, 0.2, 0.2, 0.2, 0.15))
    names(demo_cosia) <- "landcover"

    # NDVI
    ndvi <- compute_ndvi(demo_aerial)

    # Statistiques
    message("\n--- Statistiques CoSIA ---")
    lc_stats <- compute_landcover_stats(demo_cosia)

    # Visualisation
    pdf(file.path(OUTPUT_DIR, "demo_flair_hub_analysis.pdf"),
        width = 14, height = 10)

    par(mfrow = c(2, 2), mar = c(2, 2, 3, 4))
    plot_rgb(demo_aerial, title = "Aérien RGBI 0.2m (simulé)")
    plot_label_cosia(demo_cosia, title = "CoSIA (simulé)")
    plot_ndvi(ndvi, title = "NDVI (simulé)")

    # Masque végétation
    veg <- mask_vegetation(ndvi)
    plot(veg, main = "Végétation (NDVI >= 0.3)",
         col = c("white", "#1a9850"))

    dev.off()
    message("\nGraphiques: ", file.path(OUTPUT_DIR, "demo_flair_hub_analysis.pdf"))

  } else {
    message(sprintf("%d fichier(s) trouvé(s)\n", length(all_tifs)))
  }

  message("\n=== Analyse terminée ===")
}
