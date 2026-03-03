#!/usr/bin/env Rscript
# ==============================================================================
# tester_flair_hub.R
# Script de test du pipeline FLAIR-HUB multi-architecture sur une AOI
#
# Ce script exécute le pipeline complet :
#   AOI (GeoPackage) → Ortho IGN (RVB + IRC) → RGBI → Inférence FLAIR → Landcover
#
# Il permet de tester un ou plusieurs modèles FLAIR et de comparer
# les résultats visuellement et quantitativement.
#
# Usage :
#   # Depuis la racine du package :
#   Rscript inst/script/tester_flair_hub.R data/aoi.gpkg
#
#   # Avec un modèle spécifique :
#   Rscript inst/script/tester_flair_hub.R data/aoi.gpkg FLAIR-INC_rgb_15cl_resnet34-unet
#
#   # Tester tous les modèles :
#   Rscript inst/script/tester_flair_hub.R data/aoi.gpkg ALL
#
#   # Depuis n'importe où (package installé) :
#   Rscript -e "source(system.file('script', 'tester_flair_hub.R', package='flairhub'))"
#
# Prérequis :
#   - R packages : terra, sf, fs, curl, httr2, jsonlite, reticulate
#   - R packages optionnels : hfhub (téléchargement natif), ggplot2, tidyterra, patchwork
#   - Python (conda env FLAIRHUB) : torch, numpy, rasterio, segmentation_models_pytorch, timm
# ==============================================================================

# ==============================================================================
# 0. Charger le pipeline
# ==============================================================================

# Trouver et sourcer le pipeline selon le contexte d'exécution
pipeline_loaded <- FALSE

# Méthode 1 : chemin relatif depuis la racine du package (développement)
pipeline_dev <- file.path("R", "04_pipeline_aoi_to_landcover.R")
if (file.exists(pipeline_dev)) {
  source(pipeline_dev)
  pipeline_loaded <- TRUE
  message("Pipeline chargé depuis: ", pipeline_dev)
}

# Méthode 2 : depuis inst/script/ (développement, exécution directe)
if (!pipeline_loaded) {
  pipeline_inst <- file.path(dirname(sys.frame(1)$ofile %||% "."),
                              "..", "..", "R", "04_pipeline_aoi_to_landcover.R")
  pipeline_inst <- normalizePath(pipeline_inst, mustWork = FALSE)
  if (file.exists(pipeline_inst)) {
    source(pipeline_inst)
    pipeline_loaded <- TRUE
    message("Pipeline chargé depuis: ", pipeline_inst)
  }
}

# Méthode 3 : package installé
if (!pipeline_loaded) {
  pkg_pipeline <- system.file("R", "04_pipeline_aoi_to_landcover.R",
                               package = "flairhub")
  if (nzchar(pkg_pipeline)) {
    source(pkg_pipeline)
    pipeline_loaded <- TRUE
    message("Pipeline chargé depuis le package installé")
  }
}

if (!pipeline_loaded) {
  stop("Impossible de trouver le pipeline 04_pipeline_aoi_to_landcover.R\n",
       "  Lancez ce script depuis la racine du package :\n",
       "    cd flair_hub_nemeton\n",
       "    Rscript inst/script/tester_flair_hub.R data/aoi.gpkg")
}

# ==============================================================================
# 1. Configuration (arguments CLI)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# --- AOI ---
aoi_path <- if (length(args) >= 1) args[1] else file.path("data", "aoi.gpkg")

# --- Modèle(s) ---
if (length(args) >= 2) {
  model_arg <- args[2]
  if (toupper(model_arg) == "ALL") {
    models_to_test <- names(FLAIR_MODELS)
  } else {
    models_to_test <- model_arg
  }
} else {
  # Par défaut : modèle RGBI 15 classes
  models_to_test <- "FLAIR-INC_rgbie_15cl_resnet34-unet"
}

# --- Répertoire de sortie ---
output_root <- if (length(args) >= 3) args[3] else file.path(getwd(), "outputs")

# --- Buffer (optionnel, 4ème argument) ---
buffer_px <- if (length(args) >= 4) as.integer(args[4]) else BUFFER_PX

# ==============================================================================
# 2. Vérifications
# ==============================================================================

if (!file.exists(aoi_path)) {
  stop("Fichier AOI introuvable : ", aoi_path, "\n\n",
       "Usage :\n",
       "  Rscript inst/script/tester_flair_hub.R <aoi.gpkg> [modèle|ALL] [output_dir] [buffer_px]\n\n",
       "Exemples :\n",
       "  Rscript inst/script/tester_flair_hub.R data/aoi.gpkg\n",
       "  Rscript inst/script/tester_flair_hub.R data/aoi.gpkg FLAIR-INC_rgb_15cl_resnet34-unet\n",
       "  Rscript inst/script/tester_flair_hub.R data/aoi.gpkg ALL outputs/ 64\n")
}

cat("\n")
message("################################################################")
message("#  FLAIR-HUB : Test du pipeline multi-architecture             #")
message("################################################################")
message("")
message("AOI         : ", aoi_path)
message("Modèle(s)   : ", paste(models_to_test, collapse = ", "))
message("Sortie      : ", output_root)
message("Buffer (px) : ", buffer_px)
message("")

# --- Afficher les modèles disponibles ---
message("--- Modèles FLAIR disponibles ---")
for (nm in names(FLAIR_MODELS)) {
  cfg <- FLAIR_MODELS[[nm]]
  marker <- if (nm %in% models_to_test) " <<" else ""
  message(sprintf("  %-45s  %s(%s)  %dch -> %d cls%s",
                   nm, cfg$decoder, cfg$encoder,
                   cfg$in_channels, cfg$n_classes, marker))
}
message("")

# ==============================================================================
# 3. Exécution du pipeline pour chaque modèle
# ==============================================================================

results  <- list()
timings  <- list()
statuses <- list()

for (model_name in models_to_test) {

  # Vérifier que le modèle existe dans le registre
  if (!model_name %in% names(FLAIR_MODELS)) {
    message(sprintf("\n[SKIP] Modèle inconnu : '%s'", model_name))
    message("  Modèles disponibles : ", paste(names(FLAIR_MODELS), collapse = ", "))
    statuses[[model_name]] <- "SKIP"
    next
  }

  config <- FLAIR_MODELS[[model_name]]

  # Sous-répertoire de sortie par modèle (ou directement output_root si un seul)
  if (length(models_to_test) == 1) {
    output_dir <- output_root
  } else {
    output_dir <- file.path(output_root, model_name)
  }

  message(sprintf("\n################################################################"))
  message(sprintf("#  Modèle : %s", model_name))
  message(sprintf("#  Architecture : %s(%s), %d canaux -> %d classes",
                   config$decoder, config$encoder,
                   config$in_channels, config$n_classes))
  message(sprintf("################################################################\n"))

  t0 <- Sys.time()

  result <- tryCatch({
    pipeline_aoi_to_landcover(
      aoi_path   = aoi_path,
      output_dir = output_dir,
      model_name = model_name,
      buffer_px  = buffer_px
    )
  }, error = function(e) {
    message(sprintf("\n[ERREUR] %s : %s", model_name, e$message))
    NULL
  })

  dt <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  timings[[model_name]] <- dt

  if (!is.null(result)) {
    results[[model_name]] <- result
    statuses[[model_name]] <- "OK"
    message(sprintf("\n[OK] %s termine en %s min", model_name, dt))
  } else {
    statuses[[model_name]] <- "ECHEC"
    message(sprintf("\n[ECHEC] %s apres %s min", model_name, dt))
  }
}

# ==============================================================================
# 4. Rapport comparatif
# ==============================================================================

if (length(results) > 0) {
  message("\n################################################################")
  message("#  Rapport comparatif                                          #")
  message("################################################################\n")

  for (model_name in names(results)) {
    result <- results[[model_name]]
    config <- FLAIR_MODELS[[model_name]]
    dt <- timings[[model_name]]

    lc_vals <- values(result$landcover, na.rm = TRUE)
    class_counts <- table(as.integer(lc_vals))
    n_classes <- length(class_counts[names(class_counts) != "0"])

    message(sprintf("--- %s (%s min) ---", model_name, dt))
    message(sprintf("  Architecture : %s(%s)", config$decoder, config$encoder))
    message(sprintf("  Classes predites : %d / %d", n_classes, config$n_classes))

    for (i in seq_along(class_counts)) {
      cls <- as.integer(names(class_counts)[i])
      pct <- as.numeric(class_counts[i]) / length(lc_vals) * 100
      if (cls >= 1 && cls <= 15) {
        label <- COSIA_LABELS_15[cls]
      } else if (cls == 0) {
        label <- "Non classifie"
      } else {
        label <- "?"
      }
      message(sprintf("    %2d. %-15s %5.1f%%", cls, label, pct))
    }
    message("")
  }
}

# ==============================================================================
# 5. Comparaison PDF (si plusieurs modèles)
# ==============================================================================

if (length(results) > 1) {
  message("--- Generation du PDF comparatif ---")
  pdf_path <- file.path(output_root, "comparaison_modeles.pdf")
  n_models <- length(results)

  pdf(pdf_path, width = min(8 * n_models, 32), height = 10)
  par(mfrow = c(1, n_models), mar = c(2, 2, 4, 4))

  for (model_name in names(results)) {
    result <- results[[model_name]]
    config <- FLAIR_MODELS[[model_name]]
    lc <- result$landcover

    lc_present <- sort(unique(na.omit(as.integer(values(lc)))))
    lc_labels <- vapply(lc_present, function(cls) {
      if (cls == 0) "Non classifie" else COSIA_LABELS_15[cls]
    }, character(1))
    lc_colors <- vapply(lc_present, function(cls) {
      if (cls == 0) "#808080" else COSIA_COLORS_15[cls]
    }, character(1))

    lc_plot <- lc
    levels(lc_plot) <- data.frame(id = lc_present, label = lc_labels)
    plot(lc_plot,
         main = sprintf("%s\n%s(%s)", model_name, config$decoder, config$encoder),
         col = lc_colors, type = "classes",
         plg = list(legend = lc_labels, cex = 0.5, border = NA))
  }

  dev.off()
  message("PDF comparatif : ", pdf_path)
}

# ==============================================================================
# 6. Resume final
# ==============================================================================

n_ok   <- sum(statuses == "OK")
n_fail <- sum(statuses == "ECHEC")
n_skip <- sum(statuses == "SKIP")

message("\n################################################################")
message(sprintf("#  Tests termines : %d OK, %d echec(s), %d ignore(s)",
                 n_ok, n_fail, n_skip))
if (length(results) > 0) {
  message(sprintf("#  Resultats dans : %s", output_root))
  message("#")
  message("#  Fichiers generes par modele :")
  message("#    - ortho_rvb.tif, ortho_irc.tif, ortho_rgbi.tif")
  message("#    - landcover_predicted.tif")
  message("#    - ndvi.tif")
  message("#    - resultats_aoi_flair_hub.pdf")
}
message("################################################################")
