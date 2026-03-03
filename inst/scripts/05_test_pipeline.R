#!/usr/bin/env Rscript
# ==============================================================================
# 05_test_pipeline.R
# Script de test complet du pipeline multi-architecture sur aoi.gpkg
#
# Ce script teste les différents modèles FLAIR disponibles sur une AOI fournie.
# Chaque modèle produit sa propre carte d'occupation du sol, ce qui permet
# de comparer les résultats visuellement et quantitativement.
#
# Usage :
#   Rscript R/05_test_pipeline.R
#   Rscript R/05_test_pipeline.R chemin/vers/mon_aoi.gpkg
#   Rscript R/05_test_pipeline.R data/aoi.gpkg FLAIR-INC_rgb_15cl_resnet34-unet
# ==============================================================================

source("R/04_pipeline_aoi_to_landcover.R")

# ==============================================================================
# Configuration
# ==============================================================================

# Chemin vers l'AOI (argument CLI ou défaut)
args <- commandArgs(trailingOnly = TRUE)
aoi_path <- if (length(args) >= 1) args[1] else file.path("data", "aoi.gpkg")

# Modèle(s) à tester (argument CLI ou tous les modèles disponibles)
if (length(args) >= 2) {
  models_to_test <- args[2]
} else {
  models_to_test <- names(FLAIR_MODELS)
}

# Répertoire de sortie racine
output_root <- file.path(getwd(), "outputs_test")

# ==============================================================================
# Vérifications
# ==============================================================================

if (!file.exists(aoi_path)) {
  stop("Fichier AOI introuvable : ", aoi_path, "\n",
       "  Placez votre fichier aoi.gpkg dans data/ ou spécifiez le chemin en argument.\n",
       "  Usage : Rscript R/05_test_pipeline.R data/aoi.gpkg")
}

message("================================================================")
message("  Test pipeline FLAIR-HUB multi-architecture")
message("================================================================")
message("AOI         : ", aoi_path)
message("Modèles     : ", paste(models_to_test, collapse = ", "))
message("Sortie      : ", output_root)
message("Buffer (px) : ", BUFFER_PX)
message("================================================================\n")

# ==============================================================================
# Afficher les modèles disponibles
# ==============================================================================

message("--- Modèles FLAIR disponibles ---")
for (nm in names(FLAIR_MODELS)) {
  cfg <- FLAIR_MODELS[[nm]]
  message(sprintf("  %-45s  %s(%s)  %dch → %d classes",
                   nm, cfg$decoder, cfg$encoder,
                   cfg$in_channels, cfg$n_classes))
}
message("")

# ==============================================================================
# Lancer le pipeline pour chaque modèle
# ==============================================================================

results <- list()
timings <- list()

for (model_name in models_to_test) {
  # Vérifier que le modèle existe dans le registre
  if (!model_name %in% names(FLAIR_MODELS)) {
    message(sprintf("\n[SKIP] Modèle inconnu : %s", model_name))
    message("  Modèles disponibles : ", paste(names(FLAIR_MODELS), collapse = ", "))
    next
  }

  config <- FLAIR_MODELS[[model_name]]
  output_dir <- file.path(output_root, model_name)

  message(sprintf("\n================================================================"))
  message(sprintf("  Modèle : %s", model_name))
  message(sprintf("  Architecture : %s(%s), %d canaux → %d classes",
                   config$decoder, config$encoder,
                   config$in_channels, config$n_classes))
  message(sprintf("  Élévation (DEM) : %s",
                   if (config$in_channels >= 5) "OUI (5ème bande)" else "non"))
  message(sprintf("================================================================\n"))

  t0 <- Sys.time()

  # use_dem = TRUE uniquement pour les modèles avec élévation (rgbie, 5 canaux)
  needs_dem <- config$in_channels >= 5

  result <- tryCatch({
    pipeline_aoi_to_landcover(
      aoi_path   = aoi_path,
      output_dir = output_dir,
      model_name = model_name,
      use_dem    = needs_dem,
      buffer_px  = BUFFER_PX
    )
  }, error = function(e) {
    message(sprintf("\n[ERREUR] %s : %s", model_name, e$message))
    NULL
  })

  dt <- round(difftime(Sys.time(), t0, units = "mins"), 1)
  timings[[model_name]] <- dt

  if (!is.null(result)) {
    results[[model_name]] <- result
    message(sprintf("\n[OK] %s terminé en %s min", model_name, dt))
  } else {
    message(sprintf("\n[ECHEC] %s après %s min", model_name, dt))
  }
}

# ==============================================================================
# Rapport comparatif
# ==============================================================================

message("\n================================================================")
message("  Rapport comparatif")
message("================================================================\n")

for (model_name in names(results)) {
  result <- results[[model_name]]
  config <- FLAIR_MODELS[[model_name]]
  dt <- timings[[model_name]]

  lc_vals <- values(result$landcover, na.rm = TRUE)
  class_counts <- table(as.integer(lc_vals))
  n_classes <- length(class_counts[names(class_counts) != "0"])

  message(sprintf("--- %s (%s min) ---", model_name, dt))
  message(sprintf("  Architecture : %s(%s)", config$decoder, config$encoder))
  message(sprintf("  Classes prédites : %d / %d", n_classes, config$n_classes))

  for (i in seq_along(class_counts)) {
    cls <- as.integer(names(class_counts)[i])
    pct <- as.numeric(class_counts[i]) / length(lc_vals) * 100
    if (cls >= 1 && cls <= 15) {
      label <- COSIA_LABELS_15[cls]
    } else if (cls == 0) {
      label <- "Non classifié"
    } else {
      label <- "?"
    }
    message(sprintf("    %2d. %-15s %5.1f%%", cls, label, pct))
  }
  message("")
}

# ==============================================================================
# Comparaison PDF (si plusieurs modèles)
# ==============================================================================

if (length(results) > 1) {
  message("--- Génération du PDF comparatif ---")
  pdf_path <- file.path(output_root, "comparaison_modeles.pdf")
  n_models <- length(results)

  pdf(pdf_path, width = 8 * n_models, height = 10)
  par(mfrow = c(1, n_models), mar = c(2, 2, 4, 4))

  for (model_name in names(results)) {
    result <- results[[model_name]]
    config <- FLAIR_MODELS[[model_name]]
    lc <- result$landcover

    lc_present <- sort(unique(na.omit(as.integer(values(lc)))))
    lc_labels <- vapply(lc_present, function(cls) {
      if (cls == 0) "Non classifié" else COSIA_LABELS_15[cls]
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
# Résumé final
# ==============================================================================

n_ok <- length(results)
n_fail <- length(models_to_test) - n_ok

message("\n================================================================")
message(sprintf("  Tests terminés : %d OK, %d échec(s)", n_ok, n_fail))
message(sprintf("  Résultats dans : %s", output_root))
message("================================================================")
