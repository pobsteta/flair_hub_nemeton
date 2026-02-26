#!/usr/bin/env Rscript
# ==============================================================================
# 01_download_flair_hub.R
# Téléchargement du dataset FLAIR-HUB depuis Hugging Face
# + Chargement des différentes modalités
#
# Dataset HF : IGNF/FLAIR-HUB
# https://huggingface.co/datasets/IGNF/FLAIR-HUB
#
# Structure du dataset :
#   DOMAIN_SENSOR_DATATYPE / ROI / PATCH
#   Patches : 512x512 pixels à 0.2m VHR (102.4m x 102.4m au sol)
#
# 6 modalités alignées :
#   - AERIAL RGBI     : 0.2m, 4 bandes (R, G, B, NIR), UInt8
#   - AERIAL-RLT PAN  : 0.4m, 1 bande (Panchromatique historique 1950s), UInt8
#   - DEM ELEV        : 0.2m, 2 bandes (DSM, DTM), Float32
#   - SPOT RGBI       : 1.6m, 4 bandes (R, G, B, NIR), UInt16
#   - Sentinel-2 SITS : 10m, séries temporelles
#   - Sentinel-1 SITS : 10m, séries temporelles
#
# 2 supervisions :
#   - LABEL-COSIA : 19 classes occupation du sol (photo-interprétation)
#   - LABEL-LPIS  : 23 classes cultures (RPG / LPIS)
# ==============================================================================

# --- Installation des packages nécessaires ---
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("httr2")
install_if_missing("jsonlite")
install_if_missing("terra")
install_if_missing("sf")
install_if_missing("curl")
install_if_missing("fs")

library(httr2)
library(jsonlite)
library(curl)
library(fs)

# ==============================================================================
# Configuration
# ==============================================================================

# --- Hugging Face ---
HF_REPO_ID  <- "IGNF/FLAIR-HUB"
HF_API_URL  <- "https://huggingface.co/api/datasets"
HF_TOKEN    <- Sys.getenv("HF_TOKEN", unset = "")

# --- Répertoires ---
DATA_DIR    <- file.path(getwd(), "data")
DATA_DIR_HF <- file.path(DATA_DIR, "flair_hub")
dir_create(c(DATA_DIR, DATA_DIR_HF))

# --- Résolutions par modalité ---
RES_AERIAL  <- 0.2   # Aérien RGBI
RES_AERIAL_RLT <- 0.4  # Aérien historique (rééchantillonné)
RES_DEM     <- 0.2   # MNT (DSM/DTM)
RES_SPOT    <- 1.6   # SPOT RGBI (rééchantillonné depuis 1.5m)
RES_S2      <- 10    # Sentinel-2
RES_S1      <- 10    # Sentinel-1

# --- Classes d'occupation du sol (COSIA) - 19 classes ---
COSIA_CLASSES <- data.frame(
  id = 1:19,
  label_fr = c(
    "Bâtiment",              # 1
    "Serre",                 # 2
    "Piscine",               # 3
    "Surface imperméable",   # 4
    "Surface perméable",     # 5
    "Sol nu",                # 6
    "Eau",                   # 7
    "Neige",                 # 8
    "Végétation herbacée",   # 9
    "Terre agricole",        # 10
    "Terre labourée",        # 11
    "Vigne",                 # 12
    "Verger",                # 13
    "Feuillu",               # 14
    "Conifère",              # 15
    "Lande",                 # 16
    "Ligneux mélangé",       # 17
    "Fleur / Garrigue",      # 18
    "Non classé"             # 19
  ),
  label_en = c(
    "Building",
    "Greenhouse",
    "Swimming pool",
    "Impervious surface",
    "Pervious surface",
    "Bare soil",
    "Water",
    "Snow",
    "Herbaceous vegetation",
    "Agricultural land",
    "Plowed land",
    "Vineyard",
    "Orchard",
    "Deciduous",
    "Coniferous",
    "Brushwood",
    "Mixed woodland",
    "Flower / Garrigue",
    "Unclassified"
  ),
  color = c(
    "#db0e9a", "#938e7b", "#f80c00", "#a97101", "#1553ae",
    "#194a26", "#46e483", "#f3a60d", "#660082", "#55ff00",
    "#fff30d", "#e4df7c", "#3de6eb", "#ffffff", "#8ab3a0",
    "#6b714f", "#c5dc42", "#9999ff", "#000000"
  ),
  stringsAsFactors = FALSE
)

# --- Classes de cultures (LPIS) - 23 classes niveau basique ---
LPIS_CLASSES <- data.frame(
  id = 0:22,
  label_fr = c(
    "Fond / Non agricole",   # 0
    "Blé tendre",            # 1
    "Maïs",                  # 2
    "Orge",                  # 3
    "Colza",                 # 4
    "Tournesol",             # 5
    "Prairie temporaire",    # 6
    "Prairie permanente",    # 7
    "Gel / Jachère",         # 8
    "Soja",                  # 9
    "Légumineuse fourragère", # 10
    "Betterave",             # 11
    "Pomme de terre",        # 12
    "Autre culture permanente", # 13
    "Culture mélangée",      # 14
    "Autre oléagineux",      # 15
    "Protéagineux",          # 16
    "Riz",                   # 17
    "Légume / Fleur",        # 18
    "Sorgho",                # 19
    "Blé dur",               # 20
    "Avoine",                # 21
    "Triticale"              # 22
  ),
  label_en = c(
    "Background / Non-agricultural",
    "Soft wheat",
    "Corn / Maize",
    "Barley",
    "Rapeseed",
    "Sunflower",
    "Temporary grassland",
    "Permanent grassland",
    "Fallow land",
    "Soy",
    "Fodder legumes",
    "Beetroots",
    "Potatoes",
    "Other permanent crops",
    "Mixed crops",
    "Other oilseed crops",
    "Protein crops",
    "Rice",
    "Vegetable / Flower",
    "Sorghum",
    "Durum wheat",
    "Oats",
    "Triticale"
  ),
  stringsAsFactors = FALSE
)

# --- Modalités disponibles ---
FLAIR_MODALITIES <- data.frame(
  modality = c("AERIAL_RGBI", "AERIAL-RLT_PAN", "DEM_ELEV",
               "SPOT_RGBI", "S2_SITS", "S1_SITS",
               "AERIAL_LABEL-COSIA", "AERIAL_LABEL-LPIS"),
  description = c(
    "Aérien RGBI 0.2m (4 bandes R,G,B,NIR, UInt8)",
    "Aérien historique 1950s PAN 0.4m (1 bande, UInt8)",
    "MNT 0.2m (2 bandes DSM,DTM, Float32)",
    "SPOT RGBI 1.6m (4 bandes R,G,B,NIR, UInt16 BOA)",
    "Sentinel-2 séries temporelles 10m",
    "Sentinel-1 séries temporelles 10m",
    "Labels occupation du sol CoSIA (19 classes)",
    "Labels cultures LPIS (23 classes)"
  ),
  resolution_m = c(0.2, 0.4, 0.2, 1.6, 10, 10, 0.2, 0.2),
  stringsAsFactors = FALSE
)

# ==============================================================================
# A. Fonctions Hugging Face (FLAIR-HUB)
# ==============================================================================

#' Lister les fichiers du dataset sur Hugging Face
#'
#' @param repo_id Identifiant du dépôt (ex: "IGNF/FLAIR-HUB")
#' @param path Chemin dans le dépôt (ex: "" pour la racine)
#' @param token Token Hugging Face (optionnel)
#' @return data.frame avec les informations des fichiers
hf_list_files <- function(repo_id = HF_REPO_ID, path = "", token = HF_TOKEN) {
  url <- paste0(HF_API_URL, "/", repo_id, "/tree/main")
  if (nchar(path) > 0) {
    url <- paste0(url, "/", path)
  }

  req <- request(url)
  if (nchar(token) > 0) {
    req <- req |> req_headers(Authorization = paste("Bearer", token))
  }

  resp <- req |>
    req_error(is_error = function(resp) FALSE) |>
    req_perform()

  if (resp_status(resp) != 200) {
    warning("Erreur API HF: ", resp_status(resp))
    return(data.frame())
  }

  content <- resp |> resp_body_json()

  files_df <- do.call(rbind, lapply(content, function(item) {
    data.frame(
      type = item$type %||% NA,
      path = item$path %||% NA,
      size = item$size %||% NA,
      oid = item$oid %||% NA,
      stringsAsFactors = FALSE
    )
  }))

  return(files_df)
}

#' Télécharger un fichier depuis Hugging Face
#'
#' @param repo_id Identifiant du dépôt
#' @param filename Chemin du fichier dans le dépôt
#' @param dest_dir Répertoire de destination local
#' @param token Token Hugging Face (optionnel)
#' @param overwrite Écraser si le fichier existe déjà
#' @return Chemin local du fichier téléchargé
hf_download_file <- function(repo_id = HF_REPO_ID, filename, dest_dir = DATA_DIR_HF,
                              token = HF_TOKEN, overwrite = FALSE) {
  local_path <- file.path(dest_dir, filename)

  if (file.exists(local_path) && !overwrite) {
    message("Fichier déjà présent: ", local_path)
    return(local_path)
  }

  dir_create(dirname(local_path))

  url <- paste0(
    "https://huggingface.co/datasets/", repo_id,
    "/resolve/main/", filename
  )

  message("Téléchargement HF: ", filename)

  headers <- list()
  if (nchar(token) > 0) {
    headers[["Authorization"]] <- paste("Bearer", token)
  }

  tryCatch({
    curl_download(
      url = url,
      destfile = local_path,
      handle = new_handle(.list = headers),
      quiet = FALSE
    )
    message("OK: ", local_path)
    return(local_path)
  }, error = function(e) {
    warning("Échec du téléchargement: ", filename, " - ", e$message)
    return(NULL)
  })
}

#' Télécharger un ensemble de fichiers depuis Hugging Face
#'
#' @param repo_id Identifiant du dépôt
#' @param file_list Liste des fichiers à télécharger
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
#' @param overwrite Écraser les fichiers existants
#' @return Vecteur des chemins locaux
hf_download_files <- function(repo_id = HF_REPO_ID, file_list,
                               dest_dir = DATA_DIR_HF, token = HF_TOKEN,
                               overwrite = FALSE) {
  paths <- character(length(file_list))
  for (i in seq_along(file_list)) {
    paths[i] <- hf_download_file(
      repo_id = repo_id,
      filename = file_list[i],
      dest_dir = dest_dir,
      token = token,
      overwrite = overwrite
    )
  }
  return(paths)
}

#' Télécharger un sous-ensemble du dataset FLAIR-HUB
#'
#' Le dataset est structuré par : DOMAIN_SENSOR_DATATYPE / ROI / PATCH
#'
#' @param domain Domaine géographique (ex: "D001", "D032", etc.)
#' @param modality Modalité : "AERIAL_RGBI", "SPOT_RGBI", "DEM_ELEV",
#'   "S2_SITS", "S1_SITS", "AERIAL-RLT_PAN",
#'   "AERIAL_LABEL-COSIA", "AERIAL_LABEL-LPIS"
#' @param n_patches Nombre de patches à télécharger
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
#' @return Vecteur des chemins locaux
download_flair_hub_subset <- function(domain = NULL,
                                       modality = "AERIAL_RGBI",
                                       n_patches = 10,
                                       dest_dir = DATA_DIR_HF,
                                       token = HF_TOKEN) {
  message("=== Téléchargement FLAIR-HUB ===")
  message("Dataset: ", HF_REPO_ID)
  message("Modalité: ", modality)

  # Lister les dossiers de premier niveau (DOMAIN_SENSOR_DATATYPE)
  root_files <- hf_list_files(path = "", token = token)

  if (nrow(root_files) == 0) {
    message("Impossible de lister le dataset. Vérifiez votre token HF.")
    return(invisible(NULL))
  }

  # Filtrer par modalité
  dirs <- root_files[root_files$type == "directory", ]
  if (nrow(dirs) > 0) {
    # Filtrer par modalité demandée
    modality_dirs <- dirs[grep(modality, dirs$path, ignore.case = TRUE), ]

    # Filtrer par domaine si spécifié
    if (!is.null(domain)) {
      modality_dirs <- modality_dirs[grep(domain, modality_dirs$path,
                                           ignore.case = TRUE), ]
    }

    if (nrow(modality_dirs) == 0) {
      message("Aucun répertoire trouvé pour la modalité: ", modality)
      message("Répertoires disponibles:")
      print(dirs$path)
      return(invisible(NULL))
    }

    message(sprintf("%d répertoire(s) trouvé(s) pour %s",
                     nrow(modality_dirs), modality))

    # Télécharger les patches
    downloaded <- character(0)
    for (dir_path in modality_dirs$path) {
      if (length(downloaded) >= n_patches) break

      roi_files <- hf_list_files(path = dir_path, token = token)
      if (nrow(roi_files) == 0) next

      # Chercher dans les sous-répertoires (ROI)
      roi_dirs <- roi_files[roi_files$type == "directory", ]
      for (roi_dir in roi_dirs$path) {
        if (length(downloaded) >= n_patches) break

        patches <- hf_list_files(path = roi_dir, token = token)
        tif_patches <- patches[grep("\\.tif$", patches$path), ]

        if (nrow(tif_patches) > 0) {
          n_to_dl <- min(nrow(tif_patches), n_patches - length(downloaded))
          to_dl <- tif_patches$path[seq_len(n_to_dl)]

          new_paths <- hf_download_files(
            file_list = to_dl,
            dest_dir = dest_dir,
            token = token
          )
          downloaded <- c(downloaded, new_paths)
        }
      }
    }

    message(sprintf("\n=== %d patch(es) téléchargé(s) ===", length(downloaded)))
    return(downloaded)
  }

  message("Structure inattendue du dataset.")
  return(invisible(NULL))
}

#' Télécharger une modalité complète pour un domaine donné
#'
#' @param domain Code du domaine (ex: "D001_2020")
#' @param modality Modalité souhaitée
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
#' @return Vecteur des chemins téléchargés
download_flair_hub_modality <- function(domain, modality = "AERIAL_RGBI",
                                         dest_dir = DATA_DIR_HF,
                                         token = HF_TOKEN) {
  message("=== Téléchargement modalité complète ===")
  message(sprintf("Domaine: %s | Modalité: %s", domain, modality))

  # Construire le chemin : DOMAIN_SENSOR_DATATYPE
  path <- paste0(domain, "_", modality)
  files <- hf_list_files(path = path, token = token)

  if (nrow(files) == 0) {
    message("Aucun fichier trouvé pour: ", path)
    return(invisible(NULL))
  }

  # Récupérer récursivement tous les fichiers .tif
  all_tifs <- character(0)
  to_explore <- files$path[files$type == "directory"]

  while (length(to_explore) > 0) {
    current <- to_explore[1]
    to_explore <- to_explore[-1]

    sub_files <- hf_list_files(path = current, token = token)
    if (nrow(sub_files) > 0) {
      tifs <- sub_files$path[grep("\\.tif$", sub_files$path)]
      all_tifs <- c(all_tifs, tifs)

      new_dirs <- sub_files$path[sub_files$type == "directory"]
      to_explore <- c(to_explore, new_dirs)
    }
  }

  if (length(all_tifs) == 0) {
    message("Aucun fichier .tif trouvé.")
    return(invisible(NULL))
  }

  message(sprintf("%d fichier(s) .tif à télécharger", length(all_tifs)))
  downloaded <- hf_download_files(file_list = all_tifs, dest_dir = dest_dir,
                                    token = token)

  message(sprintf("\n=== %d fichier(s) téléchargé(s) ===",
                   sum(!is.na(downloaded))))
  return(downloaded)
}

#' Cloner le dataset complet via git (nécessite git-lfs)
#'
#' ATTENTION : le dataset complet fait plusieurs centaines de Go.
#'
#' @param dest_dir Répertoire de destination
#' @param token Token Hugging Face
#' @return Chemin du répertoire cloné
hf_git_clone <- function(dest_dir = file.path(DATA_DIR_HF, "FLAIR-HUB"),
                          token = HF_TOKEN) {
  if (dir.exists(dest_dir)) {
    message("Le répertoire existe déjà: ", dest_dir)
    return(invisible(dest_dir))
  }

  lfs_check <- system("git lfs version", intern = TRUE, ignore.stderr = TRUE)
  if (length(lfs_check) == 0) {
    stop("git-lfs n'est pas installé. ",
         "Installez-le avec: sudo apt install git-lfs")
  }

  message("Clonage du dataset FLAIR-HUB (ATTENTION: très volumineux)...")

  clone_url <- if (nchar(token) > 0) {
    paste0("https://", token, "@huggingface.co/datasets/", HF_REPO_ID)
  } else {
    paste0("https://huggingface.co/datasets/", HF_REPO_ID)
  }

  system2("git", args = c("clone", clone_url, dest_dir))
  message("Clonage terminé: ", dest_dir)
  return(invisible(dest_dir))
}

#' Télécharger le TOY DATASET FLAIR-HUB (petit jeu de données de test)
#'
#' Télécharge et décompresse le jeu de données jouet (~quelques Mo)
#' fourni par l'IGN pour tester le pipeline sans télécharger le dataset complet.
#'
#' @param dest_dir Répertoire de destination
#' @param overwrite Écraser si déjà présent
#' @return Chemin du répertoire décompressé
download_toy_dataset <- function(dest_dir = DATA_DIR_HF, overwrite = FALSE) {
  toy_url <- paste0(
    "https://huggingface.co/datasets/IGNF/FLAIR-HUB/resolve/main/",
    "FLAIR-HUB_TOY_DATASET.zip"
  )
  zip_path <- file.path(dest_dir, "FLAIR-HUB_TOY_DATASET.zip")
  toy_dir <- file.path(dest_dir, "FLAIR-HUB_TOY_DATASET")

  if (dir.exists(toy_dir) && !overwrite) {
    message("Toy dataset déjà présent: ", toy_dir)
    return(invisible(toy_dir))
  }

  dir_create(dest_dir)
  message("=== Téléchargement du TOY DATASET FLAIR-HUB ===")
  message("URL: ", toy_url)

  headers <- list()
  token <- Sys.getenv("HF_TOKEN", unset = "")
  if (nchar(token) > 0) {
    headers[["Authorization"]] <- paste("Bearer", token)
  }

  tryCatch({
    curl_download(
      url = toy_url,
      destfile = zip_path,
      handle = new_handle(.list = headers),
      quiet = FALSE
    )
    message("Décompression...")
    unzip(zip_path, exdir = dest_dir)
    message("TOY DATASET: ", toy_dir)

    # Lister le contenu
    toy_files <- dir_ls(toy_dir, recurse = TRUE, glob = "*.tif")
    message(sprintf("  %d fichier(s) .tif trouvé(s)", length(toy_files)))

    return(invisible(toy_dir))
  }, error = function(e) {
    warning("Échec: ", e$message)
    message("Le fichier n'est peut-être pas accessible sans token.")
    message("Configurez: Sys.setenv(HF_TOKEN = 'hf_votre_token')")
    return(NULL)
  })
}

# ==============================================================================
# B. Chargement des différentes modalités
# ==============================================================================

#' Charger un patch FLAIR-HUB (générique)
#'
#' @param file_path Chemin vers le fichier .tif
#' @param modality Type de modalité pour nommer les bandes
#' @return SpatRaster
load_flair_patch <- function(file_path, modality = "auto") {
  if (!file.exists(file_path)) {
    stop("Fichier introuvable: ", file_path)
  }

  r <- terra::rast(file_path)

  # Détection automatique de la modalité depuis le nom de fichier
  if (modality == "auto") {
    fname <- toupper(basename(file_path))
    if (grepl("LABEL-COSIA|COSIA", fname)) modality <- "LABEL_COSIA"
    else if (grepl("LABEL-LPIS|LPIS", fname)) modality <- "LABEL_LPIS"
    else if (grepl("AERIAL.*RLT|RLT.*PAN", fname)) modality <- "AERIAL_RLT"
    else if (grepl("DEM|ELEV", fname)) modality <- "DEM_ELEV"
    else if (grepl("SPOT", fname)) modality <- "SPOT_RGBI"
    else if (grepl("S2|SENTINEL-2|SENTINEL2", fname)) modality <- "S2_SITS"
    else if (grepl("S1|SENTINEL-1|SENTINEL1", fname)) modality <- "S1_SITS"
    else if (grepl("AERIAL|RGBI", fname)) modality <- "AERIAL_RGBI"
    else modality <- "unknown"
  }

  # Nommer les bandes
  if (modality == "AERIAL_RGBI" && nlyr(r) >= 4) {
    names(r)[1:4] <- c("Rouge", "Vert", "Bleu", "PIR")
  } else if (modality == "SPOT_RGBI" && nlyr(r) >= 4) {
    names(r)[1:4] <- c("Rouge", "Vert", "Bleu", "PIR")
  } else if (modality == "DEM_ELEV" && nlyr(r) >= 2) {
    names(r)[1:2] <- c("DSM", "DTM")
  } else if (modality == "AERIAL_RLT" && nlyr(r) >= 1) {
    names(r)[1] <- "PAN"
  } else if (modality == "LABEL_COSIA" && nlyr(r) >= 1) {
    names(r)[1] <- "landcover"
  } else if (modality == "LABEL_LPIS" && nlyr(r) >= 1) {
    names(r)[1] <- "crop_type"
  }

  message(sprintf("FLAIR-HUB patch chargé (%s): %s", modality, basename(file_path)))
  message(sprintf("  Dimensions: %d x %d | Bandes: %d (%s)",
                   nrow(r), ncol(r), nlyr(r), paste(names(r), collapse = ", ")))
  message(sprintf("  Résolution: %.2f x %.2f m", res(r)[1], res(r)[2]))

  return(r)
}

#' Charger une image aérienne RGBI (0.2m, 4 bandes)
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster avec bandes Rouge, Vert, Bleu, PIR
load_aerial_rgbi <- function(file_path) {
  load_flair_patch(file_path, modality = "AERIAL_RGBI")
}

#' Charger une image SPOT RGBI (1.6m, 4 bandes)
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster avec bandes Rouge, Vert, Bleu, PIR
load_spot_rgbi <- function(file_path) {
  load_flair_patch(file_path, modality = "SPOT_RGBI")
}

#' Charger une série temporelle Sentinel-2
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster multi-bandes
load_sentinel2_sits <- function(file_path) {
  load_flair_patch(file_path, modality = "S2_SITS")
}

#' Charger une série temporelle Sentinel-1
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster multi-bandes
load_sentinel1_sits <- function(file_path) {
  load_flair_patch(file_path, modality = "S1_SITS")
}

#' Charger un MNT (DSM + DTM)
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster avec bandes DSM et DTM
load_dem_elev <- function(file_path) {
  load_flair_patch(file_path, modality = "DEM_ELEV")
}

#' Charger une photo aérienne historique (1950s, panchromatique)
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster
load_aerial_rlt <- function(file_path) {
  load_flair_patch(file_path, modality = "AERIAL_RLT")
}

#' Charger les labels d'occupation du sol CoSIA (19 classes)
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster
load_label_cosia <- function(file_path) {
  load_flair_patch(file_path, modality = "LABEL_COSIA")
}

#' Charger les labels de cultures LPIS (23 classes)
#'
#' @param file_path Chemin vers le fichier .tif
#' @return SpatRaster
load_label_lpis <- function(file_path) {
  load_flair_patch(file_path, modality = "LABEL_LPIS")
}

#' Charger un ROI complet avec toutes les modalités disponibles
#'
#' @param roi_dir Répertoire contenant les fichiers du ROI
#' @return Liste de SpatRasters par modalité
load_flair_roi <- function(roi_dir) {
  if (!dir.exists(roi_dir)) {
    stop("Répertoire introuvable: ", roi_dir)
  }

  tif_files <- dir_ls(roi_dir, recurse = TRUE, glob = "*.tif")
  message(sprintf("ROI: %d fichier(s) trouvé(s) dans %s",
                   length(tif_files), basename(roi_dir)))

  result <- list()
  for (f in tif_files) {
    tryCatch({
      r <- load_flair_patch(f)
      key <- tools::file_path_sans_ext(basename(f))
      result[[key]] <- r
    }, error = function(e) {
      warning("Impossible de charger: ", f, " - ", e$message)
    })
  }

  return(result)
}

#' Scanner un répertoire pour trouver les fichiers FLAIR-HUB par modalité
#'
#' @param dir_path Répertoire à scanner
#' @param recursive Recherche récursive
#' @return data.frame avec les fichiers trouvés et leur modalité
scan_flair_files <- function(dir_path = DATA_DIR_HF, recursive = TRUE) {
  tif_files <- dir_ls(dir_path, recurse = recursive, glob = "*.tif")

  if (length(tif_files) == 0) {
    message("Aucun fichier .tif trouvé dans: ", dir_path)
    return(data.frame())
  }

  # Détection de la modalité
  detect_modality <- function(path) {
    upper_path <- toupper(path)
    if (grepl("LABEL-COSIA|COSIA", upper_path)) return("LABEL_COSIA")
    if (grepl("LABEL-LPIS|LPIS", upper_path)) return("LABEL_LPIS")
    if (grepl("AERIAL.*RLT|RLT.*PAN", upper_path)) return("AERIAL_RLT")
    if (grepl("DEM|ELEV", upper_path)) return("DEM_ELEV")
    if (grepl("SPOT", upper_path)) return("SPOT_RGBI")
    if (grepl("S2|SENTINEL-2", upper_path)) return("S2_SITS")
    if (grepl("S1|SENTINEL-1", upper_path)) return("S1_SITS")
    if (grepl("AERIAL|RGBI", upper_path)) return("AERIAL_RGBI")
    return("unknown")
  }

  df <- data.frame(
    path = tif_files,
    filename = basename(tif_files),
    modality = vapply(tif_files, detect_modality, character(1)),
    size_mb = file.size(tif_files) / 1024^2,
    stringsAsFactors = FALSE
  )

  message(sprintf("Trouvé %d fichier(s) FLAIR-HUB:", nrow(df)))
  for (mod in unique(df$modality)) {
    n <- sum(df$modality == mod)
    message(sprintf("  %s: %d fichier(s)", mod, n))
  }

  return(df)
}

# ==============================================================================
# Point d'entrée principal
# ==============================================================================

if (sys.nframe() == 0) {
  message("=== FLAIR-HUB : Chargement des données ===\n")

  # --- A. Dataset Hugging Face ---
  message("--- A. Dataset Hugging Face (IGNF/FLAIR-HUB) ---")
  message("Dataset: ", HF_REPO_ID)
  message("Token HF: ", ifelse(nchar(HF_TOKEN) > 0, "configuré", "non défini"))

  root_files <- hf_list_files(token = HF_TOKEN)
  if (nrow(root_files) > 0) {
    message("\nStructure du dataset:")
    print(root_files[, c("type", "path")])
  }

  # --- B. Modalités disponibles ---
  message("\n--- B. Modalités FLAIR-HUB ---")
  print(FLAIR_MODALITIES)

  # --- C. Classes ---
  message("\n--- C. Classes d'occupation du sol (CoSIA, 19 classes) ---")
  print(COSIA_CLASSES[, c("id", "label_fr")])

  message("\n--- D. Classes de cultures (LPIS, 23 classes) ---")
  print(LPIS_CLASSES[, c("id", "label_fr")])

  # --- D. Données locales ---
  message("\n--- E. Données locales ---")
  message("Répertoire: ", DATA_DIR_HF)
  scan_flair_files()

  message("\n=== Fin ===")
}
