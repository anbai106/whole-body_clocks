# ============================================================
# Scatter plots:
# Aging clock versus mortality-clock acceleration
# MRI + proteomics + metabolomics mortality clocks
# Generalized RStudio-ready direct-run script for 23 expected clocks
# ============================================================
.libPaths('/gpfs/fs001/cbica/home/wenju/R/x86_64-pc-linux-gnu-library/4.3')
suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(patchwork)
})

# ============================================================
# 1. Paths and settings
# ============================================================

# The script automatically chooses the first existing root directory.
clock_root_candidates <- c(
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock",
  getwd()
)

clock_root <- clock_root_candidates[file.exists(clock_root_candidates)][1]
if (is.na(clock_root) || is.null(clock_root)) {
  stop("None of the candidate WholeBodyClock root directories exists. Please set clock_root manually.")
}

# Aging-clock file(s). If multiple files exist, they are full-joined by participant_id.
# MomoBAG.tsv is expected to contain MRIBAG/ProtBAG/MetBAG columns in your workflow.
aging_clock_file_candidates <- c(
  "/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv",
  "/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv",
  file.path(clock_root, "MomoBAG.tsv"),
  file.path(clock_root, "all_BAG.tsv"),
  file.path(clock_root, "all_BAGs.tsv"),
  file.path(clock_root, "all_aging_clocks.tsv"),
  file.path(clock_root, "all_aging_clock_predictions.tsv")
)

aging_clock_files <- unique(aging_clock_file_candidates[file.exists(aging_clock_file_candidates)])
if (length(aging_clock_files) == 0) {
  stop("None of the candidate aging-clock files exists. Please set aging_clock_files manually.")
}

# Output folder for combined summaries and combined plots.
combined_outdir <- file.path(clock_root, "all_mortality_clock_vs_aging_clock")
dir.create(combined_outdir, recursive = TRUE, showWarnings = FALSE)

# If a year-scale mortality-clock acceleration is unavailable, this script can fall back to z-scale.
# Keep TRUE to avoid losing clocks with missing year-scale transforms.
allow_z_fallback <- TRUE

# For the combined/main figures, use the held-out test split only.
# This keeps the aging-clock versus mortality-clock comparison aligned with external
# model-performance reporting and avoids visually mixing training/validation data.
# Set to NULL to use all splits.
analysis_split <- "test"
analysis_split_label <- ifelse(is.null(analysis_split), "all splits", paste0(analysis_split, " split only"))

# For very large samples, downsample points only for plotting. Correlations still use all data.
# Use Inf to plot all points.
max_points_per_clock_for_plot <- Inf

message("Clock root: ", clock_root)
message("Aging-clock files:")
message(paste0("  ", aging_clock_files, collapse = "\n"))
message("Combined output folder: ", combined_outdir)
message("Analysis subset for combined scatter and bar-plot summaries: ", analysis_split_label)

# ============================================================
# 2. Clock manifest
# ============================================================
# Naming convention notes:
# - MRI folders/files use lower-case organ names, e.g. brain_mri_mortality_clock.
# - Proteomics/metabolomics folders may start with capitalized organ/system names,
#   but prediction/output files are usually lower-case prefixes.

make_manifest <- function(clock_root) {
  mri_organs <- c("brain", "heart", "adipose", "kidney", "liver", "pancreas", "spleen")
  proteomics_organs <- c(
    "Reproductive_female", "Pulmonary", "Heart", "Brain", "Eye", "Hepatic",
    "Renal", "Reproductive_male", "Endocrine", "Immune", "Skin"
  )
  metabolomics_organs <- c("Endocrine", "Digestive", "Hepatic", "Immune", "Metabolic")

  bind_rows(
    tibble(
      modality = "MRI",
      modality_key = "mri",
      organ_folder_name = mri_organs,
      organ_key = stringr::str_to_lower(mri_organs),
      folder = paste0(stringr::str_to_lower(mri_organs), "_mri_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(mri_organs), "_mri_mortality_clock")
    ),
    tibble(
      modality = "Proteomics",
      modality_key = "proteomics",
      organ_folder_name = proteomics_organs,
      organ_key = stringr::str_to_lower(proteomics_organs),
      folder = paste0(proteomics_organs, "_proteomics_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(proteomics_organs), "_proteomics_mortality_clock")
    ),
    tibble(
      modality = "Metabolomics",
      modality_key = "metabolomics",
      organ_folder_name = metabolomics_organs,
      organ_key = stringr::str_to_lower(metabolomics_organs),
      folder = paste0(metabolomics_organs, "_metabolomics_mortality_clock"),
      prefix = paste0(stringr::str_to_lower(metabolomics_organs), "_metabolomics_mortality_clock")
    )
  ) %>%
    mutate(
      organ_label = organ_folder_name %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_to_sentence(),
      clock_label = paste(organ_label, modality),
      clock_id = paste(organ_key, modality_key, sep = "__"),
      clock_dir = file.path(clock_root, folder),
      prediction_file = file.path(clock_dir, paste0(prefix, "_predictions.tsv"))
    )
}

clock_manifest <- make_manifest(clock_root)
message("Expected mortality clocks: ", nrow(clock_manifest))

# ============================================================
# 3. Aging-clock column candidates
# ============================================================

organ_aliases <- list(
  brain = c("Brain"),
  heart = c("Heart", "Cardiac"),
  eye = c("Eye", "Ocular"),
  adipose = c("Adipose", "Fat", "VAT", "SAT"),
  kidney = c("Kidney", "Renal"),
  renal = c("Renal", "Kidney"),
  liver = c("Liver", "Hepatic"),
  hepatic = c("Hepatic", "Liver"),
  pancreas = c("Pancreas"),
  spleen = c("Spleen"),
  pulmonary = c("Pulmonary", "Lung"),
  endocrine = c("Endocrine"),
  immune = c("Immune"),
  skin = c("Skin"),
  digestive = c("Digestive", "GI", "Gastrointestinal"),
  metabolic = c("Metabolic"),
  reproductive_female = c("Reproductive_female", "Reproductive female", "Female_reproductive", "Female reproductive", "ReproductiveFemale"),
  reproductive_male = c("Reproductive_male", "Reproductive male", "Male_reproductive", "Male reproductive", "ReproductiveMale")
)

modality_suffixes <- list(
  mri = c("MRIBAG", "MRI_BAG", "MRI_BAG_years", "ImagingBAG", "Imaging_BAG", "PhenoBAG", "PhenotypeBAG", "BAG"),
  proteomics = c("ProtBAG", "ProteinBAG", "Protein_BAG", "ProteomicsBAG", "Proteomics_BAG", "ProteomeBAG", "Proteome_BAG", "OlinkBAG", "Olink_BAG", "BAG"),
  metabolomics = c("MetBAG", "MetabolomicsBAG", "Metabolomics_BAG", "MetabolomeBAG", "Metabolome_BAG", "NMRBAG", "NMR_BAG", "BAG")
)

make_aging_clock_candidates <- function(organ_key, modality_key) {
  aliases <- organ_aliases[[organ_key]]
  if (is.null(aliases)) {
    aliases <- c(
      organ_key,
      stringr::str_to_title(organ_key),
      stringr::str_replace_all(organ_key, "_", " "),
      stringr::str_replace_all(stringr::str_to_title(organ_key), "_", " ")
    )
  }

  # Add common formatting variants.
  aliases <- unique(c(
    aliases,
    stringr::str_replace_all(aliases, " ", "_"),
    stringr::str_replace_all(aliases, "_", " "),
    stringr::str_to_title(aliases),
    stringr::str_to_sentence(aliases)
  ))

  suffixes <- modality_suffixes[[modality_key]]
  if (is.null(suffixes)) suffixes <- c("BAG")

  candidates <- unique(c(
    as.vector(outer(aliases, suffixes, paste, sep = "_")),
    as.vector(outer(aliases, suffixes, paste, sep = ""))
  ))

  # Manual high-priority candidates and known inconsistencies.
  manual <- list(
    brain__mri = c("Brain_PhenoBAG", "Brain_MRIBAG", "Brain_BAG"),
    heart__mri = c("Heart_MRIBAG", "Heart_BAG", "Cardiac_MRIBAG"),
    adipose__mri = c("Adipose_MRIBAG", "Adipose_BAG", "Fat_MRIBAG", "VAT_MRIBAG", "SAT_MRIBAG"),
    kidney__mri = c("Kidney_MRIBAG", "Kidney_BAG", "Renal_MRIBAG", "Renal_BAG"),
    liver__mri = c("Liver_MRIBAG", "Liver_BAG", "Hepatic_MRIBAG", "Hepatic_BAG"),
    pancreas__mri = c("Pancreas_MRIBAG", "Pancreas_BAG"),
    spleen__mri = c("Spleen_MRIBAG", "Spleen_BAG"),

    brain__proteomics = c("Brain_ProtBAG", "Brain_ProteinBAG", "Brain_ProteomicsBAG", "Brain_Proteomics_BAG"),
    heart__proteomics = c("Heart_ProtBAG", "Heart_ProteinBAG", "Heart_ProteomicsBAG", "Heart_Proteomics_BAG"),
    eye__proteomics = c("Eye_ProtBAG", "Eye_ProteinBAG", "Eye_ProteomicsBAG", "Eye_Proteomics_BAG"),
    hepatic__proteomics = c("Hepatic_ProtBAG", "Liver_ProtBAG", "Hepatic_ProteinBAG", "Hepatic_ProteomicsBAG"),
    renal__proteomics = c("Renal_ProtBAG", "Kidney_ProtBAG", "Renal_ProteinBAG", "Renal_ProteomicsBAG"),
    pulmonary__proteomics = c("Pulmonary_ProtBAG", "Lung_ProtBAG", "Pulmonary_ProteinBAG", "Pulmonary_ProteomicsBAG"),
    endocrine__proteomics = c("Endocrine_ProtBAG", "Endocrine_ProteinBAG", "Endocrine_ProteomicsBAG"),
    immune__proteomics = c("Immune_ProtBAG", "Immune_ProteinBAG", "Immune_ProteomicsBAG"),
    skin__proteomics = c("Skin_ProtBAG", "Skin_ProteinBAG", "Skin_ProteomicsBAG"),
    reproductive_female__proteomics = c("Reproductive_female_ProtBAG", "Female_reproductive_ProtBAG", "ReproductiveFemale_ProtBAG", "Reproductive_female_ProteinBAG"),
    reproductive_male__proteomics = c("Reproductive_male_ProtBAG", "Male_reproductive_ProtBAG", "ReproductiveMale_ProtBAG", "Reproductive_male_ProteinBAG"),

    endocrine__metabolomics = c("Endocrine_MetBAG", "Endocrine_MetabolomicsBAG", "Endocrine_Metabolomics_BAG"),
    digestive__metabolomics = c("Digestive_MetBAG", "Digestive_MetabolomicsBAG", "Digestive_Metabolomics_BAG"),
    hepatic__metabolomics = c("Hepatic_MetBAG", "Liver_MetBAG", "Hepatic_MetabolomicsBAG", "Hepatic_Metabolomics_BAG"),
    immune__metabolomics = c("Immune_MetBAG", "Immune_MetabolomicsBAG", "Immune_Metabolomics_BAG"),
    metabolic__metabolomics = c("Metabolic_MetBAG", "Metabolic_MetabolomicsBAG", "Metabolic_Metabolomics_BAG")
  )

  key <- paste(organ_key, modality_key, sep = "__")
  unique(c(manual[[key]], candidates))
}

# ============================================================
# 4. Helper functions
# ============================================================

read_table_auto <- function(path) {
  if (!file.exists(path)) stop("File does not exist: ", path)
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
}

harmonize_participant_id <- function(df, label = "data") {
  if ("participant_id" %in% colnames(df)) return(df)
  if ("eid" %in% colnames(df)) return(df %>% rename(participant_id = eid))
  if ("id" %in% colnames(df)) return(df %>% rename(participant_id = id))
  if ("IID" %in% colnames(df)) return(df %>% rename(participant_id = IID))
  if ("FID" %in% colnames(df)) return(df %>% rename(participant_id = FID))
  stop("Could not find participant_id/eid/id/IID/FID column in ", label, ".")
}

find_first_existing_col <- function(df, candidates, clock_label) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  # Exact match first.
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[[1]])

  # Case-insensitive exact match.
  lower_map <- setNames(colnames(df), tolower(colnames(df)))
  lower_candidates <- tolower(candidates)
  hit_lower <- lower_candidates[lower_candidates %in% names(lower_map)]
  if (length(hit_lower) > 0) return(lower_map[[hit_lower[[1]]]])

  # Helpful suggestions.
  bag_like <- colnames(df)[stringr::str_detect(tolower(colnames(df)), "bag|age_gap|agegap|clock")]
  message("Available BAG-like columns in aging-clock file for debugging:")
  print(bag_like)

  stop(
    "Could not find aging-clock column for ", clock_label, ". Tried: ",
    paste(candidates, collapse = ", ")
  )
}

format_p <- function(p) {
  format.pval(p, digits = 3, eps = 1e-300)
}

safe_range <- function(x_min, x_max) {
  r <- x_max - x_min
  if (!is.finite(r) || r == 0) {
    r <- max(abs(c(x_min, x_max)), 1)
  }
  r
}

theme_elegant <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      axis.title = element_text(face = "bold", color = "#222222"),
      axis.text = element_text(color = "#333333"),
      axis.line = element_line(linewidth = 0.8, color = "#222222"),
      axis.ticks = element_line(linewidth = 0.7, color = "#222222"),
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(12, 16, 10, 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

point_color <- "#4B006E"   # deep purple
line_color  <- "#000000"   # black
stat_color  <- "#3E82A8"   # blue annotation

sanitize_for_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# ============================================================
# 5. Read and merge aging-clock data
# ============================================================

aging_clock_list <- purrr::map(aging_clock_files, function(path) {
  read_table_auto(path) %>%
    harmonize_participant_id(label = path) %>%
    mutate(participant_id = as.character(participant_id))
})

aging_all <- purrr::reduce(aging_clock_list, full_join, by = "participant_id")

message("Aging-clock dataset:")
message("  N = ", nrow(aging_all))
message("  Columns = ", ncol(aging_all))

# ============================================================
# 6. Main function for one mortality clock
# ============================================================

make_clock_scatter <- function(meta_row) {
  meta <- as.list(meta_row)

  if (!file.exists(meta$prediction_file)) {
    warning("Skipping ", meta$clock_label, ": prediction file not found: ", meta$prediction_file)
    return(NULL)
  }

  aging_candidates <- make_aging_clock_candidates(meta$organ_key, meta$modality_key)
  aging_col <- find_first_existing_col(
    aging_all,
    candidates = aging_candidates,
    clock_label = meta$clock_label
  )

  feature_token <- paste0(meta$organ_key, "_", meta$modality_key)
  mortality_year_col <- paste0(feature_token, "_mortality_clock_acceleration_years")
  mortality_z_col <- paste0(feature_token, "_mortality_clock_acceleration_z")
  mortality_risk_col <- paste0(feature_token, "_mortality_risk_score")

  mort <- read_table_auto(meta$prediction_file) %>%
    harmonize_participant_id(label = paste0(meta$clock_label, " mortality prediction file")) %>%
    mutate(participant_id = as.character(participant_id))

  # If the expected risk/acceleration columns are not found, identify equivalent columns by suffix.
  if (!mortality_year_col %in% colnames(mort)) {
    fallback_year <- grep("_mortality_clock_acceleration_years$", colnames(mort), value = TRUE)
    if (length(fallback_year) == 1) mortality_year_col <- fallback_year[[1]]
  }
  if (!mortality_z_col %in% colnames(mort)) {
    fallback_z <- grep("_mortality_clock_acceleration_z$", colnames(mort), value = TRUE)
    if (length(fallback_z) == 1) mortality_z_col <- fallback_z[[1]]
  }
  if (!mortality_risk_col %in% colnames(mort)) {
    fallback_risk <- grep("_mortality_risk_score$", colnames(mort), value = TRUE)
    if (length(fallback_risk) == 1) mortality_risk_col <- fallback_risk[[1]]
  }

  y_col <- NA_character_
  y_label_suffix <- NA_character_
  y_scale <- NA_character_

  if (mortality_year_col %in% colnames(mort)) {
    y_col <- mortality_year_col
    y_label_suffix <- "mortality clock, years"
    y_scale <- "years"
  } else if (allow_z_fallback && mortality_z_col %in% colnames(mort)) {
    y_col <- mortality_z_col
    y_label_suffix <- "mortality clock acceleration, z"
    y_scale <- "z"
  } else {
    warning(
      "Skipping ", meta$clock_label, ": missing mortality-clock acceleration columns. Tried: ",
      mortality_year_col, " and ", mortality_z_col
    )
    return(NULL)
  }

  optional_cols <- c("split", "age_at_imaging", "age_at_baseline", "sex", mortality_z_col, mortality_risk_col)
  optional_cols <- optional_cols[optional_cols %in% colnames(mort)]

  df <- mort %>%
    transmute(
      participant_id,
      across(all_of(optional_cols)),
      mortality_clock_value = as.numeric(.data[[y_col]])
    ) %>%
    inner_join(
      aging_all %>%
        transmute(
          participant_id,
          aging_clock_value = as.numeric(.data[[aging_col]])
        ),
      by = "participant_id"
    ) %>%
    mutate(
      clock_id = meta$clock_id,
      modality = meta$modality,
      modality_key = meta$modality_key,
      organ_key = meta$organ_key,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      aging_clock_column = aging_col,
      mortality_clock_column = y_col,
      mortality_clock_scale = y_scale,
      split = if ("split" %in% colnames(.)) {
        factor(split, levels = c("train", "validation", "test"))
      } else {
        factor(NA_character_, levels = c("train", "validation", "test"))
      }
    ) %>%
    filter(
      is.finite(aging_clock_value),
      is.finite(mortality_clock_value)
    )

  # Use only the held-out test split for the main combined scatter/bar outputs.
  # The original split is retained in the output data table.
  n_all_splits <- nrow(df)
  if (!is.null(analysis_split) && "split" %in% colnames(df)) {
    df <- df %>% filter(split == analysis_split)
  }

  if (nrow(df) < 10) {
    warning(
      "Skipping ", meta$clock_label, ": analysis subset N < 10 after applying analysis_split=",
      ifelse(is.null(analysis_split), "NULL", analysis_split), "."
    )
    return(NULL)
  }

  message("============================================================")
  message("Clock: ", meta$clock_label)
  message("Aging-clock column: ", aging_col)
  message("Mortality-clock column: ", y_col)
  message("Merged N across all splits = ", n_all_splits)
  message("Analysis N used for statistics/plotting = ", nrow(df), " (", analysis_split_label, ")")
  if ("split" %in% colnames(df)) {
    print(table(df$split, useNA = "ifany"))
  }

  # ------------------------------------------------------------
  # Correlation statistics
  # ------------------------------------------------------------

  pearson_test <- cor.test(df$aging_clock_value, df$mortality_clock_value, method = "pearson")
  spearman_test <- cor.test(df$aging_clock_value, df$mortality_clock_value, method = "spearman", exact = FALSE)
  lm_fit <- lm(mortality_clock_value ~ aging_clock_value, data = df)
  lm_summary <- summary(lm_fit)

  pearson_r <- unname(pearson_test$estimate)
  pearson_p <- pearson_test$p.value
  spearman_rho <- unname(spearman_test$estimate)
  spearman_p <- spearman_test$p.value
  lm_beta <- coef(lm_fit)[["aging_clock_value"]]
  lm_r2 <- lm_summary$r.squared

  cor_tbl <- tibble(
    clock_id = meta$clock_id,
    modality = meta$modality,
    organ_key = meta$organ_key,
    organ_label = meta$organ_label,
    clock_label = meta$clock_label,
    n = nrow(df),
    n_all_splits_before_subset = n_all_splits,
    analysis_split = ifelse(is.null(analysis_split), "all", analysis_split),
    aging_clock_column = aging_col,
    mortality_clock_column = y_col,
    mortality_clock_scale = y_scale,
    pearson_r = pearson_r,
    pearson_p = pearson_p,
    spearman_rho = spearman_rho,
    spearman_p = spearman_p,
    lm_beta_mortality_clock_per_aging_clock_year = lm_beta,
    lm_r2 = lm_r2
  )

  print(cor_tbl)

  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------

  stat_text <- paste0(
    "R = ", round(pearson_r, 2),
    "; P = ", format_p(pearson_p),
    "; R\u00b2 = ", round(lm_r2, 3)
  )

  x_min <- min(df$aging_clock_value, na.rm = TRUE)
  x_max <- max(df$aging_clock_value, na.rm = TRUE)
  y_min <- min(df$mortality_clock_value, na.rm = TRUE)
  y_max <- max(df$mortality_clock_value, na.rm = TRUE)

  x_range <- safe_range(x_min, x_max)
  y_range <- safe_range(y_min, y_max)

  y_upper_plot <- y_max + 0.20 * y_range
  y_lower_plot <- y_min - 0.04 * y_range

  annot_x <- x_min + 0.04 * x_range
  annot_y <- y_max + 0.14 * y_range

  plot_df <- df
  if (is.finite(max_points_per_clock_for_plot) && nrow(plot_df) > max_points_per_clock_for_plot) {
    set.seed(2026)
    plot_df <- plot_df %>% slice_sample(n = max_points_per_clock_for_plot)
  }

  p <- ggplot(plot_df, aes(x = aging_clock_value, y = mortality_clock_value)) +
    geom_point(
      color = point_color,
      alpha = 0.55,
      size = 2.1,
      shape = 16,
      stroke = 0
    ) +
    geom_smooth(
      data = df,
      method = "lm",
      se = FALSE,
      linewidth = 1.6,
      color = line_color
    ) +
    annotate(
      "text",
      x = annot_x,
      y = annot_y,
      label = stat_text,
      hjust = 0,
      vjust = 1,
      size = 4.9,
      color = stat_color,
      fontface = "italic"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(
      limits = c(y_lower_plot, y_upper_plot),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = paste0(meta$organ_label, " ", meta$modality, " aging clock, years"),
      y = paste0(meta$organ_label, " ", meta$modality, " ", y_label_suffix)
    ) +
    theme_elegant(base_size = 14)

  # ------------------------------------------------------------
  # Save per-clock outputs
  # ------------------------------------------------------------

  # Keep the old MRI brain naming convention: source column can be Brain_PhenoBAG,
  # but output filenames use Brain_MRIBAG for consistency with other MRI clocks.
  aging_label_for_file <- if (meta$clock_id == "brain__mri" && aging_col == "Brain_PhenoBAG") {
    "Brain_MRIBAG"
  } else {
    sanitize_for_filename(aging_col)
  }
  split_suffix <- ifelse(is.null(analysis_split), "", paste0("_", analysis_split))
  out_prefix <- file.path(
    meta$clock_dir,
    paste0("scatter_", aging_label_for_file, "_vs_", meta$prefix, "_acceleration", split_suffix)
  )

  out_pdf <- paste0(out_prefix, ".pdf")
  out_png <- paste0(out_prefix, ".png")
  out_tsv <- paste0(out_prefix, "_data.tsv")
  out_stat <- paste0(out_prefix, "_correlation_stats.tsv")

  readr::write_tsv(df, out_tsv)
  readr::write_tsv(cor_tbl, out_stat)

  if (capabilities("cairo")) {
    ggsave(filename = out_pdf, plot = p, width = 5.4, height = 4.6, device = cairo_pdf)
  } else {
    ggsave(filename = out_pdf, plot = p, width = 5.4, height = 4.6)
  }
  ggsave(filename = out_png, plot = p, width = 5.4, height = 4.6, dpi = 400)

  message("Saved:")
  message("  ", out_pdf)
  message("  ", out_png)
  message("  ", out_tsv)
  message("  ", out_stat)

  list(
    clock_id = meta$clock_id,
    plot = p,
    data = df,
    stats = cor_tbl,
    files = tibble(
      clock_id = meta$clock_id,
      modality = meta$modality,
      organ_label = meta$organ_label,
      clock_label = meta$clock_label,
      plot_pdf = out_pdf,
      plot_png = out_png,
      data_tsv = out_tsv,
      stat_tsv = out_stat
    )
  )
}

# ============================================================
# 7. Run all mortality clocks
# ============================================================

results <- purrr::map(seq_len(nrow(clock_manifest)), function(i) {
  tryCatch(
    make_clock_scatter(clock_manifest[i, ]),
    error = function(e) {
      warning("Failed for clock '", clock_manifest$clock_label[[i]], "': ", conditionMessage(e))
      NULL
    }
  )
})

results <- purrr::compact(results)

if (length(results) == 0) {
  stop("No scatter plots were generated. Please check prediction files and aging-clock column names.")
}

all_stats <- purrr::map_dfr(results, "stats") %>%
  mutate(
    pearson_q = p.adjust(pearson_p, method = "BH"),
    spearman_q = p.adjust(spearman_p, method = "BH")
  )
all_data <- purrr::map_dfr(results, "data")
all_files <- purrr::map_dfr(results, "files")

stats_file <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_correlation_stats.tsv")
data_file <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_data.tsv")
files_file <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_output_files.tsv")

readr::write_tsv(all_stats, stats_file)
readr::write_tsv(all_data, data_file)
readr::write_tsv(all_files, files_file)

message("============================================================")
message("Saved combined outputs:")
message("  ", stats_file)
message("  ", data_file)
message("  ", files_file)

# ============================================================
# 8. Combined multi-panel figure
# ============================================================

# Order by modality and then by absolute Pearson R within modality.
clock_order_tbl <- all_stats %>%
  mutate(modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics"))) %>%
  arrange(modality, desc(abs(pearson_r))) %>%
  mutate(clock_plot = paste0(clock_label, "  "))

clock_levels <- clock_order_tbl$clock_plot

plot_list <- purrr::map(results, "plot")
plot_meta <- purrr::map_dfr(results, "stats") %>%
  left_join(clock_order_tbl %>% select(clock_id, clock_plot), by = "clock_id") %>%
  mutate(clock_plot = factor(clock_plot, levels = clock_levels))

# Reorder the individual plot list according to clock_order_tbl.
plot_lookup <- setNames(plot_list, purrr::map_chr(results, "clock_id"))
plot_list_ordered <- purrr::map(clock_order_tbl$clock_id, ~ plot_lookup[[.x]])

plot_list_named <- purrr::map2(plot_list_ordered, seq_len(nrow(clock_order_tbl)), function(p, idx) {
  label <- clock_order_tbl$clock_label[[idx]]
  p + ggtitle(label) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.02, size = 12, color = "#17202A")
    )
})

combined_plot <- wrap_plots(plot_list_named, ncol = 4) +
  plot_annotation(
    title = "Aging clocks versus mortality-clock acceleration across modalities",
    subtitle = glue("Held-out {analysis_split_label}. MRI uses MRIBAG columns, with Brain MRI read from Brain_PhenoBAG when present; proteomics uses ProtBAG-like columns; metabolomics uses MetBAG-like columns."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, color = "#17202A"),
      plot.subtitle = element_text(size = 11, color = "#566573")
    )
  )

combined_pdf <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_scatter.pdf")
combined_png <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_scatter.png")

print(combined_plot)

combined_height <- max(13, ceiling(length(plot_list_named) / 4) * 4.2)
if (capabilities("cairo")) {
  ggsave(filename = combined_pdf, plot = combined_plot, width = 20, height = combined_height, device = cairo_pdf)
} else {
  ggsave(filename = combined_pdf, plot = combined_plot, width = 20, height = combined_height)
}
ggsave(filename = combined_png, plot = combined_plot, width = 20, height = combined_height, dpi = 350)

message("Saved combined scatter figure:")
message("  ", combined_pdf)
message("  ", combined_png)

# ============================================================
# 9. Compact correlation overview figure
# ============================================================

overview_tbl <- all_stats %>%
  mutate(
    modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),
    clock_label = factor(clock_label, levels = clock_order_tbl$clock_label),
    direction = if_else(pearson_r >= 0, "Positive", "Negative")
  )

p_overview <- overview_tbl %>%
  ggplot(aes(x = reorder(clock_label, pearson_r), y = pearson_r, fill = modality)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7B7D7D", linewidth = 0.5) +
  geom_col(width = 0.74, alpha = 0.92) +
  coord_flip() +
  scale_fill_manual(values = c("MRI" = "#2E86AB", "Proteomics" = "#8E44AD", "Metabolomics" = "#D35400")) +
  labs(
    title = "Test-set correlation between aging clocks and mortality-clock acceleration",
    subtitle = glue("Pearson R in the held-out {analysis_split_label}; positive values indicate aging-like phenotype aligns with mortality-proximity phenotype"),
    x = NULL,
    y = "Pearson R",
    fill = "Modality"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", color = "#17202A"),
    plot.subtitle = element_text(color = "#566573"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

overview_pdf <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_correlation_overview.pdf")
overview_png <- file.path(combined_outdir, "all_mortality_clocks_vs_aging_clocks_test_correlation_overview.png")

if (capabilities("cairo")) {
  ggsave(filename = overview_pdf, plot = p_overview, width = 9, height = 8.5, device = cairo_pdf)
} else {
  ggsave(filename = overview_pdf, plot = p_overview, width = 9, height = 8.5)
}
ggsave(filename = overview_png, plot = p_overview, width = 9, height = 8.5, dpi = 350)

message("Saved correlation overview figure:")
message("  ", overview_pdf)
message("  ", overview_png)

message("\n===== Combined correlation summary =====")
print(
  all_stats %>%
    arrange(factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")), desc(abs(pearson_r))) %>%
    select(
      modality,
      organ_label,
      n,
      aging_clock_column,
      mortality_clock_column,
      mortality_clock_scale,
      pearson_r,
      pearson_p,
      pearson_q,
      spearman_rho,
      lm_beta_mortality_clock_per_aging_clock_year,
      lm_r2
    )
)

message("\nInterpretation guide:")
message("  - Pearson R quantifies linear agreement between the aging clock and mortality-clock acceleration in the selected analysis subset.")
message("  - R^2 near 0 means the age-prediction clock and mortality-proximity clock are largely distinct.")
message("  - Positive R means aging-like acceleration aligns with higher mortality-clock acceleration.")
message("  - Negative R means mortality-proximity biology moves opposite to the normative aging-clock axis.")
message("  - Significant P values can occur with very large N even when the effect size is modest; prioritize R and R^2.")
