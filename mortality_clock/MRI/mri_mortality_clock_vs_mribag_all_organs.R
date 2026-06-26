# ============================================================
# Scatter plots:
# Organ MRIBAG vs organ MRI mortality clock acceleration
# Generalized RStudio-ready direct-run script for 7 organs
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(patchwork)
})

# ============================================================
# 1. Paths and organ settings
# ============================================================

# The script automatically chooses the first existing root directory.
root_candidates <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"
)

clock_root <- root_candidates[file.exists(root_candidates)][1]
if (is.na(clock_root)) {
  stop("None of the candidate WholeBodyClock root directories exists. Please set clock_root manually.")
}

# Aging-clock/MRIBAG file. The script automatically chooses the first existing file.
bag_file_candidates <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/SleepAging/data/MomoBAG.tsv",
  "/cbica/home/wenju/Reproducibile_paper/SleepAging/data/MomoBAG.tsv"
)

bag_file <- bag_file_candidates[file.exists(bag_file_candidates)][1]
if (is.na(bag_file)) {
  stop("None of the candidate MRIBAG files exists. Please set bag_file manually.")
}

# Output folder for combined cross-organ summaries and combined plot.
combined_outdir <- file.path(clock_root, "mri_mortality_clock_vs_mribag")
dir.create(combined_outdir, recursive = TRUE, showWarnings = FALSE)

# Organs expected from your WholeBodyClock folders.
organs <- c("brain", "heart", "adipose", "kidney", "liver", "pancreas", "spleen")

# Display labels.
organ_labels <- c(
  brain = "Brain",
  heart = "Heart",
  adipose = "Adipose",
  kidney = "Kidney",
  liver = "Liver",
  pancreas = "Pancreas",
  spleen = "Spleen"
)

# MRIBAG column candidates in the aging-clock file.
# Important inconsistency: brain MRIBAG is stored as Brain_PhenoBAG.
bag_column_candidates <- list(
  brain = c("Brain_PhenoBAG", "Brain_MRIBAG", "Brain_BAG"),
  heart = c("Heart_MRIBAG", "Heart_BAG"),
  adipose = c("Adipose_MRIBAG", "Adipose_BAG", "Fat_MRIBAG", "VAT_MRIBAG", "SAT_MRIBAG"),
  kidney = c("Kidney_MRIBAG", "Kidney_BAG", "Renal_MRIBAG", "Renal_BAG"),
  liver = c("Liver_MRIBAG", "Liver_BAG", "Hepatic_MRIBAG", "Hepatic_BAG"),
  pancreas = c("Pancreas_MRIBAG", "Pancreas_BAG"),
  spleen = c("Spleen_MRIBAG", "Spleen_BAG")
)

message("Clock root: ", clock_root)
message("MRIBAG file: ", bag_file)
message("Combined output folder: ", combined_outdir)

# ============================================================
# 2. Helper functions
# ============================================================

read_table_auto <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
  }
}

harmonize_participant_id <- function(df, label = "data") {
  if ("participant_id" %in% colnames(df)) {
    return(df)
  }
  if ("eid" %in% colnames(df)) {
    return(df %>% rename(participant_id = eid))
  }
  if ("id" %in% colnames(df)) {
    return(df %>% rename(participant_id = id))
  }
  stop("Could not find participant_id/eid/id column in ", label, ".")
}

find_first_existing_col <- function(df, candidates, organ) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    message("Available columns in MRIBAG file:")
    print(colnames(df))
    stop(
      "Could not find MRIBAG column for organ '", organ, "'. Tried: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[[1]]
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

# ============================================================
# 3. Read MRIBAG data
# ============================================================

bag_all <- read_table_auto(bag_file) %>%
  harmonize_participant_id(label = "MRIBAG file")

# Make participant_id type consistent for joining.
bag_all <- bag_all %>%
  mutate(participant_id = as.character(participant_id))

# ============================================================
# 4. Main function for one organ
# ============================================================

make_organ_scatter <- function(organ) {
  organ_label <- organ_labels[[organ]]
  if (is.null(organ_label)) {
    organ_label <- stringr::str_to_title(organ)
  }

  clock_dir <- file.path(clock_root, paste0(organ, "_mri_mortality_clock"))
  prediction_file <- file.path(clock_dir, paste0(organ, "_mri_mortality_clock_predictions.tsv"))

  if (!file.exists(prediction_file)) {
    warning("Skipping ", organ, ": prediction file not found: ", prediction_file)
    return(NULL)
  }

  bag_col <- find_first_existing_col(
    bag_all,
    candidates = bag_column_candidates[[organ]],
    organ = organ
  )

  mortality_year_col <- paste0(organ, "_mri_mortality_clock_acceleration_years")
  mortality_z_col <- paste0(organ, "_mri_mortality_clock_acceleration_z")
  mortality_risk_col <- paste0(organ, "_mri_mortality_risk_score")

  mort <- read_table_auto(prediction_file) %>%
    harmonize_participant_id(label = paste0(organ, " mortality prediction file")) %>%
    mutate(participant_id = as.character(participant_id))

  required_mort_cols <- c("participant_id", mortality_year_col)
  missing_mort_cols <- setdiff(required_mort_cols, colnames(mort))
  if (length(missing_mort_cols) > 0) {
    warning(
      "Skipping ", organ, ": missing mortality columns: ",
      paste(missing_mort_cols, collapse = ", ")
    )
    return(NULL)
  }

  # Keep useful columns when present.
  optional_cols <- c(
    "split", "age_at_imaging", "sex",
    mortality_z_col,
    mortality_risk_col
  )
  optional_cols <- optional_cols[optional_cols %in% colnames(mort)]

  df <- mort %>%
    select(
      participant_id,
      all_of(optional_cols),
      mortality_clock_years = all_of(mortality_year_col)
    ) %>%
    inner_join(
      bag_all %>%
        select(participant_id, mribag = all_of(bag_col)),
      by = "participant_id"
    ) %>%
    mutate(
      organ = organ,
      organ_label = organ_label,
      bag_column = bag_col,
      mribag = as.numeric(mribag),
      mortality_clock_years = as.numeric(mortality_clock_years),
      split = if ("split" %in% colnames(.)) {
        factor(split, levels = c("train", "validation", "test"))
      } else {
        factor(NA_character_, levels = c("train", "validation", "test"))
      }
    ) %>%
    filter(
      is.finite(mribag),
      is.finite(mortality_clock_years)
    )

  if (nrow(df) < 10) {
    warning("Skipping ", organ, ": merged N < 10.")
    return(NULL)
  }

  message("============================================================")
  message("Organ: ", organ_label)
  message("MRIBAG column: ", bag_col)
  message("Merged N = ", nrow(df))
  if ("split" %in% colnames(df)) {
    message("Participants by split:")
    print(table(df$split, useNA = "ifany"))
  }

  # ------------------------------------------------------------
  # Correlation statistics
  # ------------------------------------------------------------

  pearson_test <- cor.test(df$mribag, df$mortality_clock_years, method = "pearson")
  spearman_test <- cor.test(df$mribag, df$mortality_clock_years, method = "spearman", exact = FALSE)
  lm_fit <- lm(mortality_clock_years ~ mribag, data = df)
  lm_summary <- summary(lm_fit)

  pearson_r <- unname(pearson_test$estimate)
  pearson_p <- pearson_test$p.value
  spearman_rho <- unname(spearman_test$estimate)
  spearman_p <- spearman_test$p.value
  lm_beta <- coef(lm_fit)[["mribag"]]
  lm_r2 <- lm_summary$r.squared

  cor_tbl <- tibble(
    organ = organ,
    organ_label = organ_label,
    n = nrow(df),
    bag_column = bag_col,
    mortality_clock_column = mortality_year_col,
    pearson_r = pearson_r,
    pearson_p = pearson_p,
    spearman_rho = spearman_rho,
    spearman_p = spearman_p,
    lm_beta_mortality_clock_years_per_MRIBAG_year = lm_beta,
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

  x_min <- min(df$mribag, na.rm = TRUE)
  x_max <- max(df$mribag, na.rm = TRUE)
  y_min <- min(df$mortality_clock_years, na.rm = TRUE)
  y_max <- max(df$mortality_clock_years, na.rm = TRUE)

  x_range <- safe_range(x_min, x_max)
  y_range <- safe_range(y_min, y_max)

  y_upper_plot <- y_max + 0.20 * y_range
  y_lower_plot <- y_min - 0.04 * y_range

  annot_x <- x_min + 0.04 * x_range
  annot_y <- y_max + 0.14 * y_range

  p <- ggplot(
    df,
    aes(
      x = mribag,
      y = mortality_clock_years
    )
  ) +
    geom_point(
      color = point_color,
      alpha = 0.55,
      size = 2.3,
      shape = 16,
      stroke = 0
    ) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      linewidth = 1.8,
      color = line_color
    ) +
    annotate(
      "text",
      x = annot_x,
      y = annot_y,
      label = stat_text,
      hjust = 0,
      vjust = 1,
      size = 5.2,
      color = stat_color,
      fontface = "italic"
    ) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(y_lower_plot, y_upper_plot),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = paste0(organ_label, " MRIBAG, years"),
      y = paste0(organ_label, " MRI mortality clock, years")
    ) +
    theme_elegant(base_size = 15)

  # ------------------------------------------------------------
  # Save per-organ outputs
  # ------------------------------------------------------------

  # Use Brain_MRIBAG in output names even though the source column is Brain_PhenoBAG.
  out_prefix <- file.path(
    clock_dir,
    paste0("scatter_", organ_label, "_MRIBAG_vs_", organ, "_mri_mortality_clock_acceleration_years")
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
    organ = organ,
    plot = p,
    data = df,
    stats = cor_tbl,
    files = tibble(
      organ = organ,
      plot_pdf = out_pdf,
      plot_png = out_png,
      data_tsv = out_tsv,
      stat_tsv = out_stat
    )
  )
}

# ============================================================
# 5. Run all organs
# ============================================================

results <- purrr::map(organs, function(x) {
  tryCatch(
    make_organ_scatter(x),
    error = function(e) {
      warning("Failed for organ '", x, "': ", conditionMessage(e))
      NULL
    }
  )
})

results <- purrr::compact(results)

if (length(results) == 0) {
  stop("No organ scatter plots were generated. Please check input files and column names.")
}

all_stats <- purrr::map_dfr(results, "stats")
all_data <- purrr::map_dfr(results, "data")
all_files <- purrr::map_dfr(results, "files")

stats_file <- file.path(combined_outdir, "all_organs_mribag_vs_mri_mortality_clock_correlation_stats.tsv")
data_file <- file.path(combined_outdir, "all_organs_mribag_vs_mri_mortality_clock_data.tsv")
files_file <- file.path(combined_outdir, "all_organs_mribag_vs_mri_mortality_clock_output_files.tsv")

readr::write_tsv(all_stats, stats_file)
readr::write_tsv(all_data, data_file)
readr::write_tsv(all_files, files_file)

message("============================================================")
message("Saved combined outputs:")
message("  ", stats_file)
message("  ", data_file)
message("  ", files_file)

# ============================================================
# 6. Combined multi-panel figure
# ============================================================

plot_list <- purrr::map(results, "plot")
plot_names <- purrr::map_chr(results, "organ")

# Add small organ labels as panel titles through patchwork annotation wrappers.
plot_list_named <- purrr::map2(plot_list, plot_names, function(p, organ) {
  organ_label <- organ_labels[[organ]]
  p + ggtitle(organ_label) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.02,
        size = 13,
        color = "#17202A"
      )
    )
})

combined_plot <- wrap_plots(plot_list_named, ncol = 3) +
  plot_annotation(
    title = "MRIBAG versus MRI mortality-clock acceleration across organs",
    subtitle = "Brain MRIBAG is read from Brain_PhenoBAG; all other organs use organ-specific MRIBAG columns when available.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, color = "#17202A"),
      plot.subtitle = element_text(size = 11, color = "#566573")
    )
  )

combined_pdf <- file.path(combined_outdir, "all_organs_mribag_vs_mri_mortality_clock_scatter.pdf")
combined_png <- file.path(combined_outdir, "all_organs_mribag_vs_mri_mortality_clock_scatter.png")

print(combined_plot)

if (capabilities("cairo")) {
  ggsave(filename = combined_pdf, plot = combined_plot, width = 16, height = 13, device = cairo_pdf)
} else {
  ggsave(filename = combined_pdf, plot = combined_plot, width = 16, height = 13)
}

ggsave(filename = combined_png, plot = combined_plot, width = 16, height = 13, dpi = 350)

message("Saved combined figure:")
message("  ", combined_pdf)
message("  ", combined_png)

message("\n===== Combined correlation summary =====")
print(all_stats)

message("\nInterpretation guide:")
message("  - Pearson R quantifies linear agreement between MRIBAG and mortality-clock acceleration.")
message("  - R^2 near 0 means the age-prediction clock and mortality-proximity clock are largely distinct.")
message("  - Significant P values can occur with very large N even when the effect size is biologically small.")
