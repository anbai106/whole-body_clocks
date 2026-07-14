#!/usr/bin/env Rscript

# ============================================================
# Supplementary trumpet plots for 22 mortality EPOCH clocks
# Fix: TrumpetPlots requires data.table input because it uses :=
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(TrumpetPlots)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# ============================================================
# 1. Paths
# ============================================================

BASE_DIR_CANDIDATES <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  "/cbica/home/wenju/Reproducibile_paper/WholeBodyClock"
)

BASE_DIR_CANDIDATES <- BASE_DIR_CANDIDATES[dir.exists(BASE_DIR_CANDIDATES)]

if (length(BASE_DIR_CANDIDATES) == 0) {
  stop("Cannot find WholeBodyClock base directory on Mac or CUBIC.")
}

BASE_DIR <- BASE_DIR_CANDIDATES[1]

INPUT_DIR <- file.path(
  BASE_DIR,
  "Result",
  "TrumpetPlots_mortality_epoch"
)

MANIFEST_FILE <- file.path(
  INPUT_DIR,
  "TrumpetPlots_mortality_epoch_input_manifest.tsv"
)

COMBINED_INPUT_FILE <- file.path(
  INPUT_DIR,
  "TrumpetPlots_input_22_mortality_epoch_combined.tsv"
)

OUTDIR <- file.path(
  INPUT_DIR,
  "figures_22_mortality_epoch"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("Base directory: ", BASE_DIR)
message("Input directory: ", INPUT_DIR)
message("Output directory: ", OUTDIR)

# ============================================================
# 2. Expected order
# ============================================================

EXPECTED_CLOCKS <- c(
  "brain_mri_mortality_clock",
  "adipose_mri_mortality_clock",
  "heart_mri_mortality_clock",
  "kidney_mri_mortality_clock",
  "liver_mri_mortality_clock",
  "pancreas_mri_mortality_clock",
  "spleen_mri_mortality_clock",
  
  "Brain_proteomics_mortality_clock",
  "Eye_proteomics_mortality_clock",
  "Heart_proteomics_mortality_clock",
  "Hepatic_proteomics_mortality_clock",
  "Renal_proteomics_mortality_clock",
  "Pulmonary_proteomics_mortality_clock",
  "Endocrine_proteomics_mortality_clock",
  "Immune_proteomics_mortality_clock",
  "Skin_proteomics_mortality_clock",
  "Reproductive_female_proteomics_mortality_clock",
  "Reproductive_male_proteomics_mortality_clock",
  
  "Endocrine_metabolomics_mortality_clock",
  "Digestive_metabolomics_mortality_clock",
  "Hepatic_metabolomics_mortality_clock",
  "Immune_metabolomics_mortality_clock"
)

MODALITY_ORDER <- c("MRI", "Proteomics", "Metabolomics")

# ============================================================
# 3. Colors
# ============================================================

organ_colors <- c(
  "Brain"               = "#0072B2",
  "Adipose"             = "#A6761D",
  "Heart"               = "#D55E00",
  "Renal"               = "#8B5A2B",
  "Kidney"              = "#8B5A2B",
  "Hepatic"             = "#009E73",
  "Liver"               = "#009E73",
  "Pancreas"            = "#66A61E",
  "Spleen"              = "#666666",
  "Eye"                 = "#56B4E9",
  "Pulmonary"           = "#7B61A8",
  "Endocrine"           = "#E69F00",
  "Immune"              = "#CC79A7",
  "Skin"                = "#F0E442",
  "Digestive"           = "#1B9E77",
  "Metabolic"           = "#999933",
  "Reproductive female" = "#E78AC3",
  "Reproductive male"   = "#4D4D4D",
  "Other"               = "#999999"
)

# ============================================================
# 4. Helper functions
# ============================================================

infer_modality <- function(clock_folder) {
  x <- tolower(clock_folder)
  
  if (grepl("_mri_mortality_clock$", x)) {
    return("MRI")
  }
  
  if (grepl("_proteomics_mortality_clock$", x)) {
    return("Proteomics")
  }
  
  if (grepl("_metabolomics_mortality_clock$", x)) {
    return("Metabolomics")
  }
  
  return("Unknown")
}

infer_organ_raw <- function(clock_folder) {
  x <- clock_folder
  x <- sub("_mortality_clock$", "", x)
  x <- sub("_mri$", "", x, ignore.case = TRUE)
  x <- sub("_proteomics$", "", x, ignore.case = TRUE)
  x <- sub("_metabolomics$", "", x, ignore.case = TRUE)
  x
}

format_organ_label <- function(organ_raw) {
  x <- as.character(organ_raw)
  x <- gsub("_", " ", x)
  xl <- tolower(x)
  
  if (xl == "brain") return("Brain")
  if (xl == "adipose") return("Adipose")
  if (xl == "heart") return("Heart")
  if (xl %in% c("kidney", "renal")) return("Renal")
  if (xl %in% c("liver", "hepatic")) return("Hepatic")
  if (xl == "pancreas") return("Pancreas")
  if (xl == "spleen") return("Spleen")
  if (xl == "eye") return("Eye")
  if (xl == "pulmonary") return("Pulmonary")
  if (xl == "endocrine") return("Endocrine")
  if (xl == "immune") return("Immune")
  if (xl == "skin") return("Skin")
  if (xl == "digestive") return("Digestive")
  if (xl == "metabolic") return("Metabolic")
  if (xl == "reproductive female") return("Reproductive female")
  if (xl == "reproductive male") return("Reproductive male")
  
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

make_clock_title <- function(organ_label, modality, n_snps) {
  paste0(
    organ_label,
    " ",
    modality,
    "\n",
    "n = ",
    n_snps,
    " lead SNPs"
  )
}

standardize_trumpet_input <- function(dt) {
  dt <- as.data.table(dt)
  
  required_cols <- c("rsID", "freq", "A1_beta", "N", "Gene", "Analysis")
  missing_cols <- setdiff(required_cols, names(dt))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required TrumpetPlots columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dt[, rsID := as.character(rsID)]
  dt[, freq := as.numeric(freq)]
  dt[, A1_beta := as.numeric(A1_beta)]
  dt[, N := as.numeric(N)]
  dt[, Gene := as.character(Gene)]
  dt[, Analysis := as.character(Analysis)]
  
  dt <- dt[
    !is.na(rsID) &
      is.finite(freq) &
      is.finite(A1_beta) &
      is.finite(N) &
      freq > 0 &
      freq < 1
  ]
  
  dt[is.na(Gene) | Gene == "" | Gene == "NaN", Gene := "NA"]
  dt[is.na(Analysis) | Analysis == "", Analysis := "GWAS"]
  
  # Critical fix:
  # Return a real data.table because TrumpetPlots uses data.table :=
  setDT(dt)
  
  dt
}

safe_plot_trumpets <- function(dt, organ_col) {
  # Critical fix:
  # Use copy() because plot_trumpets modifies the data.table by reference.
  dt_for_plot <- copy(as.data.table(dt))
  
  plot_trumpets(
    dataset = dt_for_plot,
    calculate_power = FALSE,
    show_power_curves = FALSE,
    analysis_color_palette = c("GWAS" = organ_col)
  )
}

# ============================================================
# 5. Read manifest
# ============================================================

if (!file.exists(MANIFEST_FILE)) {
  stop("Manifest file does not exist: ", MANIFEST_FILE)
}

manifest <- fread(MANIFEST_FILE)

if (!"clock_folder" %in% names(manifest)) {
  stop("Manifest must contain column: clock_folder")
}

if (!"modality" %in% names(manifest)) {
  manifest[, modality := vapply(clock_folder, infer_modality, character(1))]
}

if (!"organ_raw" %in% names(manifest)) {
  manifest[, organ_raw := vapply(clock_folder, infer_organ_raw, character(1))]
}

if (!"organ_label" %in% names(manifest)) {
  manifest[, organ_label := vapply(organ_raw, format_organ_label, character(1))]
}

manifest[, clock_order := match(clock_folder, EXPECTED_CLOCKS)]
manifest[is.na(clock_order), clock_order := 999]

manifest[, modality_order := match(modality, MODALITY_ORDER)]
manifest[is.na(modality_order), modality_order := 999]

setorder(manifest, clock_order, modality_order, organ_label, clock_folder)

message("Found ", nrow(manifest), " clocks in manifest.")
message("Expected 22 clocks.")

if (nrow(manifest) != 22) {
  warning("Manifest does not contain exactly 22 rows.")
}

# ============================================================
# 6. Generate individual plots
# ============================================================

plot_list <- list()
summary_list <- list()

for (i in seq_len(nrow(manifest))) {
  clock_folder <- manifest$clock_folder[i]
  organ_label <- manifest$organ_label[i]
  modality <- manifest$modality[i]
  
  input_file <- file.path(
    INPUT_DIR,
    paste0("TrumpetPlots_input_", clock_folder, ".tsv")
  )
  
  message("------------------------------------------------------------")
  message("Clock: ", clock_folder)
  message("Input: ", input_file)
  
  if (!file.exists(input_file)) {
    warning("Input file not found for clock: ", clock_folder)
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      clock_folder = clock_folder,
      organ_label = organ_label,
      modality = modality,
      status = "missing_file",
      n_snps = 0,
      input_file = input_file,
      individual_pdf = NA_character_,
      individual_png = NA_character_
    )
    
    next
  }
  
  dt <- fread(input_file)
  dt <- standardize_trumpet_input(dt)
  
  if (nrow(dt) == 0) {
    warning("No valid rows after filtering for clock: ", clock_folder)
    
    summary_list[[length(summary_list) + 1]] <- data.table(
      clock_folder = clock_folder,
      organ_label = organ_label,
      modality = modality,
      status = "no_valid_rows",
      n_snps = 0,
      input_file = input_file,
      individual_pdf = NA_character_,
      individual_png = NA_character_
    )
    
    next
  }
  
  organ_col <- organ_colors[[organ_label]]
  
  if (is.null(organ_col) || is.na(organ_col)) {
    organ_col <- organ_colors[["Other"]]
  }
  
  p <- safe_plot_trumpets(dt, organ_col) +
    ggtitle(make_clock_title(organ_label, modality, nrow(dt))) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9.2),
      axis.title = element_text(face = "bold", size = 8),
      axis.text = element_text(size = 7),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, color = "grey88"),
      plot.margin = margin(3, 3, 3, 3)
    )
  
  plot_list[[clock_folder]] <- p
  
  individual_pdf <- file.path(
    OUTDIR,
    paste0("TrumpetPlot_", clock_folder, ".pdf")
  )
  
  individual_png <- file.path(
    OUTDIR,
    paste0("TrumpetPlot_", clock_folder, ".png")
  )
  
  ggsave(
    filename = individual_pdf,
    plot = p,
    width = 4.0,
    height = 3.4,
    device = "pdf"
  )
  
  ggsave(
    filename = individual_png,
    plot = p,
    width = 4.0,
    height = 3.4,
    dpi = 300
  )
  
  summary_list[[length(summary_list) + 1]] <- data.table(
    clock_folder = clock_folder,
    organ_label = organ_label,
    modality = modality,
    status = "ok",
    n_snps = nrow(dt),
    input_file = input_file,
    individual_pdf = individual_pdf,
    individual_png = individual_png
  )
}

summary_df <- rbindlist(summary_list, fill = TRUE)

summary_out <- file.path(
  OUTDIR,
  "TrumpetPlots_22_mortality_epoch_plot_summary.tsv"
)

fwrite(summary_df, summary_out, sep = "\t")

# ============================================================
# 7. Combined figure: all 22 clocks
# ============================================================

if (length(plot_list) == 0) {
  stop("No trumpet plots were generated.")
}

combined_pdf <- file.path(
  OUTDIR,
  "TrumpetPlots_22_mortality_epoch_all_clocks.pdf"
)

combined_png <- file.path(
  OUTDIR,
  "TrumpetPlots_22_mortality_epoch_all_clocks.png"
)

pdf(
  combined_pdf,
  width = 16.5,
  height = 18.0,
  onefile = TRUE
)

grid.arrange(
  grobs = plot_list,
  ncol = 4,
  top = grid::textGrob(
    "Trumpet plots for 22 mortality EPOCH clocks",
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )
)

dev.off()

png(
  combined_png,
  width = 16.5,
  height = 18.0,
  units = "in",
  res = 300
)

grid.arrange(
  grobs = plot_list,
  ncol = 4,
  top = grid::textGrob(
    "Trumpet plots for 22 mortality EPOCH clocks",
    gp = grid::gpar(fontsize = 18, fontface = "bold")
  )
)

dev.off()

# ============================================================
# 8. Modality-specific combined figures
# ============================================================

for (mod in MODALITY_ORDER) {
  clock_ids <- summary_df[
    status == "ok" & modality == mod,
    clock_folder
  ]
  
  if (length(clock_ids) == 0) {
    next
  }
  
  mod_plots <- plot_list[clock_ids]
  
  mod_pdf <- file.path(
    OUTDIR,
    paste0("TrumpetPlots_mortality_epoch_", mod, "_clocks.pdf")
  )
  
  mod_png <- file.path(
    OUTDIR,
    paste0("TrumpetPlots_mortality_epoch_", mod, "_clocks.png")
  )
  
  ncol_mod <- if (mod == "MRI") {
    4
  } else if (mod == "Proteomics") {
    4
  } else if (mod == "Metabolomics") {
    2
  } else {
    3
  }
  
  nrow_mod <- ceiling(length(mod_plots) / ncol_mod)
  
  pdf(
    mod_pdf,
    width = 4.1 * ncol_mod,
    height = 3.7 * nrow_mod + 0.7,
    onefile = TRUE
  )
  
  grid.arrange(
    grobs = mod_plots,
    ncol = ncol_mod,
    top = grid::textGrob(
      paste0("Trumpet plots for mortality EPOCH: ", mod),
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )
  )
  
  dev.off()
  
  png(
    mod_png,
    width = 4.1 * ncol_mod,
    height = 3.7 * nrow_mod + 0.7,
    units = "in",
    res = 300
  )
  
  grid.arrange(
    grobs = mod_plots,
    ncol = ncol_mod,
    top = grid::textGrob(
      paste0("Trumpet plots for mortality EPOCH: ", mod),
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )
  )
  
  dev.off()
}

# ============================================================
# 9. Optional overlay by modality
#    Fix: do not use facet_wrap(~ modality) on plot_trumpets(),
#    because plot_trumpets() layers do not retain modality.
#    Instead, split by modality and generate one overlay plot per modality.
# ============================================================

if (file.exists(COMBINED_INPUT_FILE)) {
  combined_dt <- fread(COMBINED_INPUT_FILE)
  combined_dt <- standardize_trumpet_input(combined_dt)
  
  if (!"organ_label" %in% names(combined_dt)) {
    if ("organ_raw" %in% names(combined_dt)) {
      combined_dt[, organ_label := vapply(organ_raw, format_organ_label, character(1))]
    } else {
      combined_dt[, organ_label := "Other"]
    }
  }
  
  if (!"modality" %in% names(combined_dt)) {
    if ("clock_folder" %in% names(combined_dt)) {
      combined_dt[, modality := vapply(clock_folder, infer_modality, character(1))]
    } else {
      combined_dt[, modality := "Unknown"]
    }
  }
  
  combined_dt[, modality := as.character(modality)]
  combined_dt[, organ_label := as.character(organ_label)]
  
  overlay_plot_list <- list()
  overlay_summary_list <- list()
  
  for (mod in MODALITY_ORDER) {
    dt_mod <- combined_dt[modality == mod]
    
    if (nrow(dt_mod) == 0) {
      next
    }
    
    # Color points by organ/system within this modality.
    dt_mod[, Analysis := organ_label]
    
    organ_levels_mod <- sort(unique(dt_mod$organ_label))
    organ_palette_mod <- organ_colors[organ_levels_mod]
    organ_palette_mod[is.na(organ_palette_mod)] <- organ_colors[["Other"]]
    
    # Critical: pass a real copied data.table because TrumpetPlots modifies by reference.
    p_mod <- plot_trumpets(
      dataset = copy(as.data.table(dt_mod)),
      calculate_power = FALSE,
      show_power_curves = FALSE,
      analysis_color_palette = organ_palette_mod
    ) +
      ggtitle(paste0(mod, " mortality EPOCH")) +
      theme_bw(base_size = 10.5) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 8),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.18, color = "grey88"),
        plot.margin = margin(4, 4, 4, 4)
      )
    
    overlay_plot_list[[mod]] <- p_mod
    
    overlay_summary_list[[length(overlay_summary_list) + 1]] <- data.table(
      modality = mod,
      n_rows = nrow(dt_mod),
      n_clocks = uniqueN(dt_mod$clock_folder),
      n_organs = uniqueN(dt_mod$organ_label),
      organs = paste(sort(unique(dt_mod$organ_label)), collapse = ";")
    )
    
    mod_overlay_pdf <- file.path(
      OUTDIR,
      paste0("TrumpetPlots_22_mortality_epoch_overlay_", mod, ".pdf")
    )
    
    mod_overlay_png <- file.path(
      OUTDIR,
      paste0("TrumpetPlots_22_mortality_epoch_overlay_", mod, ".png")
    )
    
    ggsave(
      filename = mod_overlay_pdf,
      plot = p_mod,
      width = 5.8,
      height = 4.8,
      device = "pdf"
    )
    
    ggsave(
      filename = mod_overlay_png,
      plot = p_mod,
      width = 5.8,
      height = 4.8,
      dpi = 300
    )
  }
  
  if (length(overlay_summary_list) > 0) {
    overlay_summary_dt <- rbindlist(overlay_summary_list, fill = TRUE)
    
    overlay_summary_out <- file.path(
      OUTDIR,
      "TrumpetPlots_22_mortality_epoch_overlay_by_modality_summary.tsv"
    )
    
    fwrite(overlay_summary_dt, overlay_summary_out, sep = "\t")
  }
  
  if (length(overlay_plot_list) > 0) {
    overlay_pdf <- file.path(
      OUTDIR,
      "TrumpetPlots_22_mortality_epoch_overlay_by_modality.pdf"
    )
    
    overlay_png <- file.path(
      OUTDIR,
      "TrumpetPlots_22_mortality_epoch_overlay_by_modality.png"
    )
    
    pdf(
      overlay_pdf,
      width = 16.5,
      height = 5.8,
      onefile = TRUE
    )
    
    grid.arrange(
      grobs = overlay_plot_list,
      ncol = length(overlay_plot_list),
      top = grid::textGrob(
        "Trumpet plot overlay for 22 mortality EPOCH clocks by modality",
        gp = grid::gpar(fontsize = 16, fontface = "bold")
      )
    )
    
    dev.off()
    
    png(
      overlay_png,
      width = 16.5,
      height = 5.8,
      units = "in",
      res = 300
    )
    
    grid.arrange(
      grobs = overlay_plot_list,
      ncol = length(overlay_plot_list),
      top = grid::textGrob(
        "Trumpet plot overlay for 22 mortality EPOCH clocks by modality",
        gp = grid::gpar(fontsize = 16, fontface = "bold")
      )
    )
    
    dev.off()
  }
}

# ============================================================
# 10. Console output
# ============================================================

message("============================================================")
message("Finished TrumpetPlots for 22 mortality EPOCH clocks.")
message("Output directory:")
message(OUTDIR)
message("")
message("Summary:")
message(summary_out)
message("")
message("Main combined figure:")
message(combined_pdf)
message(combined_png)
message("============================================================")

print("Stop")