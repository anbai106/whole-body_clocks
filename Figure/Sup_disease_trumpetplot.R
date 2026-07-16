#!/usr/bin/env Rscript

# ============================================================
# Compact combined TrumpetPlots for 47 disease EPOCH clocks
# Creates p_combine as a single patchwork object for RStudio
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(TrumpetPlots)
  library(ggplot2)
  library(grid)
  library(patchwork)
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
  "TrumpetPlots_47_disease_epoch"
)

COMBINED_INPUT_FILE <- file.path(
  INPUT_DIR,
  "TrumpetPlots_input_47_disease_epoch_combined.tsv"
)

MANIFEST_FILE <- file.path(
  INPUT_DIR,
  "TrumpetPlots_47_disease_epoch_input_manifest.tsv"
)

OUTDIR <- file.path(
  INPUT_DIR,
  "figures_47_disease_epoch"
)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("Base directory: ", BASE_DIR)
message("Input directory: ", INPUT_DIR)
message("Combined input: ", COMBINED_INPUT_FILE)
message("Output directory: ", OUTDIR)

if (!file.exists(COMBINED_INPUT_FILE)) {
  stop("Combined TrumpetPlots input file does not exist: ", COMBINED_INPUT_FILE)
}

# ============================================================
# 2. Orders and colors
# ============================================================

DISEASE_ORDER <- c("asthma", "copd", "dementia", "mi", "stroke")

DISEASE_LABELS <- c(
  "asthma" = "Asthma",
  "copd" = "COPD",
  "dementia" = "Dementia",
  "mi" = "MI",
  "stroke" = "Stroke"
)

MODALITY_ORDER <- c("MRI", "Proteomics", "Metabolomics")

ORGAN_ORDER <- c(
  "Brain",
  "Eye",
  "Heart",
  "Hepatic",
  "Renal",
  "Pulmonary",
  "Endocrine",
  "Immune",
  "Skin",
  "Digestive",
  "Metabolic",
  "Adipose",
  "Pancreas",
  "Spleen",
  "Reproductive female",
  "Reproductive male",
  "Other"
)

organ_colors <- c(
  "Brain"               = "#0072B2",
  "Eye"                 = "#56B4E9",
  "Heart"               = "#D55E00",
  "Hepatic"             = "#009E73",
  "Liver"               = "#009E73",
  "Renal"               = "#8B5A2B",
  "Kidney"              = "#8B5A2B",
  "Pulmonary"           = "#7B61A8",
  "Endocrine"           = "#E69F00",
  "Immune"              = "#CC79A7",
  "Skin"                = "#F0E442",
  "Digestive"           = "#1B9E77",
  "Metabolic"           = "#999933",
  "Adipose"             = "#A6761D",
  "Pancreas"            = "#66A61E",
  "Spleen"              = "#666666",
  "Reproductive female" = "#E78AC3",
  "Reproductive male"   = "#4D4D4D",
  "Other"               = "#999999"
)

disease_colors <- c(
  "Asthma"   = "#0072B2",
  "COPD"     = "#7B61A8",
  "Dementia" = "#009E73",
  "MI"       = "#D55E00",
  "Stroke"   = "#CC79A7"
)

# ============================================================
# 3. Helper functions
# ============================================================

clean_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

canonical_organ_label <- function(x) {
  x <- clean_chr(x)
  xl <- tolower(gsub("_", " ", x))

  if (xl == "brain") return("Brain")
  if (xl == "eye") return("Eye")
  if (xl == "heart") return("Heart")
  if (xl %in% c("hepatic", "liver")) return("Hepatic")
  if (xl %in% c("renal", "kidney")) return("Renal")
  if (xl %in% c("pulmonary", "lung")) return("Pulmonary")
  if (xl == "endocrine") return("Endocrine")
  if (xl == "immune") return("Immune")
  if (xl == "skin") return("Skin")
  if (xl == "digestive") return("Digestive")
  if (xl == "metabolic") return("Metabolic")
  if (xl == "adipose") return("Adipose")
  if (xl == "pancreas") return("Pancreas")
  if (xl == "spleen") return("Spleen")
  if (xl %in% c("reproductive female", "female reproductive")) return("Reproductive female")
  if (xl %in% c("reproductive male", "male reproductive")) return("Reproductive male")

  if (xl == "") return("Other")

  paste0(toupper(substr(xl, 1, 1)), substr(xl, 2, nchar(xl)))
}

parse_clock_folder <- function(clock_folder) {
  x <- as.character(clock_folder)

  m <- regexec(
    pattern = "^(.*)_(mri|proteomics|metabolomics)_(asthma|copd|dementia|mi|stroke)_clock$",
    text = x,
    ignore.case = TRUE
  )

  z <- regmatches(x, m)[[1]]

  if (length(z) != 4) {
    return(data.table(
      clock_folder = x,
      parse_ok = FALSE,
      organ_raw_parsed = "",
      organ_label_parsed = "Other",
      modality_parsed = "Unknown",
      disease_parsed = "",
      disease_label_parsed = "Unknown"
    ))
  }

  organ_raw <- z[2]
  modality_key <- tolower(z[3])
  disease_key <- tolower(z[4])

  modality_label <- ifelse(
    modality_key == "mri",
    "MRI",
    ifelse(
      modality_key == "proteomics",
      "Proteomics",
      ifelse(modality_key == "metabolomics", "Metabolomics", "Unknown")
    )
  )

  data.table(
    clock_folder = x,
    parse_ok = TRUE,
    organ_raw_parsed = organ_raw,
    organ_label_parsed = canonical_organ_label(organ_raw),
    modality_parsed = modality_label,
    disease_parsed = disease_key,
    disease_label_parsed = DISEASE_LABELS[[disease_key]]
  )
}

make_clock_title <- function(disease_label, organ_label, modality, n_snps) {
  modality_short <- ifelse(
    modality == "Proteomics",
    "Prot.",
    ifelse(modality == "Metabolomics", "Metab.", modality)
  )

  organ_short <- organ_label
  organ_short <- gsub("Reproductive female", "Repro. F", organ_short)
  organ_short <- gsub("Reproductive male", "Repro. M", organ_short)

  paste0(
    disease_label,
    " | ",
    organ_short,
    " ",
    modality_short,
    "\n",
    "n = ",
    n_snps
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

  setDT(dt)

  dt
}

safe_plot_trumpets <- function(dt, organ_col) {
  dt_for_plot <- copy(as.data.table(dt))

  plot_trumpets(
    dataset = dt_for_plot,
    calculate_power = FALSE,
    show_power_curves = FALSE,
    analysis_color_palette = c("GWAS" = organ_col)
  )
}

make_empty_plot <- function(disease_label, organ_label, modality) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No lead SNPs",
      size = 3,
      fontface = "bold",
      color = "grey35"
    ) +
    ggtitle(make_clock_title(disease_label, organ_label, modality, 0)) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 6.8,
        lineheight = 0.88,
        margin = margin(b = 1)
      ),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5)
    )
}

compact_trumpet_theme <- function() {
  theme_bw(base_size = 7.2) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 6.8,
        lineheight = 0.88,
        margin = margin(b = 1)
      ),
      axis.title = element_text(face = "bold", size = 6.1),
      axis.text = element_text(size = 5.3),
      axis.ticks = element_line(linewidth = 0.16),
      panel.grid.major = element_line(linewidth = 0.10, color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1.4, 1.4, 1.4, 1.4),
      legend.position = "none"
    )
}

# ============================================================
# 4. Read combined data
# ============================================================

combined_dt <- fread(COMBINED_INPUT_FILE)
combined_dt <- standardize_trumpet_input(combined_dt)

if (!"clock_folder" %in% names(combined_dt)) {
  stop("Combined input must contain column: clock_folder")
}

parsed_from_combined <- rbindlist(
  lapply(unique(combined_dt$clock_folder), parse_clock_folder),
  fill = TRUE
)

# Add or repair metadata columns in the combined data.
combined_dt <- merge(
  combined_dt,
  parsed_from_combined,
  by = "clock_folder",
  all.x = TRUE
)

if (!"disease" %in% names(combined_dt)) {
  combined_dt[, disease := disease_parsed]
} else {
  combined_dt[is.na(disease) | disease == "", disease := disease_parsed]
}

if (!"disease_label" %in% names(combined_dt)) {
  combined_dt[, disease_label := disease_label_parsed]
} else {
  combined_dt[is.na(disease_label) | disease_label == "", disease_label := disease_label_parsed]
}

if (!"modality" %in% names(combined_dt)) {
  combined_dt[, modality := modality_parsed]
} else {
  combined_dt[is.na(modality) | modality == "", modality := modality_parsed]
}

if (!"organ_label" %in% names(combined_dt)) {
  combined_dt[, organ_label := organ_label_parsed]
} else {
  combined_dt[is.na(organ_label) | organ_label == "", organ_label := organ_label_parsed]
}

if (!"organ_raw" %in% names(combined_dt)) {
  combined_dt[, organ_raw := organ_raw_parsed]
} else {
  combined_dt[is.na(organ_raw) | organ_raw == "", organ_raw := organ_raw_parsed]
}

combined_dt[, organ_label := vapply(organ_label, canonical_organ_label, character(1))]
combined_dt[, modality := as.character(modality)]
combined_dt[, disease := tolower(as.character(disease))]
combined_dt[, disease_label := as.character(disease_label)]
combined_dt[is.na(disease_label) | disease_label == "", disease_label := DISEASE_LABELS[disease]]

message("Rows in combined input after QC: ", nrow(combined_dt))
message("Unique clocks in combined input: ", uniqueN(combined_dt$clock_folder))

# ============================================================
# 5. Build plot manifest
#    Prefer Python-generated 47-clock manifest if available,
#    so clocks with no rows can still be represented.
# ============================================================

if (file.exists(MANIFEST_FILE)) {
  manifest <- fread(MANIFEST_FILE)

  if (!"clock_folder" %in% names(manifest)) {
    stop("Manifest exists but does not contain column: clock_folder")
  }

  if ("parse_ok" %in% names(manifest)) {
    manifest <- manifest[
      parse_ok %in% c(TRUE, "TRUE", "True", "true", 1, "1")
    ]
  }

  parsed_from_manifest <- rbindlist(
    lapply(unique(manifest$clock_folder), parse_clock_folder),
    fill = TRUE
  )

  manifest <- merge(
    manifest,
    parsed_from_manifest,
    by = "clock_folder",
    all.x = TRUE,
    suffixes = c("", "_parsed2")
  )

  if (!"disease" %in% names(manifest)) {
    manifest[, disease := disease_parsed]
  } else {
    manifest[is.na(disease) | disease == "", disease := disease_parsed]
  }

  if (!"disease_label" %in% names(manifest)) {
    manifest[, disease_label := disease_label_parsed]
  } else {
    manifest[is.na(disease_label) | disease_label == "", disease_label := disease_label_parsed]
  }

  if (!"modality" %in% names(manifest)) {
    manifest[, modality := modality_parsed]
  } else {
    manifest[is.na(modality) | modality == "", modality := modality_parsed]
  }

  if (!"organ_label" %in% names(manifest)) {
    manifest[, organ_label := organ_label_parsed]
  } else {
    manifest[is.na(organ_label) | organ_label == "", organ_label := organ_label_parsed]
  }

  if (!"organ_raw" %in% names(manifest)) {
    manifest[, organ_raw := organ_raw_parsed]
  } else {
    manifest[is.na(organ_raw) | organ_raw == "", organ_raw := organ_raw_parsed]
  }

  plot_manifest <- unique(
    manifest[
      ,
      .(
        clock_folder,
        disease,
        disease_label,
        organ_raw,
        organ_label,
        modality
      )
    ]
  )

} else {
  plot_manifest <- unique(
    combined_dt[
      ,
      .(
        clock_folder,
        disease,
        disease_label,
        organ_raw,
        organ_label,
        modality
      )
    ]
  )
}

plot_manifest[, disease := tolower(as.character(disease))]
plot_manifest[, disease_label := as.character(disease_label)]
plot_manifest[is.na(disease_label) | disease_label == "", disease_label := DISEASE_LABELS[disease]]

plot_manifest[, organ_label := vapply(organ_label, canonical_organ_label, character(1))]
plot_manifest[, modality := as.character(modality)]

plot_manifest[, disease_order := match(disease, DISEASE_ORDER)]
plot_manifest[is.na(disease_order), disease_order := 999]

plot_manifest[, modality_order := match(modality, MODALITY_ORDER)]
plot_manifest[is.na(modality_order), modality_order := 999]

plot_manifest[, organ_order := match(organ_label, ORGAN_ORDER)]
plot_manifest[is.na(organ_order), organ_order := 999]

setorder(
  plot_manifest,
  disease_order,
  modality_order,
  organ_order,
  organ_label,
  clock_folder
)

message("Clocks in plot manifest: ", nrow(plot_manifest))

if (nrow(plot_manifest) != 47) {
  warning(
    "Plot manifest does not contain exactly 47 clocks. It contains ",
    nrow(plot_manifest),
    " clocks. The script will continue."
  )
}

# ============================================================
# 6. Generate one trumpet plot per disease EPOCH clock
# ============================================================

plot_list <- list()
summary_list <- list()

for (i in seq_len(nrow(plot_manifest))) {
  
  clock_id <- plot_manifest$clock_folder[i]
  disease_label_i <- plot_manifest$disease_label[i]
  disease_i <- plot_manifest$disease[i]
  organ_label_i <- plot_manifest$organ_label[i]
  modality_i <- plot_manifest$modality[i]
  
  message("------------------------------------------------------------")
  message("Clock:    ", clock_id)
  message("Disease:  ", disease_label_i)
  message("Organ:    ", organ_label_i)
  message("Modality: ", modality_i)
  
  # IMPORTANT:
  # Do not use !! inside data.table.
  # Use an external variable with a different name from the column.
  dt_clock <- combined_dt[clock_folder == clock_id]
  
  if (nrow(dt_clock) > 0) {
    dt_clock <- standardize_trumpet_input(dt_clock)
  }
  
  if (nrow(dt_clock) == 0) {
    
    p <- make_empty_plot(
      disease_label = disease_label_i,
      organ_label = organ_label_i,
      modality = modality_i
    )
    
    status <- "no_valid_rows"
    n_snps <- 0
    
  } else {
    
    organ_col <- organ_colors[[organ_label_i]]
    
    if (is.null(organ_col) || is.na(organ_col)) {
      organ_col <- organ_colors[["Other"]]
    }
    
    dt_clock[, Analysis := "GWAS"]
    
    p <- safe_plot_trumpets(dt_clock, organ_col) +
      ggtitle(
        make_clock_title(
          disease_label = disease_label_i,
          organ_label = organ_label_i,
          modality = modality_i,
          n_snps = nrow(dt_clock)
        )
      ) +
      compact_trumpet_theme()
    
    status <- "ok"
    n_snps <- nrow(dt_clock)
  }
  
  plot_list[[clock_id]] <- p
  
  individual_pdf <- file.path(
    OUTDIR,
    paste0("TrumpetPlot_", clock_id, ".pdf")
  )
  
  individual_png <- file.path(
    OUTDIR,
    paste0("TrumpetPlot_", clock_id, ".png")
  )
  
  ggsave(
    filename = individual_pdf,
    plot = p,
    width = 3.2,
    height = 2.75,
    device = "pdf"
  )
  
  ggsave(
    filename = individual_png,
    plot = p,
    width = 3.2,
    height = 2.75,
    dpi = 300
  )
  
  summary_list[[length(summary_list) + 1]] <- data.table(
    clock_folder = clock_id,
    disease = disease_i,
    disease_label = disease_label_i,
    organ_label = organ_label_i,
    modality = modality_i,
    status = status,
    n_snps = n_snps,
    individual_pdf = individual_pdf,
    individual_png = individual_png
  )
}

summary_df <- rbindlist(summary_list, fill = TRUE)

summary_out <- file.path(
  OUTDIR,
  "TrumpetPlots_47_disease_epoch_plot_summary.tsv"
)

fwrite(summary_df, summary_out, sep = "\t")

# ============================================================
# 7. Compact combined figure: all 47 disease EPOCH clocks
#    This creates p_combine as the main RStudio object.
# ============================================================

plot_ids_ordered <- summary_df$clock_folder
plot_list_ordered <- plot_list[plot_ids_ordered]

# Use 6 columns to keep all 47 clocks in one compact figure.
N_COL_COMBINE <- 6

p_combine <- patchwork::wrap_plots(
  plot_list_ordered,
  ncol = N_COL_COMBINE,
  byrow = TRUE
) +
  patchwork::plot_annotation(
    title = "Trumpet plots for 47 disease EPOCH clocks",
    subtitle = "Lead SNPs from FUMA independent significant SNPs merged with fastGWA effect size, allele frequency, and sample size",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 17,
        margin = margin(b = 3)
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 9.5,
        color = "grey30",
        margin = margin(b = 6)
      )
    )
  )

# Print/view in RStudio.
p_combine

combined_pdf <- file.path(
  OUTDIR,
  "TrumpetPlots_47_disease_epoch_all_clocks_compact_p_combine.pdf"
)

combined_png <- file.path(
  OUTDIR,
  "TrumpetPlots_47_disease_epoch_all_clocks_compact_p_combine.png"
)

combined_rds <- file.path(
  OUTDIR,
  "TrumpetPlots_47_disease_epoch_all_clocks_compact_p_combine.rds"
)

ggsave(
  filename = combined_pdf,
  plot = p_combine,
  width = 18.0,
  height = 21.5,
  device = "pdf"
)

ggsave(
  filename = combined_png,
  plot = p_combine,
  width = 18.0,
  height = 21.5,
  dpi = 350
)

saveRDS(
  p_combine,
  file = combined_rds
)

# ============================================================
# 8. Disease-specific compact figures
# ============================================================

p_combine_by_disease <- list()

for (dis in DISEASE_ORDER) {
  ids_dis <- summary_df[disease == dis, clock_folder]

  if (length(ids_dis) == 0) {
    next
  }

  plots_dis <- plot_list[ids_dis]
  disease_label <- DISEASE_LABELS[[dis]]

  p_dis <- patchwork::wrap_plots(
    plots_dis,
    ncol = 4,
    byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = paste0("Trumpet plots for ", disease_label, " EPOCH clocks"),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 15,
          margin = margin(b = 6),
          color = disease_colors[[disease_label]]
        )
      )
    )

  p_combine_by_disease[[disease_label]] <- p_dis

  disease_pdf <- file.path(
    OUTDIR,
    paste0("TrumpetPlots_47_disease_epoch_", dis, "_compact.pdf")
  )

  disease_png <- file.path(
    OUTDIR,
    paste0("TrumpetPlots_47_disease_epoch_", dis, "_compact.png")
  )

  nrow_dis <- ceiling(length(plots_dis) / 4)
  height_dis <- max(4.8, 2.75 * nrow_dis + 0.75)

  ggsave(
    filename = disease_pdf,
    plot = p_dis,
    width = 12.8,
    height = height_dis,
    device = "pdf"
  )

  ggsave(
    filename = disease_png,
    plot = p_dis,
    width = 12.8,
    height = height_dis,
    dpi = 350
  )
}

disease_rds <- file.path(
  OUTDIR,
  "TrumpetPlots_47_disease_epoch_p_combine_by_disease.rds"
)

saveRDS(
  p_combine_by_disease,
  file = disease_rds
)

# ============================================================
# 9. Modality-specific compact figures
# ============================================================

p_combine_by_modality <- list()

for (mod in MODALITY_ORDER) {
  ids_mod <- summary_df[modality == mod, clock_folder]

  if (length(ids_mod) == 0) {
    next
  }

  plots_mod <- plot_list[ids_mod]

  ncol_mod <- if (mod == "MRI") {
    2
  } else if (mod == "Proteomics") {
    5
  } else if (mod == "Metabolomics") {
    4
  } else {
    4
  }

  p_mod <- patchwork::wrap_plots(
    plots_mod,
    ncol = ncol_mod,
    byrow = TRUE
  ) +
    patchwork::plot_annotation(
      title = paste0("Trumpet plots for disease EPOCH clocks: ", mod),
      theme = theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 15,
          margin = margin(b = 6)
        )
      )
    )

  p_combine_by_modality[[mod]] <- p_mod

  mod_pdf <- file.path(
    OUTDIR,
    paste0("TrumpetPlots_47_disease_epoch_", mod, "_compact.pdf")
  )

  mod_png <- file.path(
    OUTDIR,
    paste0("TrumpetPlots_47_disease_epoch_", mod, "_compact.png")
  )

  nrow_mod <- ceiling(length(plots_mod) / ncol_mod)
  height_mod <- max(4.8, 2.75 * nrow_mod + 0.75)

  ggsave(
    filename = mod_pdf,
    plot = p_mod,
    width = 3.25 * ncol_mod,
    height = height_mod,
    device = "pdf"
  )

  ggsave(
    filename = mod_png,
    plot = p_mod,
    width = 3.25 * ncol_mod,
    height = height_mod,
    dpi = 350
  )
}

modality_rds <- file.path(
  OUTDIR,
  "TrumpetPlots_47_disease_epoch_p_combine_by_modality.rds"
)

saveRDS(
  p_combine_by_modality,
  file = modality_rds
)

# ============================================================
# 10. Console output
# ============================================================

message("============================================================")
message("Finished TrumpetPlots for 47 disease EPOCH clocks.")
message("Output directory:")
message(OUTDIR)
message("")
message("Summary:")
message(summary_out)
message("")
message("Main p_combine object saved:")
message(combined_rds)
message("")
message("Main combined figure:")
message(combined_pdf)
message(combined_png)
message("")
message("Disease-specific p_combine list:")
message(disease_rds)
message("")
message("Modality-specific p_combine list:")
message(modality_rds)
message("============================================================")

print("Stop")