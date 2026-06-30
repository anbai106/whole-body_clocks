#!/usr/bin/env Rscript

# ============================================================
# Multi-organ visualization:
# Mortality L'EPOCH clocks vs first-generation biological aging clocks
#
# Input:
#   clock_vs_BAG_powered_ok_Ncase_ge_50.tsv
#
# Logic requested:
#   1) Bonferroni threshold = 0.05 / n_diseases / n_clock_pairs
#   2) Keep only disease-clock-pair rows where either:
#        mortality L'EPOCH P < threshold OR BAG P < threshold
#   3) For every kept row, plot both L'EPOCH and BAG HRs for direct comparison
#   4) Annotate statistically significant effect-size difference:
#        joint_p_diff < threshold
#        joint_beta_diff_clock_minus_bag > 0 => L'EPOCH stronger
#        joint_beta_diff_clock_minus_bag < 0 => aging clock stronger
#
# Important plotting fix:
#   With ~1000 disease endpoints, a single 400-dpi PNG exceeds ragg's
#   50,000-pixel dimension limit. This script therefore saves readable
#   paginated figures and body-system figures by default.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# =========================
# 1. Paths
# =========================

stat_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/clock_vs_BAG_survival_summary/clock_vs_BAG_powered_ok_Ncase_ge_50.tsv"

# Optional disease annotation file.
# If available, create a TSV with columns:
# disease_id    disease_name    body_system    body_system_order
# If unavailable, the script uses ICD code as disease_name and orders by ICD chapter.
disease_map_file <- "/Users/hao/Downloads/disease_annotation.tsv"

out_dir <- "/Users/hao/Downloads/LEPOCH_vs_BAG_multiorgan_significant_only"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 2. Figure options
# =========================

min_case <- 50
alpha <- 0.05

# Show only disease-clock-pair rows where L'EPOCH or BAG is Bonferroni-significant.
show_only_significant_pairs <- TRUE

# Main fix for the ragg error:
# Do not save a single giant all-system PNG. Use paginated pages instead.
save_full_all_systems_pdf <- FALSE
save_full_all_systems_png <- FALSE

# Recommended outputs.
save_paginated_all_systems <- TRUE
save_by_body_system <- TRUE

# Number of disease rows per paginated figure.
# 40-60 is usually readable for manuscript/supplement panels.
max_rows_per_page <- 45

# PNG settings. The script automatically lowers DPI if needed to keep below 50,000 pixels.
png_dpi_target <- 300
max_png_pixels <- 49000

# Set to Inf to keep all retained disease endpoints. Use a finite number for a compact main figure.
max_diseases_to_plot <- Inf

# =========================
# 3. Helper functions
# =========================

to_num <- function(x) suppressWarnings(as.numeric(x))

winsorize <- function(x, probs = c(0.01, 0.99)) {
  x <- as.numeric(x)
  qs <- quantile(x, probs = probs, na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}

format_range <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("NA")
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  if (abs(xmin - xmax) < 1e-8) {
    return(format(round(xmin, digits), nsmall = digits))
  }
  paste0(
    format(round(xmin, digits), nsmall = digits),
    "-",
    format(round(xmax, digits), nsmall = digits)
  )
}

icd_body_system <- function(id) {
  id <- as.character(id)
  first_letter <- str_sub(id, 1, 1)

  case_when(
    first_letter == "G" ~ "Brain/neurologic",
    first_letter == "F" ~ "Mental/behavioral",
    first_letter == "H" ~ "Eye/ear",
    first_letter == "I" ~ "Cardiovascular",
    first_letter == "J" ~ "Respiratory",
    first_letter == "K" ~ "Digestive/hepatic",
    first_letter == "E" ~ "Endocrine/metabolic",
    first_letter == "C" ~ "Cancer",
    first_letter == "D" ~ "Blood/immune",
    first_letter == "N" ~ "Renal/urogenital",
    first_letter == "L" ~ "Skin",
    first_letter == "M" ~ "Musculoskeletal",
    first_letter %in% c("S", "T") ~ "Injury/trauma",
    first_letter == "R" ~ "Symptoms/signs",
    first_letter == "Z" ~ "Health factors",
    TRUE ~ "Other"
  )
}

icd_body_order <- function(id) {
  system <- icd_body_system(id)

  case_when(
    system == "Brain/neurologic" ~ 1,
    system == "Mental/behavioral" ~ 2,
    system == "Eye/ear" ~ 3,
    system == "Cardiovascular" ~ 4,
    system == "Respiratory" ~ 5,
    system == "Digestive/hepatic" ~ 6,
    system == "Endocrine/metabolic" ~ 7,
    system == "Cancer" ~ 8,
    system == "Blood/immune" ~ 9,
    system == "Renal/urogenital" ~ 10,
    system == "Skin" ~ 11,
    system == "Musculoskeletal" ~ 12,
    system == "Injury/trauma" ~ 13,
    system == "Symptoms/signs" ~ 14,
    system == "Health factors" ~ 15,
    TRUE ~ 99
  )
}

make_safe_filename <- function(x) {
  x %>%
    str_replace_all("[/ ]+", "_") %>%
    str_replace_all("[^A-Za-z0-9_\\-]", "")
}

save_pdf_png <- function(plot, file_prefix, width, height, save_pdf = TRUE, save_png = TRUE) {
  if (save_pdf) {
    ggsave(
      filename = paste0(file_prefix, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      limitsize = FALSE,
      device = cairo_pdf
    )
  }

  if (save_png) {
    dpi_allowed <- floor(max_png_pixels / max(width, height))
    dpi_use <- min(png_dpi_target, dpi_allowed)

    if (dpi_use < 72) {
      message("Skipping PNG for ", basename(file_prefix),
              " because required DPI would be too low: ", dpi_use)
    } else {
      ggsave(
        filename = paste0(file_prefix, ".png"),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = dpi_use,
        limitsize = FALSE
      )
    }
  }
}

# =========================
# 4. Load data
# =========================

df <- fread(stat_file)

needed_cols <- c(
  "disease_id", "pair_id", "organ", "modality",
  "N", "N_case", "N_noncase",
  "followup_years_min", "followup_years_max",
  "event_followup_years_min", "event_followup_years_max",
  "clock_hr", "clock_ci_lo", "clock_ci_hi", "clock_p",
  "bag_hr", "bag_ci_lo", "bag_ci_hi", "bag_p",
  "clock_joint_beta", "bag_joint_beta",
  "joint_beta_diff_clock_minus_bag", "joint_p_diff",
  "delta_cindex_clock_minus_bag"
)

missing_cols <- setdiff(needed_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input file: ", paste(missing_cols, collapse = ", "))
}

num_cols <- setdiff(needed_cols, c("disease_id", "pair_id", "organ", "modality"))
for (cc in num_cols) df[[cc]] <- to_num(df[[cc]])

df <- df %>%
  filter(!is.na(disease_id), !is.na(pair_id)) %>%
  filter(N_case >= min_case) %>%
  mutate(
    disease_id = as.character(disease_id),
    pair_id = as.character(pair_id),
    organ = as.character(organ),
    modality = as.character(modality)
  )

n_diseases <- n_distinct(df$disease_id)
n_clock_pairs <- n_distinct(df$pair_id)
p_bonf <- alpha / n_diseases / n_clock_pairs

message("n_diseases = ", n_diseases)
message("n_clock_pairs = ", n_clock_pairs)
message("Bonferroni P threshold = ", signif(p_bonf, 4))

# =========================
# 5. Define significance and keep rows
# =========================

df <- df %>%
  mutate(
    clock_sig_bonf = !is.na(clock_p) & clock_p < p_bonf,
    bag_sig_bonf = !is.na(bag_p) & bag_p < p_bonf,
    either_clock_sig_bonf = clock_sig_bonf | bag_sig_bonf,

    diff_sig_bonf = !is.na(joint_p_diff) & joint_p_diff < p_bonf,
    diff_direction = case_when(
      diff_sig_bonf & joint_beta_diff_clock_minus_bag > 0 ~ "L'EPOCH stronger",
      diff_sig_bonf & joint_beta_diff_clock_minus_bag < 0 ~ "Aging clock stronger",
      TRUE ~ "No significant difference"
    ),
    diff_label = case_when(
      diff_sig_bonf & joint_beta_diff_clock_minus_bag > 0 ~ paste0("Delta beta +", sprintf("%.2f", joint_beta_diff_clock_minus_bag)),
      diff_sig_bonf & joint_beta_diff_clock_minus_bag < 0 ~ paste0("Delta beta ", sprintf("%.2f", joint_beta_diff_clock_minus_bag)),
      TRUE ~ ""
    ),
    diff_symbol = case_when(
      diff_sig_bonf & joint_beta_diff_clock_minus_bag > 0 ~ "Delta+",
      diff_sig_bonf & joint_beta_diff_clock_minus_bag < 0 ~ "Delta-",
      TRUE ~ ""
    )
  )

if (show_only_significant_pairs) {
  plot_base <- df %>% filter(either_clock_sig_bonf)
} else {
  selected_diseases <- df %>%
    group_by(disease_id) %>%
    summarise(any_sig = any(either_clock_sig_bonf, na.rm = TRUE), .groups = "drop") %>%
    filter(any_sig) %>%
    pull(disease_id)

  plot_base <- df %>% filter(disease_id %in% selected_diseases)
}

if (nrow(plot_base) == 0) {
  stop("No rows passed the Bonferroni filter.")
}

message("Significant disease-clock-pair rows retained = ", nrow(plot_base))
message("Disease endpoints retained = ", n_distinct(plot_base$disease_id))

# =========================
# 6. Add disease annotation
# =========================

if (file.exists(disease_map_file)) {
  disease_map <- fread(disease_map_file)
  required_map_cols <- c("disease_id", "disease_name", "body_system", "body_system_order")
  missing_map_cols <- setdiff(required_map_cols, names(disease_map))
  if (length(missing_map_cols) > 0) {
    stop("Disease annotation file is missing: ", paste(missing_map_cols, collapse = ", "))
  }

  disease_map <- disease_map %>%
    mutate(
      disease_id = as.character(disease_id),
      disease_name = as.character(disease_name),
      body_system = as.character(body_system),
      body_system_order = to_num(body_system_order)
    )
} else {
  disease_map <- data.frame(disease_id = unique(df$disease_id)) %>%
    mutate(
      disease_name = disease_id,
      body_system = icd_body_system(disease_id),
      body_system_order = icd_body_order(disease_id)
    )
}

plot_base <- plot_base %>%
  left_join(disease_map, by = "disease_id") %>%
  mutate(
    disease_name = ifelse(is.na(disease_name), disease_id, disease_name),
    body_system = ifelse(is.na(body_system), icd_body_system(disease_id), body_system),
    body_system_order = ifelse(is.na(body_system_order), icd_body_order(disease_id), body_system_order),
    disease_label = ifelse(disease_name == disease_id, disease_id, paste0(disease_name, " (", disease_id, ")"))
  )

# =========================
# 7. Disease and pair ordering
# =========================

disease_rank <- plot_base %>%
  group_by(disease_id, disease_label, body_system, body_system_order) %>%
  summarise(
    min_p = min(c(clock_p, bag_p), na.rm = TRUE),
    min_diff_p = min(joint_p_diff, na.rm = TRUE),
    n_sig_pairs = n(),
    n_sig_diff = sum(diff_sig_bonf, na.rm = TRUE),
    max_abs_log_hr = max(abs(log(c(clock_hr, bag_hr))), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(body_system_order, min_p, desc(n_sig_pairs), desc(max_abs_log_hr))

if (is.finite(max_diseases_to_plot)) {
  disease_rank <- disease_rank %>% slice_head(n = max_diseases_to_plot)
  plot_base <- plot_base %>% filter(disease_id %in% disease_rank$disease_id)
}

organ_order <- c(
  "Brain", "Eye", "Heart", "Pulmonary",
  "Liver", "Hepatic", "Digestive",
  "Kidney", "Renal",
  "Endocrine", "Metabolic", "Immune",
  "Spleen", "Pancreas", "Adipose",
  "Skin", "Reproductive_female", "Reproductive_male"
)

modality_order <- c("MRI", "Proteomics", "Metabolomics")

pair_rank <- plot_base %>%
  distinct(pair_id, organ, modality) %>%
  mutate(
    organ_order_id = match(organ, organ_order),
    organ_order_id = ifelse(is.na(organ_order_id), 999, organ_order_id),
    modality_order_id = match(modality, modality_order),
    modality_order_id = ifelse(is.na(modality_order_id), 999, modality_order_id),
    pair_label = paste0(organ, "\n", modality)
  ) %>%
  arrange(modality_order_id, organ_order_id, pair_id)

# =========================
# 8. Build long plot tables
# =========================

lepoch_long <- plot_base %>%
  transmute(
    disease_id, disease_label, body_system, body_system_order,
    pair_id, organ, modality,
    biomarker = "Mortality L'EPOCH",
    hr = clock_hr,
    ci_lo = clock_ci_lo,
    ci_hi = clock_ci_hi,
    p_value = clock_p,
    sig_bonf = clock_sig_bonf,
    diff_sig_bonf,
    diff_direction,
    diff_symbol,
    diff_label,
    joint_beta_diff_clock_minus_bag,
    joint_p_diff,
    N, N_case, N_noncase,
    followup_years_min, followup_years_max,
    event_followup_years_min, event_followup_years_max
  )

bag_long <- plot_base %>%
  transmute(
    disease_id, disease_label, body_system, body_system_order,
    pair_id, organ, modality,
    biomarker = "Biological aging clock",
    hr = bag_hr,
    ci_lo = bag_ci_lo,
    ci_hi = bag_ci_hi,
    p_value = bag_p,
    sig_bonf = bag_sig_bonf,
    diff_sig_bonf,
    diff_direction,
    diff_symbol = "",
    diff_label = "",
    joint_beta_diff_clock_minus_bag,
    joint_p_diff,
    N, N_case, N_noncase,
    followup_years_min, followup_years_max,
    event_followup_years_min, event_followup_years_max
  )

base_disease_levels <- disease_rank$disease_label

long_df <- bind_rows(lepoch_long, bag_long) %>%
  left_join(pair_rank %>% select(pair_id, pair_label), by = "pair_id") %>%
  mutate(
    disease_label = factor(disease_label, levels = rev(base_disease_levels)),
    pair_label = factor(pair_label, levels = pair_rank$pair_label),
    biomarker = factor(biomarker, levels = c("Mortality L'EPOCH", "Biological aging clock")),
    log_hr = log(hr),
    log_hr_plot = winsorize(log_hr),
    hr_label = case_when(
      is.na(hr) ~ "",
      sig_bonf ~ paste0(sprintf("%.2f", hr), "*"),
      TRUE ~ sprintf("%.2f", hr)
    ),
    diff_tile_label = case_when(
      biomarker == "Mortality L'EPOCH" & diff_sig_bonf & joint_beta_diff_clock_minus_bag > 0 ~ "Delta+",
      biomarker == "Mortality L'EPOCH" & diff_sig_bonf & joint_beta_diff_clock_minus_bag < 0 ~ "Delta-",
      TRUE ~ ""
    )
  )

anno_df <- plot_base %>%
  group_by(disease_id, disease_label, body_system, body_system_order) %>%
  summarise(
    N_min = min(N, na.rm = TRUE),
    N_max = max(N, na.rm = TRUE),
    case_min = min(N_case, na.rm = TRUE),
    case_max = max(N_case, na.rm = TRUE),
    fu_min = min(followup_years_min, na.rm = TRUE),
    fu_max = max(followup_years_max, na.rm = TRUE),
    event_fu_min = min(event_followup_years_min, na.rm = TRUE),
    event_fu_max = max(event_followup_years_max, na.rm = TRUE),
    n_sig_pairs = n(),
    n_lepoch_sig = sum(clock_sig_bonf, na.rm = TRUE),
    n_bag_sig = sum(bag_sig_bonf, na.rm = TRUE),
    n_diff_sig = sum(diff_sig_bonf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    disease_label = factor(disease_label, levels = rev(base_disease_levels)),
    N_label = ifelse(N_min == N_max, paste0("N=", comma(N_min)), paste0("N=", comma(N_min), "-", comma(N_max))),
    case_label = ifelse(case_min == case_max, paste0("cases=", comma(case_min)), paste0("cases=", comma(case_min), "-", comma(case_max))),
    followup_label = paste0("FU=", format_range(c(fu_min, fu_max), digits = 1), " y"),
    sig_label = paste0("sig pairs=", n_sig_pairs, "; L=", n_lepoch_sig, "; BAG=", n_bag_sig, "; diff=", n_diff_sig),
    annotation = paste(N_label, case_label, followup_label, sig_label, sep = " | ")
  )

# =========================
# 9. Plotting functions
# =========================

reset_disease_levels <- function(plot_long, plot_anno, disease_labels_ordered) {
  disease_labels_ordered <- as.character(disease_labels_ordered)

  plot_long2 <- plot_long %>%
    filter(as.character(disease_label) %in% disease_labels_ordered) %>%
    mutate(disease_label = factor(as.character(disease_label), levels = rev(disease_labels_ordered)))

  plot_anno2 <- plot_anno %>%
    filter(as.character(disease_label) %in% disease_labels_ordered) %>%
    mutate(disease_label = factor(as.character(disease_label), levels = rev(disease_labels_ordered)))

  list(long = plot_long2, anno = plot_anno2)
}

make_plot <- function(plot_long, plot_anno, title_suffix = "selected disease endpoints") {

  n_rows <- n_distinct(plot_long$disease_id)
  max_abs <- max(abs(plot_long$log_hr_plot), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  p_heat <- ggplot(plot_long, aes(x = pair_label, y = disease_label)) +
    geom_tile(aes(fill = log_hr_plot), color = "grey86", linewidth = 0.15) +
    geom_text(aes(label = hr_label), size = 1.65, color = "black", na.rm = TRUE) +
    geom_text(
      data = plot_long %>% filter(diff_tile_label != ""),
      aes(label = diff_tile_label),
      size = 1.45,
      color = "black",
      fontface = "bold",
      vjust = 1.85,
      na.rm = TRUE
    ) +
    facet_grid(. ~ biomarker, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-max_abs, max_abs),
      oob = squish,
      name = "log(HR)"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Mortality L'EPOCH vs biological aging clocks: ", title_suffix),
      subtitle = paste0(
        "Only disease-clock-pair rows where L'EPOCH or BAG P < 0.05 / ",
        n_diseases, " / ", n_clock_pairs,
        " are shown. Tile text = HR per 1-SD score; * = significant association; ",
        "Delta+/Delta- = significant joint beta difference favoring L'EPOCH/aging clock. ",
        "Bonferroni P < ", signif(p_bonf, 3), "."
      )
    ) +
    theme_classic(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.6),
      axis.text.y = element_text(size = 6.0),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7),
      panel.grid = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )

  p_anno <- ggplot(plot_anno, aes(y = disease_label, x = 1)) +
    geom_text(aes(label = annotation), hjust = 0, size = 2.0) +
    scale_x_continuous(limits = c(1, 5.7), expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = "Sample size / follow-up / significance") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 8, hjust = 0),
      plot.margin = margin(t = 28, r = 5, b = 5, l = 4)
    )

  p_combined <- p_heat + p_anno + plot_layout(widths = c(5.6, 2.3))

  fig_height <- max(7, 0.23 * n_rows + 3.0)
  fig_width <- 22

  list(plot = p_combined, width = fig_width, height = fig_height, n_rows = n_rows)
}

# =========================
# 10. Optional full all-systems figure
# =========================

if (save_full_all_systems_pdf || save_full_all_systems_png) {
  full_plot <- make_plot(long_df, anno_df, title_suffix = "all systems")
  save_pdf_png(
    plot = full_plot$plot,
    file_prefix = file.path(out_dir, "LEPOCH_vs_BAG_significant_only_all_systems"),
    width = full_plot$width,
    height = full_plot$height,
    save_pdf = save_full_all_systems_pdf,
    save_png = save_full_all_systems_png
  )
}

# =========================
# 11. Paginated all-system figures
# =========================

if (save_paginated_all_systems) {
  page_dir <- file.path(out_dir, "paginated_all_systems")
  dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)

  ordered_diseases <- disease_rank$disease_label
  n_pages <- ceiling(length(ordered_diseases) / max_rows_per_page)
  message("Saving paginated all-system figures: ", n_pages, " pages")

  page_index <- data.frame()

  for (page_i in seq_len(n_pages)) {
    i1 <- (page_i - 1) * max_rows_per_page + 1
    i2 <- min(page_i * max_rows_per_page, length(ordered_diseases))
    disease_subset <- ordered_diseases[i1:i2]

    tmp <- reset_disease_levels(long_df, anno_df, disease_subset)
    p_obj <- make_plot(tmp$long, tmp$anno, title_suffix = paste0("all systems, page ", page_i, " of ", n_pages))

    prefix <- file.path(page_dir, sprintf("LEPOCH_vs_BAG_significant_only_all_systems_page_%03d", page_i))
    save_pdf_png(p_obj$plot, prefix, p_obj$width, p_obj$height, save_pdf = TRUE, save_png = TRUE)

    page_index <- bind_rows(
      page_index,
      data.frame(
        page = page_i,
        disease_label = as.character(disease_subset),
        stringsAsFactors = FALSE
      )
    )
  }

  fwrite(page_index, file.path(page_dir, "page_index.tsv"), sep = "\t")
}

# =========================
# 12. Body-system-specific figures
# =========================

if (save_by_body_system) {
  body_dir <- file.path(out_dir, "by_body_system")
  dir.create(body_dir, recursive = TRUE, showWarnings = FALSE)

  systems <- disease_rank %>%
    arrange(body_system_order, body_system) %>%
    pull(body_system) %>%
    unique()

  message("Saving body-system figures: ", length(systems), " systems")

  for (sys in systems) {
    sys_diseases <- disease_rank %>%
      filter(body_system == sys) %>%
      pull(disease_label)

    if (length(sys_diseases) == 0) next

    # If a body system is very large, also paginate it.
    n_sys_pages <- ceiling(length(sys_diseases) / max_rows_per_page)

    for (page_i in seq_len(n_sys_pages)) {
      i1 <- (page_i - 1) * max_rows_per_page + 1
      i2 <- min(page_i * max_rows_per_page, length(sys_diseases))
      disease_subset <- sys_diseases[i1:i2]

      tmp <- reset_disease_levels(long_df, anno_df, disease_subset)
      if (nrow(tmp$long) == 0) next

      title_txt <- if (n_sys_pages == 1) sys else paste0(sys, ", page ", page_i, " of ", n_sys_pages)
      p_obj <- make_plot(tmp$long, tmp$anno, title_suffix = title_txt)

      prefix <- file.path(
        body_dir,
        sprintf("LEPOCH_vs_BAG_significant_only_%s_page_%03d", make_safe_filename(sys), page_i)
      )

      save_pdf_png(p_obj$plot, prefix, p_obj$width, p_obj$height, save_pdf = TRUE, save_png = TRUE)
    }
  }
}

# =========================
# 13. Export exactly what was plotted
# =========================

fwrite(plot_base, file.path(out_dir, "LEPOCH_vs_BAG_significant_pair_rows.tsv"), sep = "\t")
fwrite(long_df, file.path(out_dir, "LEPOCH_vs_BAG_significant_pair_rows_long_for_plot.tsv"), sep = "\t")
fwrite(anno_df, file.path(out_dir, "LEPOCH_vs_BAG_disease_annotations_for_plot.tsv"), sep = "\t")
fwrite(disease_rank, file.path(out_dir, "LEPOCH_vs_BAG_disease_ordering.tsv"), sep = "\t")

summary_df <- plot_base %>%
  summarise(
    n_diseases_powered = n_diseases,
    n_clock_pairs = n_clock_pairs,
    bonferroni_p = p_bonf,
    n_significant_pair_rows = n(),
    n_disease_endpoints_retained = n_distinct(disease_id),
    n_pairs_retained = n_distinct(pair_id),
    n_diff_significant_rows = sum(diff_sig_bonf, na.rm = TRUE),
    max_rows_per_page = max_rows_per_page,
    save_full_all_systems_pdf = save_full_all_systems_pdf,
    save_full_all_systems_png = save_full_all_systems_png
  )

fwrite(summary_df, file.path(out_dir, "LEPOCH_vs_BAG_significant_only_summary.tsv"), sep = "\t")

message("Done. Outputs saved to: ", out_dir)
