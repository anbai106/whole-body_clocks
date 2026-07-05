# ============================================================
# Supplementary sensitivity analysis:
# Respiratory medication subtypes and immune metabolomics
# mortality-clock delta
#
# Main question:
#   Does the respiratory medication cluster remain associated
#   with immune delta clock after adjustment for baseline
#   asthma and COPD?
#
# Style:
#   Compatible with main medication-cluster plots.
#
# Significance encoding:
#   Solid circle        = Bonferroni-significant, p_bonferroni < 0.05
#   Light filled circle = nominal only, raw P < 0.05 but Bonferroni P >= 0.05
#   Empty circle        = non-significant, raw P >= 0.05
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(forcats)
})

# -----------------------------
# 1. Paths
# -----------------------------

results_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/delta_metabolomics_algorithmic_disease_onset/respiratory_medication_subtype_immune_delta_results"

out_dir <- file.path(results_dir, "plots_supplement")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lm_results_tsv <- file.path(
  results_dir,
  "immune_delta_respiratory_lm_results.tsv"
)

top_codes_tsv <- file.path(
  results_dir,
  "respiratory_cluster_top_respiratory_medication_codes.tsv"
)

subtype_counts_tsv <- file.path(
  results_dir,
  "respiratory_subtype_counts_in_immune_delta_dataset.tsv"
)

# -----------------------------
# 2. Settings
# -----------------------------

plot_bg <- "#FFF9E8"
strip_bg <- "#F3E7C8"

# Compatible Van Gogh-inspired palette
exposure_palette <- c(
  "Respiratory cluster vs no/minimal medication" = "#4E6FAE",
  "Inhaled corticosteroid" = "#3F7F5F",
  "Beta-agonist / bronchodilator" = "#C9A227",
  "Anticholinergic" = "#7E6AA2",
  "Leukotriene modifier" = "#8AA85A",
  "Systemic corticosteroid" = "#B65E16",
  "Other respiratory medication" = "#6B6B6B"
)

exposure_palette_light <- c(
  "Respiratory cluster vs no/minimal medication" = "#A9BCE3",
  "Inhaled corticosteroid" = "#A7C7AE",
  "Beta-agonist / bronchodilator" = "#E7D48A",
  "Anticholinergic" = "#B8AED5",
  "Leukotriene modifier" = "#C8D6A3",
  "Systemic corticosteroid" = "#E2AA72",
  "Other respiratory medication" = "#BDBDBD"
)

exposure_order <- c(
  "Respiratory cluster vs no/minimal medication",
  "Inhaled corticosteroid",
  "Beta-agonist / bronchodilator",
  "Anticholinergic",
  "Leukotriene modifier",
  "Systemic corticosteroid",
  "Other respiratory medication"
)

sig_levels <- c(
  "Bonferroni-significant",
  "Nominal only",
  "Not significant"
)

# -----------------------------
# 3. Helper functions
# -----------------------------

fmt_p <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(
      p < 1e-300, "<1e-300",
      ifelse(
        p < 0.001,
        formatC(p, format = "e", digits = 2),
        formatC(p, format = "f", digits = 3)
      )
    )
  )
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

clean_exposure_label <- function(x, model_name = NULL) {
  out <- dplyr::case_when(
    x == "respiratory_cluster_binary" ~ "Respiratory cluster vs no/minimal medication",
    x == "inhaled_corticosteroid" ~ "Inhaled corticosteroid",
    x == "beta_agonist_bronchodilator" ~ "Beta-agonist / bronchodilator",
    x == "anticholinergic" ~ "Anticholinergic",
    x == "leukotriene_modifier" ~ "Leukotriene modifier",
    x == "systemic_corticosteroid" ~ "Systemic corticosteroid",
    x == "other_respiratory" ~ "Other respiratory medication",
    TRUE ~ x
  )

  # Label the joint-model rows explicitly.
  if (!is.null(model_name)) {
    out <- ifelse(
      model_name == "joint_respiratory_subtype_model",
      paste0(out, " | joint model"),
      out
    )
  }

  out
}

# -----------------------------
# 4. Read linear-model sensitivity results
# -----------------------------

if (!file.exists(lm_results_tsv)) {
  stop("Missing LM results file: ", lm_results_tsv)
}

res <- fread(lm_results_tsv, fill = TRUE)

# Some accidental concatenated column names can occur if a TSV was copied
# through terminal output. The key columns below are the only ones required.
required_cols <- c(
  "model_name", "outcome", "exposure", "N", "N_exposed",
  "beta", "se", "ci_lo", "ci_hi", "p", "status"
)

missing_cols <- setdiff(required_cols, names(res))
if (length(missing_cols) > 0) {
  stop(
    "The LM results file is missing required columns:\n",
    paste(missing_cols, collapse = ", ")
  )
}

if (!"p_bonferroni" %in% names(res)) {
  res$p_bonferroni <- NA_real_
}

if (!"p_fdr_bh" %in% names(res)) {
  res$p_fdr_bh <- NA_real_
}

num_cols <- c(
  "N", "N_exposed", "N_unexposed",
  "mean_outcome_exposed", "mean_outcome_unexposed",
  "beta", "se", "ci_lo", "ci_hi", "p",
  "p_fdr_bh", "p_bonferroni",
  "r_squared", "adj_r_squared", "aic", "bic"
)

for (cc in intersect(num_cols, names(res))) {
  res[[cc]] <- safe_num(res[[cc]])
}

# Compute Bonferroni P if missing.
valid_p_idx <- which(res$status == "ok" & !is.na(res$p))
if (length(valid_p_idx) > 0 && all(is.na(res$p_bonferroni[valid_p_idx]))) {
  res$p_bonferroni[valid_p_idx] <- p.adjust(res$p[valid_p_idx], method = "bonferroni")
}

# -----------------------------
# 5. Prepare plot table
# -----------------------------

plot_tbl <- res %>%
  filter(status == "ok") %>%
  mutate(
    exposure_label = clean_exposure_label(exposure, model_name),

    model_type = case_when(
      model_name == "primary_respiratory_cluster_vs_no_minimal" ~ "Primary indication-adjusted model",
      model_name == "joint_respiratory_subtype_model" ~ "Joint subtype model",
      str_detect(model_name, "^subtype_") ~ "Single subtype model",
      TRUE ~ "Other model"
    ),

    bonf_sig = !is.na(p_bonferroni) & p_bonferroni < 0.05,
    nominal_only = !bonf_sig & !is.na(p) & p < 0.05,

    sig_group = case_when(
      bonf_sig ~ "Bonferroni-significant",
      nominal_only ~ "Nominal only",
      TRUE ~ "Not significant"
    ),
    sig_group = factor(sig_group, levels = sig_levels),

    exposure_base = clean_exposure_label(exposure),
    exposure_base = factor(exposure_base, levels = rev(exposure_order)),

    model_type = factor(
      model_type,
      levels = c(
        "Primary indication-adjusted model",
        "Single subtype model",
        "Joint subtype model"
      )
    ),

    annotation = paste0(
      "N=", comma(N),
      "; exposed=", comma(N_exposed),
      "; beta=", fmt_num(beta, 2),
      "; P=", fmt_p(p),
      "; Bonf P=", fmt_p(p_bonferroni)
    )
  )

small_n_tbl <- res %>%
  filter(status != "ok") %>%
  mutate(
    exposure_label = clean_exposure_label(exposure, model_name)
  )

fwrite(
  plot_tbl,
  file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_plot_table.tsv"),
  sep = "\t"
)

fwrite(
  small_n_tbl,
  file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_small_n_table.tsv"),
  sep = "\t"
)

# -----------------------------
# 6. Supplementary forest plot
# -----------------------------

if (nrow(plot_tbl) > 0) {

  x_min <- min(plot_tbl$ci_lo, na.rm = TRUE)
  x_max <- max(plot_tbl$ci_hi, na.rm = TRUE)
  x_pad <- 0.22 * (x_max - x_min)

  plot_tbl <- plot_tbl %>%
    mutate(
      x_text = x_max + 0.06 * (x_max - x_min)
    )

  p_sens <- ggplot(
    plot_tbl,
    aes(y = exposure_base)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "#4D4D4D",
      linewidth = 0.7
    ) +
    geom_segment(
      aes(
        x = ci_lo,
        xend = ci_hi,
        yend = exposure_base,
        color = exposure_base
      ),
      linewidth = 1.0,
      alpha = 0.95
    ) +

    # Bonferroni-significant: solid circle
    geom_point(
      data = plot_tbl %>% filter(sig_group == "Bonferroni-significant"),
      aes(x = beta, color = exposure_base),
      shape = 16,
      size = 3.9,
      stroke = 1.0,
      show.legend = FALSE
    ) +

    # Nominal-only: light filled circle
    geom_point(
      data = plot_tbl %>% filter(sig_group == "Nominal only"),
      aes(x = beta, color = exposure_base, fill = exposure_base),
      shape = 21,
      size = 3.7,
      stroke = 1.1,
      alpha = 0.55,
      show.legend = FALSE
    ) +

    # Non-significant: empty circle
    geom_point(
      data = plot_tbl %>% filter(sig_group == "Not significant"),
      aes(x = beta, color = exposure_base),
      shape = 1,
      size = 3.7,
      stroke = 1.1,
      show.legend = FALSE
    ) +

    # Shape legend
    geom_point(
      aes(x = beta, shape = sig_group),
      color = "#333333",
      fill = "#BDBDBD",
      size = 3.5,
      stroke = 1.0,
      alpha = 0,
      show.legend = TRUE
    ) +

    geom_text(
      aes(
        x = x_text,
        label = annotation,
        color = exposure_base
      ),
      hjust = 0,
      size = 3.0,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ model_type,
      ncol = 1,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = exposure_palette,
      name = "Respiratory medication exposure"
    ) +
    scale_fill_manual(
      values = exposure_palette_light,
      name = "Respiratory medication exposure"
    ) +
    scale_shape_manual(
      values = c(
        "Bonferroni-significant" = 16,
        "Nominal only" = 21,
        "Not significant" = 1
      ),
      breaks = sig_levels,
      drop = FALSE,
      name = "Significance"
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          alpha = 1,
          color = "#333333",
          fill = "#BDBDBD",
          size = 3.8
        )
      )
    ) +
    coord_cartesian(
      xlim = c(x_min - x_pad, x_max + 2.7 * x_pad),
      clip = "off"
    ) +
    labs(
      x = paste0(
        "Adjusted difference in immune delta clock age years\n",
        "versus comparison group"
      ),
      y = NULL,
      title = "Sensitivity analysis of respiratory medications and immune delta clock age",
      subtitle = paste0(
        "Models adjusted for baseline immune clock acceleration, chronological age, sex, ",
        "smoking, BMI, blood pressure, baseline asthma, and baseline COPD. ",
        "Solid = Bonferroni-significant; light-filled = nominal only; empty = not significant."
      )
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.background = element_rect(fill = plot_bg, color = NA),
      panel.background = element_rect(fill = plot_bg, color = NA),
      legend.background = element_rect(fill = plot_bg, color = NA),
      legend.box.background = element_rect(fill = plot_bg, color = NA),
      strip.background = element_rect(fill = strip_bg, color = NA),
      strip.text = element_text(face = "bold", size = 12, hjust = 0),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5),
      axis.text.y = element_text(face = "bold", size = 10),
      axis.title.x = element_text(size = 12),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(10, 230, 10, 10)
    )

  ggsave(
    file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_forest_plot.pdf"),
    p_sens,
    width = 14.8,
    height = 7.2
  )

  ggsave(
    file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_forest_plot.png"),
    p_sens,
    width = 14.8,
    height = 7.2,
    dpi = 300
  )

  ggsave(
    file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_forest_plot.svg"),
    p_sens,
    width = 14.8,
    height = 7.2
  )
}

# -----------------------------
# 7. Optional: top respiratory medication codes
# -----------------------------

if (file.exists(top_codes_tsv)) {

  top_codes <- fread(top_codes_tsv, fill = TRUE)

  if (all(c("meaning", "n_participants") %in% names(top_codes))) {

    top_codes_plot <- top_codes %>%
      mutate(
        n_participants = safe_num(n_participants),
        medication_label = str_to_sentence(meaning),
        medication_label = str_replace_all(medication_label, "\\s+", " "),
        medication_label = str_trunc(medication_label, width = 55)
      ) %>%
      arrange(desc(n_participants)) %>%
      slice_head(n = 20) %>%
      mutate(
        medication_label = fct_reorder(medication_label, n_participants)
      )

    fwrite(
      top_codes_plot,
      file.path(out_dir, "supp_respiratory_top_codes_for_plotting.tsv"),
      sep = "\t"
    )

    p_top <- ggplot(
      top_codes_plot,
      aes(x = n_participants, y = medication_label)
    ) +
      geom_col(
        fill = "#4E6FAE",
        alpha = 0.85,
        width = 0.72
      ) +
      geom_text(
        aes(label = comma(n_participants)),
        hjust = -0.15,
        size = 3.0
      ) +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.12)),
        labels = comma
      ) +
      labs(
        x = "Number of participants",
        y = NULL,
        title = "Top respiratory medication codes among respiratory-cluster participants",
        subtitle = "Baseline UK Biobank medication field 20003, instance 0"
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.background = element_rect(fill = plot_bg, color = NA),
        panel.background = element_rect(fill = plot_bg, color = NA),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10.5),
        axis.text.y = element_text(size = 9.5),
        axis.title.x = element_text(size = 11)
      )

    ggsave(
      file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_top_codes.pdf"),
      p_top,
      width = 9.5,
      height = 6.4
    )

    ggsave(
      file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_top_codes.png"),
      p_top,
      width = 9.5,
      height = 6.4,
      dpi = 300
    )

    ggsave(
      file.path(out_dir, "supp_respiratory_immune_delta_sensitivity_top_codes.svg"),
      p_top,
      width = 9.5,
      height = 6.4
    )
  }
}

# -----------------------------
# 8. Console summary
# -----------------------------

cat("\n============================================================\n")
cat("Finished respiratory medication sensitivity plots.\n")
cat("Input directory:\n", results_dir, "\n\n")
cat("Output directory:\n", out_dir, "\n\n")

cat("Bonferroni-significant sensitivity signals:\n")
print(
  plot_tbl %>%
    filter(bonf_sig) %>%
    arrange(p_bonferroni) %>%
    select(
      model_type,
      exposure_label,
      N,
      N_exposed,
      beta,
      ci_lo,
      ci_hi,
      p,
      p_bonferroni
    )
)

cat("\nNominal-only sensitivity signals:\n")
print(
  plot_tbl %>%
    filter(nominal_only) %>%
    arrange(p) %>%
    select(
      model_type,
      exposure_label,
      N,
      N_exposed,
      beta,
      ci_lo,
      ci_hi,
      p,
      p_bonferroni
    )
)

cat("\nSmall-N / skipped respiratory subtype rows:\n")
print(
  small_n_tbl %>%
    select(
      model_name,
      exposure_label,
      status,
      error,
      N_exposed,
      min_n_exposure
    )
)

cat("============================================================\n")