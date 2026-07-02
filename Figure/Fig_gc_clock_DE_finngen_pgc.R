# ============================================================
# Panel D: Mortality clocks vs FinnGen / PGC disease endpoints
# LDSC genetic-correlation disease-system landscape
#
# Each bubble summarizes Bonferroni-significant disease endpoints:
#   size = number of significant disease endpoints
#   fill = mean LDSC genetic correlation among significant endpoints
#
# Significance threshold:
#   P < 0.05 / number of unique diseases / 22 mortality clocks
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(grid)
})

# ----------------------------
# 1. Input / output
# ----------------------------
infile <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/LDSC_gc_mortality_clocks_all_targets.tsv"

out_prefix <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/panelD_mortality_clocks_FinnGen_PGC_disease_LDSC"

dat <- readr::read_tsv(infile, show_col_types = FALSE)

# ----------------------------
# 2. Organ/system ordering for mortality clocks
# ----------------------------
organ_levels <- c(
  "Brain",
  "Eye",
  "Pulmonary",
  "Heart",
  "Liver / hepatic",
  "Kidney / renal",
  "Pancreas",
  "Spleen",
  "Immune",
  "Endocrine",
  "Digestive",
  "Metabolic",
  "Adipose",
  "Skin",
  "Reproductive"
)

clock_key <- tribble(
  ~mortality_clock,                    ~clock_block,        ~clock_label,                         ~clock_suborder,
  "brain_mri",                         "Brain",             "Brain MRI",                          1,
  "Brain_proteomics",                  "Brain",             "Brain protein",                      2,

  "Eye_proteomics",                    "Eye",               "Eye protein",                        1,

  "Pulmonary_proteomics",              "Pulmonary",         "Pulmonary protein",                  1,

  "heart_mri",                         "Heart",             "Heart MRI",                          1,
  "Heart_proteomics",                  "Heart",             "Heart protein",                      2,

  "liver_mri",                         "Liver / hepatic",   "Liver MRI",                          1,
  "Hepatic_proteomics",                "Liver / hepatic",   "Hepatic protein",                    2,
  "Hepatic_metabolomics",              "Liver / hepatic",   "Hepatic metabolite",                 3,

  "kidney_mri",                        "Kidney / renal",    "Kidney MRI",                         1,
  "Renal_proteomics",                  "Kidney / renal",    "Renal protein",                      2,

  "pancreas_mri",                      "Pancreas",          "Pancreas MRI",                       1,

  "spleen_mri",                        "Spleen",            "Spleen MRI",                         1,

  "Immune_proteomics",                 "Immune",            "Immune protein",                     1,
  "Immune_metabolomics",               "Immune",            "Immune metabolite",                  2,

  "Endocrine_proteomics",              "Endocrine",         "Endocrine protein",                  1,
  "Endocrine_metabolomics",            "Endocrine",         "Endocrine metabolite",               2,

  "Digestive_metabolomics",            "Digestive",         "Digestive metabolite",               1,

  "Metabolic_metabolomics",            "Metabolic",         "Metabolic metabolite",               1,

  "adipose_mri",                       "Adipose",           "Adipose MRI",                        1,

  "Skin_proteomics",                   "Skin",              "Skin protein",                       1,

  "Reproductive_female_proteomics",    "Reproductive",      "Reproductive female protein",        1,
  "Reproductive_male_proteomics",      "Reproductive",      "Reproductive male protein",          2
) %>%
  mutate(clock_block = factor(clock_block, levels = organ_levels))

clock_order <- clock_key %>%
  arrange(clock_block, clock_suborder) %>%
  pull(clock_label)

# ----------------------------
# 3. Disease endpoint data
# ----------------------------
disease <- dat %>%
  filter(
    analysis_group == "Disease_endpoint",
    target_source %in% c("FinnGen", "PGC")
  ) %>%
  mutate(
    endpoint_id = paste0(target_source, "::", target_id),
    endpoint_label_raw = target_display,
    valid_ldsc = !is.na(P) & !is.na(gc_mean)
  ) %>%
  distinct(mortality_clock, endpoint_id, .keep_all = TRUE)

if (nrow(disease) == 0) {
  stop("No FinnGen/PGC disease-endpoint rows found in the input file.")
}

# ----------------------------
# 4. Disease chapter assignment
# ----------------------------
assign_disease_chapter <- function(target_source, target_display) {

  case_when(
    target_source == "PGC" ~ "PGC brain/psychiatric",

    str_detect(target_display, "^(AD_LO|G6|NEURODEG|MIGRAINE|TRAUMBRAIN)") ~ "Neurologic",
    str_detect(target_display, "^(F5|KRA|AUD|ALCOHOL|ANTIDEPRESSANTS|SLEEP)") ~ "Psychiatric / sleep",

    str_detect(target_display, "^H7") ~ "Eye",
    str_detect(target_display, "^(H8|FE)") ~ "Ear / sensory",

    str_detect(target_display, "^(I9|FG|CARDIAC|OTHER_SYSTCON|RX_ANTIHYP|RX_STATIN)") ~ "Cardiovascular",

    str_detect(target_display, "^(J10|ASTHMA|COPD|BRONCHITIS|CPAP|NIV|PULM|VOCALCORD|INFLUENZA)") ~ "Respiratory",

    str_detect(target_display, "^(K11|ABDOM|APPEND|DENTAL|TEMPOROMANDIB|RX_CROHN)") ~ "Digestive / oral",

    str_detect(target_display, "^(E4|DM|T1D|T2D|KELA_DIAB|HYPOTHY|THYRO|GOUT|AUTOIMMUNE)") ~ "Endocrine / immune",

    str_detect(target_display, "^(D3|BLEEDING)") ~ "Blood",

    str_detect(target_display, "^(C3|C$|CD2)") ~ "Cancer / benign neoplasm",

    str_detect(target_display, "^(L12|ALLERG|NONALLERG|POLLEN|DRY)") ~ "Skin / allergy",

    str_detect(target_display, "^(M13|RHEU|RHEUMA|SPONDYLO|JOINTPAIN|PAIN|PRIM|FALLS|ST19|VWXY20|RX_CODEINE|RX_PARACETAMOL|RX_GLUCO)") ~ "MSK / injury",

    str_detect(target_display, "^(N14|R18|RX_INFERTILITY|Z21)") ~ "GU / reproductive",

    str_detect(target_display, "^(O15|GEST)") ~ "Pregnancy",

    str_detect(target_display, "^(AB1)") ~ "Infectious",

    TRUE ~ "Other"
  )
}

chapter_levels <- c(
  "PGC brain/psychiatric",
  "Neurologic",
  "Psychiatric / sleep",
  "Eye",
  "Ear / sensory",
  "Cardiovascular",
  "Respiratory",
  "Digestive / oral",
  "Endocrine / immune",
  "Blood",
  "Cancer / benign neoplasm",
  "Skin / allergy",
  "MSK / injury",
  "GU / reproductive",
  "Pregnancy",
  "Infectious",
  "Other"
)

pretty_endpoint_label <- function(target_source, target_display) {
  out <- ifelse(
    target_source == "PGC",
    paste0("PGC ", target_display),
    str_replace_all(target_display, "_", " ")
  )

  out <- str_replace_all(out, "EXALLC", "")
  out <- str_replace_all(out, "INCLAVO", "")
  out <- str_replace_all(out, "EXMORE", "")
  out <- str_squish(out)

  out
}

disease <- disease %>%
  mutate(
    disease_chapter = assign_disease_chapter(target_source, target_display),
    disease_chapter = factor(disease_chapter, levels = chapter_levels),
    endpoint_label = pretty_endpoint_label(target_source, target_display)
  )

# ----------------------------
# 5. Multiple-comparison threshold
# ----------------------------
n_unique_diseases <- n_distinct(disease$endpoint_id)

# User-specified denominator uses 22 mortality clocks.
n_mortality_clocks_for_correction <- 22

n_mortality_clocks_observed <- n_distinct(disease$mortality_clock)

if (n_mortality_clocks_observed != 22) {
  warning(
    "Detected ", n_mortality_clocks_observed,
    " mortality clocks in the disease endpoint rows, but using 22 for correction as requested."
  )
}

p_disease_bonf <- 0.05 / n_unique_diseases / n_mortality_clocks_for_correction
logp_disease_bonf <- -log10(p_disease_bonf)

message("Number of unique disease endpoints: ", n_unique_diseases)
message("Disease endpoint Bonferroni threshold: ", signif(p_disease_bonf, 4))
message("-log10 threshold: ", signif(logp_disease_bonf, 4))

threshold_tbl <- tibble(
  n_unique_diseases = n_unique_diseases,
  n_mortality_clocks_for_correction = n_mortality_clocks_for_correction,
  n_tests = n_unique_diseases * n_mortality_clocks_for_correction,
  p_threshold = p_disease_bonf,
  neg_log10_p_threshold = logp_disease_bonf
)

readr::write_tsv(
  threshold_tbl,
  paste0(out_prefix, "_threshold.tsv")
)

# ----------------------------
# 6. Significant endpoint-level results
# ----------------------------
disease_sig <- disease %>%
  left_join(clock_key, by = "mortality_clock") %>%
  mutate(
    clock_label = factor(clock_label, levels = clock_order),
    clock_block = factor(clock_block, levels = organ_levels),
    neg_log10_P = -log10(P),
    sig_disease_bonf = valid_ldsc & P < p_disease_bonf
  ) %>%
  filter(sig_disease_bonf) %>%
  arrange(P)

readr::write_tsv(
  disease_sig,
  paste0(out_prefix, "_significant_endpoint_level_results.tsv")
)

message("Number of significant mortality clock x disease endpoint pairs: ", nrow(disease_sig))
message("Number of unique significant disease endpoints: ", n_distinct(disease_sig$endpoint_id))
message("Number of mortality clocks with at least one significant disease endpoint: ", n_distinct(disease_sig$mortality_clock))

# ----------------------------
# 7. Summarize for compact Panel D
# ----------------------------
panel_df <- disease_sig %>%
  group_by(
    clock_block,
    clock_label,
    disease_chapter
  ) %>%
  summarise(
    n_sig = n(),
    mean_rg = mean(gc_mean, na.rm = TRUE),
    median_rg = median(gc_mean, na.rm = TRUE),
    max_abs_rg = max(abs(gc_mean), na.rm = TRUE),
    min_P = min(P, na.rm = TRUE),
    max_neg_log10_P = max(neg_log10_P, na.rm = TRUE),
    top_endpoint = endpoint_label[which.min(P)],
    .groups = "drop"
  ) %>%
  mutate(
    disease_chapter = factor(disease_chapter, levels = rev(chapter_levels)),
    clock_label = factor(clock_label, levels = clock_order),
    clock_block = factor(clock_block, levels = organ_levels),
    label_color = if_else(abs(mean_rg) >= 0.45, "white", "black")
  )

readr::write_tsv(
  panel_df,
  paste0(out_prefix, "_panelD_summary_table.tsv")
)

# ----------------------------
# 8. Build Panel D
# ----------------------------
if (nrow(panel_df) == 0) {

  pD <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.15,
      label = paste0(
        "No FinnGen/PGC disease endpoint associations passed\n",
        "P < 0.05 / ", n_unique_diseases, " diseases / 22 clocks = ",
        signif(p_disease_bonf, 3)
      ),
      size = 5,
      lineheight = 1.05,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = 0,
      y = -0.15,
      label = "Check whether disease-endpoint LDSC estimates are available in the input file.",
      size = 3.5,
      color = "grey35"
    ) +
    labs(tag = "D") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.tag = element_text(face = "bold", size = 22),
      plot.tag.position = c(0.01, 0.98),
      plot.margin = margin(10, 10, 10, 10)
    )

} else {

  rg_lim <- max(abs(panel_df$mean_rg), na.rm = TRUE)
  rg_lim <- min(1, max(0.20, ceiling(rg_lim * 10) / 10))

  pD <- ggplot(panel_df, aes(x = clock_label, y = disease_chapter)) +

    geom_point(
      aes(size = n_sig, fill = mean_rg),
      shape = 21,
      color = "grey15",
      stroke = 0.28,
      alpha = 0.96
    ) +

    geom_text(
      aes(label = n_sig, color = label_color),
      size = 2.25,
      fontface = "bold",
      show.legend = FALSE
    ) +

    scale_color_identity() +

    scale_fill_gradient2(
      name = "Mean LDSC rg",
      low = "#2166AC",
      mid = "#FAFAFA",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-rg_lim, rg_lim),
      breaks = scales::pretty_breaks(n = 5),
      oob = scales::squish
    ) +

    scale_size_area(
      name = "# significant\nendpoints",
      max_size = 11,
      breaks = scales::pretty_breaks(n = 4)
    ) +

    facet_grid(
      cols = vars(clock_block),
      scales = "free_x",
      space = "free_x",
      switch = "x"
    ) +

    labs(
      tag = "D",
      title = "Disease-endpoint genetic correlation landscape",
      subtitle = paste0(
        "FinnGen and PGC endpoints passing P < 0.05 / ",
        n_unique_diseases, " diseases / 22 mortality clocks = ",
        signif(p_disease_bonf, 3),
        ". Bubble size shows the number of significant endpoints; color shows mean signed genetic correlation."
      ),
      x = NULL,
      y = NULL,
      caption = "Numbers inside bubbles indicate significant mortality clock x disease endpoint pairs within each disease system. Endpoint-level results are saved to the accompanying TSV file."
    ) +

    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(38, "mm"),
        barwidth = grid::unit(5, "mm")
      ),
      size = guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(fill = "grey75")
      )
    ) +

    theme_minimal(base_size = 9) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),

      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.22),

      panel.spacing.x = grid::unit(0.12, "lines"),

      strip.placement = "outside",
      strip.background = element_rect(
        fill = "grey95",
        color = "grey78",
        linewidth = 0.25
      ),
      strip.text.x = element_text(
        face = "bold",
        size = 7.4,
        margin = margin(t = 3, b = 3)
      ),

      axis.text.x = element_text(
        angle = 58,
        hjust = 1,
        vjust = 1,
        size = 6.6,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = 7.3,
        face = "bold",
        color = "grey20"
      ),

      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 8.5),
      legend.text = element_text(size = 7.5),

      plot.tag = element_text(face = "bold", size = 22),
      plot.tag.position = c(0.005, 0.995),

      plot.title = element_text(face = "bold", size = 14.5),
      plot.subtitle = element_text(size = 8.8, color = "grey25", lineheight = 1.05),
      plot.caption = element_text(size = 7.3, color = "grey35", hjust = 0, lineheight = 1.05),

      plot.margin = margin(8, 12, 8, 8)
    )
}

print(pD)

# ----------------------------
# 9. Save Panel D
# ----------------------------
ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = pD,
  width = 17.8,
  height = 6.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = pD,
  width = 17.8,
  height = 6.2,
  units = "in",
  dpi = 500,
  bg = "white"
)

# Optional editable SVG output
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, ".svg"),
    plot = pD,
    width = 17.8,
    height = 6.2,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}