# ============================================================
# Mortality clocks vs AI biomarkers: LDSC genetic correlation
# Organ-resolved heatmap
#
# Revision:
#   1. DNE subtypes are treated as brain subtypes.
#   2. MAE biomarkers are labeled as subtypes.
#   3. Significant cells are annotated with gc_mean and stars.
#   4. Multiple comparison correction:
#        *   nominal: P < 0.05
#        **  Bonferroni over mortality clocks: P < 0.05 / number of mortality clocks
#        *** Bonferroni over all AI-biomarker tests:
#            P < 0.05 / number of valid mortality clock x AI biomarker tests
#   5. Axis labels are forced to one-line text.
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

out_prefix <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/mortality_clocks_vs_AI_biomarkers_LDSC_rg"

dat <- readr::read_tsv(infile, show_col_types = FALSE)

# Keep only previous AI biomarkers, not FinnGen/PGC disease endpoints
ai <- dat %>%
  filter(analysis_group == "AI_biomarker")

message(
  "Using ", nrow(ai), " rows: ",
  n_distinct(ai$mortality_clock), " mortality clocks x ",
  n_distinct(ai$target_display), " AI biomarkers"
)

# ----------------------------
# 2. Organ/system ordering
# ----------------------------
# DNE is no longer a separate organ block.
# DNE rows are assigned to Brain because they are brain subtypes.
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

# Organ-first ordering of mortality clocks.
# All labels are one-line text.
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

missing_clocks <- setdiff(unique(ai$mortality_clock), clock_key$mortality_clock)

if (length(missing_clocks) > 0) {
  stop(
    "These mortality clocks are missing from clock_key: ",
    paste(missing_clocks, collapse = ", ")
  )
}

clock_order <- clock_key %>%
  arrange(clock_block, clock_suborder) %>%
  pull(clock_label)

# ----------------------------
# 3. AI biomarker one-line labeling
# ----------------------------
pretty_target_label <- function(target_display, target_source) {
  
  out <- target_display
  
  # Pan-disease MAE: explicitly label as subtype, one line
  out <- ifelse(
    target_source == "Pan_disease_MAE",
    out %>%
      str_replace("^R([0-9]+) brain MAE$", "MAE subtype R\\1 brain") %>%
      str_replace("^R([0-9]+) eye MAE$",   "MAE subtype R\\1 eye") %>%
      str_replace("^R([0-9]+) heart MAE$", "MAE subtype R\\1 heart"),
    out
  )
  
  # DNE: explicitly label as brain subtype, one line
  out <- ifelse(
    target_source == "DNE",
    out %>%
      str_replace("^AD SurrealGAN ([0-9]+)$", "DNE brain subtype AD S-GAN \\1") %>%
      str_replace("^ASD ([0-9]+)$",           "DNE brain subtype ASD \\1") %>%
      str_replace("^LLD ([0-9]+)$",           "DNE brain subtype LLD \\1") %>%
      str_replace("^SCZ ([0-9]+)$",           "DNE brain subtype SCZ \\1"),
    out
  )
  
  # Brain MRIBAG components
  out <- out %>%
    str_replace("^GM brain MRIBAG$", "Brain MRI BAG GM") %>%
    str_replace("^WM brain MRIBAG$", "Brain MRI BAG WM") %>%
    str_replace("^FC brain MRIBAG$", "Brain MRI BAG FC")
  
  # Multi-organ MRIBAG
  out <- out %>%
    str_replace("^MRIBAG adipose$",  "MRI BAG adipose") %>%
    str_replace("^MRIBAG brain$",    "MRI BAG brain") %>%
    str_replace("^MRIBAG heart$",    "MRI BAG heart") %>%
    str_replace("^MRIBAG kidney$",   "MRI BAG kidney") %>%
    str_replace("^MRIBAG liver$",    "MRI BAG liver") %>%
    str_replace("^MRIBAG pancreas$", "MRI BAG pancreas") %>%
    str_replace("^MRIBAG spleen$",   "MRI BAG spleen")
  
  # ProtBAG
  out <- out %>%
    str_replace("^ProtBAG Brain$",                 "ProtBAG brain") %>%
    str_replace("^ProtBAG Endocrine$",             "ProtBAG endocrine") %>%
    str_replace("^ProtBAG Eye$",                   "ProtBAG eye") %>%
    str_replace("^ProtBAG Heart$",                 "ProtBAG heart") %>%
    str_replace("^ProtBAG Hepatic$",               "ProtBAG hepatic") %>%
    str_replace("^ProtBAG Immune$",                "ProtBAG immune") %>%
    str_replace("^ProtBAG Pulmonary$",             "ProtBAG pulmonary") %>%
    str_replace("^ProtBAG Renal$",                 "ProtBAG renal") %>%
    str_replace("^ProtBAG Reproductive_female$",   "ProtBAG reproductive female") %>%
    str_replace("^ProtBAG Reproductive_male$",     "ProtBAG reproductive male") %>%
    str_replace("^ProtBAG Skin$",                  "ProtBAG skin")
  
  # MetBAG
  out <- out %>%
    str_replace("^MetBAG Digestive$",  "MetBAG digestive") %>%
    str_replace("^MetBAG Endocrine$",  "MetBAG endocrine") %>%
    str_replace("^MetBAG Hepatic$",    "MetBAG hepatic") %>%
    str_replace("^MetBAG Immune$",     "MetBAG immune") %>%
    str_replace("^MetBAG Metabolic$",  "MetBAG metabolic")
  
  # Safety: remove any remaining accidental line breaks
  out <- str_replace_all(out, "\\n", " ")
  out <- str_squish(out)
  
  out
}

# ----------------------------
# 4. AI biomarker organ assignment and ordering
# ----------------------------
target_source_levels <- c(
  "MRIBAG",
  "Brain_MRIBAG",
  "ProtBAG",
  "MetBAG",
  "Pan_disease_MAE",
  "DNE"
)

target_order_tbl <- ai %>%
  mutate(first_row = row_number()) %>%
  distinct(target_display, target_source, target_family, .keep_all = TRUE) %>%
  mutate(
    target_block = case_when(
      # DNE belongs to brain
      target_source == "DNE" ~ "Brain",
      
      # MAE subtype organs
      target_source == "Pan_disease_MAE" &
        str_detect(target_display, regex("brain", ignore_case = TRUE)) ~ "Brain",
      target_source == "Pan_disease_MAE" &
        str_detect(target_display, regex("eye", ignore_case = TRUE)) ~ "Eye",
      target_source == "Pan_disease_MAE" &
        str_detect(target_display, regex("heart", ignore_case = TRUE)) ~ "Heart",
      
      # General organ assignment
      str_detect(target_display, regex("brain|\\bWM\\b|\\bFC\\b|\\bGM\\b", ignore_case = TRUE)) ~ "Brain",
      str_detect(target_display, regex("eye", ignore_case = TRUE)) ~ "Eye",
      str_detect(target_display, regex("pulmonary", ignore_case = TRUE)) ~ "Pulmonary",
      str_detect(target_display, regex("heart", ignore_case = TRUE)) ~ "Heart",
      str_detect(target_display, regex("liver|hepatic", ignore_case = TRUE)) ~ "Liver / hepatic",
      str_detect(target_display, regex("kidney|renal", ignore_case = TRUE)) ~ "Kidney / renal",
      str_detect(target_display, regex("pancreas", ignore_case = TRUE)) ~ "Pancreas",
      str_detect(target_display, regex("spleen", ignore_case = TRUE)) ~ "Spleen",
      str_detect(target_display, regex("immune", ignore_case = TRUE)) ~ "Immune",
      str_detect(target_display, regex("endocrine", ignore_case = TRUE)) ~ "Endocrine",
      str_detect(target_display, regex("digestive", ignore_case = TRUE)) ~ "Digestive",
      str_detect(target_display, regex("metabolic", ignore_case = TRUE)) ~ "Metabolic",
      str_detect(target_display, regex("adipose", ignore_case = TRUE)) ~ "Adipose",
      str_detect(target_display, regex("skin", ignore_case = TRUE)) ~ "Skin",
      str_detect(target_display, regex("reproductive", ignore_case = TRUE)) ~ "Reproductive",
      TRUE ~ NA_character_
    ),
    
    target_block = factor(target_block, levels = organ_levels),
    
    target_source_rank = match(target_source, target_source_levels),
    target_source_rank = replace_na(target_source_rank, 999),
    
    target_suborder = case_when(
      target_display == "MRIBAG brain" ~ 1,
      target_display == "GM brain MRIBAG" ~ 2,
      target_display == "WM brain MRIBAG" ~ 3,
      target_display == "FC brain MRIBAG" ~ 4,
      
      # MAE subtype ordering
      str_detect(target_display, "^R[0-9]+") ~ readr::parse_number(target_display),
      
      # DNE subtype ordering
      str_detect(target_display, "^AD SurrealGAN") ~ 100 + readr::parse_number(target_display),
      str_detect(target_display, "^ASD") ~ 200 + readr::parse_number(target_display),
      str_detect(target_display, "^LLD") ~ 300 + readr::parse_number(target_display),
      str_detect(target_display, "^SCZ") ~ 400 + readr::parse_number(target_display),
      
      TRUE ~ suppressWarnings(readr::parse_number(target_display))
    ),
    
    target_suborder = replace_na(target_suborder, 999),
    target_label = pretty_target_label(target_display, target_source)
  )

unassigned_targets <- target_order_tbl %>%
  filter(is.na(target_block)) %>%
  pull(target_display)

if (length(unassigned_targets) > 0) {
  stop(
    "These AI biomarkers could not be assigned to an organ block: ",
    paste(unassigned_targets, collapse = ", ")
  )
}

target_order_tbl <- target_order_tbl %>%
  arrange(target_block, target_source_rank, target_suborder, first_row)

# Reverse factor levels so the desired organ order appears top-to-bottom
target_order <- rev(target_order_tbl$target_label)

# ----------------------------
# 5. Multiple comparison correction
# ----------------------------
# User-requested correction levels:
#   1. Nominal
#   2. Bonferroni over mortality clocks
#   3. Bonferroni over all AI-biomarker tests in the heatmap

n_mortality_clocks <- n_distinct(ai$mortality_clock)
n_ai_biomarkers <- n_distinct(ai$target_display)

if (n_mortality_clocks != 22) {
  warning(
    "Expected 22 mortality clocks, but detected ",
    n_mortality_clocks,
    ". The Bonferroni clock-level threshold will use the detected number."
  )
}

# Correction across full heatmap.
# This uses mortality clocks x AI biomarkers as the test family.
n_all_ai_tests <- n_mortality_clocks * n_ai_biomarkers

p_nominal <- 0.05
p_bonf_clock <- 0.05 / n_mortality_clocks
p_bonf_all_ai <- 0.05 / n_all_ai_tests

correction_tbl <- tibble(
  correction = c(
    "Nominal",
    "Bonferroni across mortality clocks",
    "Bonferroni across all mortality clock x AI biomarker tests"
  ),
  denominator = c(
    1,
    n_mortality_clocks,
    n_all_ai_tests
  ),
  p_threshold = c(
    p_nominal,
    p_bonf_clock,
    p_bonf_all_ai
  )
)

print(correction_tbl)

readr::write_tsv(
  correction_tbl,
  paste0(out_prefix, "_multiple_testing_thresholds.tsv")
)

# ----------------------------
# 6. Plotting table
# ----------------------------
plot_df <- ai %>%
  left_join(clock_key, by = "mortality_clock") %>%
  left_join(
    target_order_tbl %>%
      select(target_display, target_block, target_label),
    by = "target_display"
  ) %>%
  mutate(
    clock_label = factor(clock_label, levels = clock_order),
    target_label = factor(target_label, levels = target_order),
    clock_block = factor(clock_block, levels = organ_levels),
    target_block = factor(target_block, levels = organ_levels),
    
    valid_ldsc = !is.na(P) & !is.na(gc_mean),
    
    sig_nominal = valid_ldsc & P < p_nominal,
    sig_clock = valid_ldsc & P < p_bonf_clock,
    sig_all_ai = valid_ldsc & P < p_bonf_all_ai,
    
    sig_symbol = case_when(
      sig_all_ai ~ "***",
      sig_clock ~ "**",
      sig_nominal ~ "*",
      TRUE ~ ""
    ),
    
    # Annotate gc_mean only for significant cells
    gc_label = case_when(
      sig_nominal ~ paste0(sprintf("%.2f", gc_mean), "\n", sig_symbol),
      TRUE ~ ""
    ),
    
    label_color = case_when(
      is.na(gc_mean) ~ "black",
      abs(gc_mean) >= 0.45 ~ "white",
      TRUE ~ "black"
    )
  )

# Symmetric color range for LDSC genetic correlation
rg_lim <- max(abs(plot_df$gc_mean), na.rm = TRUE)
rg_lim <- min(1, max(0.25, ceiling(rg_lim * 10) / 10))

# ----------------------------
# 7. Fancy organ-resolved heatmap
# ----------------------------
caption_txt <- paste0(
  "* nominal P < 0.05; ",
  "** Bonferroni over ", n_mortality_clocks, " mortality clocks, P < ",
  signif(p_bonf_clock, 3),
  "; ",
  "*** Bonferroni over all AI-biomarker tests, P < ",
  signif(p_bonf_all_ai, 3),
  ". Numbers inside significant cells show LDSC mean genetic correlation. ",
  "Grey cells indicate unavailable LDSC estimates."
)

p <- ggplot(plot_df, aes(x = clock_label, y = target_label)) +
  
  # Main heatmap
  geom_tile(
    aes(fill = gc_mean),
    width = 0.96,
    height = 0.96,
    color = "white",
    linewidth = 0.18
  ) +
  
  # Thin border for nominally significant cells
  geom_tile(
    data = plot_df %>% filter(sig_nominal),
    fill = NA,
    color = "grey20",
    linewidth = 0.22,
    width = 0.96,
    height = 0.96
  ) +
  
  # Stronger border for global Bonferroni-significant cells
  geom_tile(
    data = plot_df %>% filter(sig_all_ai),
    fill = NA,
    color = "black",
    linewidth = 0.55,
    width = 0.96,
    height = 0.96
  ) +
  
  # Annotate significant cells with gc_mean and significance stars
  geom_text(
    aes(label = gc_label, color = label_color),
    size = 1.85,
    lineheight = 0.78,
    fontface = "bold",
    na.rm = TRUE
  ) +
  
  scale_color_identity() +
  
  scale_fill_gradient2(
    name = "LDSC rg",
    low = "#2166AC",
    mid = "#FAFAFA",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-rg_lim, rg_lim),
    breaks = scales::pretty_breaks(n = 5),
    oob = scales::squish,
    na.value = "#E6E6E6"
  ) +
  
  facet_grid(
    rows = vars(target_block),
    cols = vars(clock_block),
    scales = "free",
    space = "free",
    switch = "y"
  ) +
  
  labs(
    title = "Genetic correlation between mortality clocks and AI biomarkers",
    subtitle = paste0(
      "Rows are AI biomarkers, including biological aging clocks, MAE subtypes, and DNE brain subtypes; ",
      "columns are mortality clocks. Both axes are ordered by organ/system."
    ),
    x = NULL,
    y = NULL,
    caption = caption_txt
  ) +
  
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = grid::unit(60, "mm"),
      barwidth = grid::unit(5, "mm")
    )
  ) +
  
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    
    panel.spacing.x = grid::unit(0.14, "lines"),
    panel.spacing.y = grid::unit(0.30, "lines"),
    
    strip.placement = "outside",
    strip.background = element_rect(
      fill = "grey95",
      color = "grey75",
      linewidth = 0.25
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 7.5,
      margin = margin(t = 3, b = 3)
    ),
    strip.text.y.left = element_text(
      face = "bold",
      angle = 0,
      size = 7.5,
      margin = margin(l = 3, r = 3)
    ),
    
    # One-line x-axis labels
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 6.8,
      face = "bold"
    ),
    
    # One-line y-axis labels
    axis.text.y = element_text(
      size = 6.3,
      lineheight = 0.95
    ),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.6, color = "grey25"),
    plot.caption = element_text(
      hjust = 0,
      size = 7.6,
      color = "grey35",
      lineheight = 1.05
    ),
    
    plot.margin = margin(8, 12, 8, 8)
  )

print(p)

# ----------------------------
# 8. Save high-resolution outputs
# ----------------------------
ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = p,
  width = 22,
  height = 15,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = p,
  width = 22,
  height = 15,
  units = "in",
  dpi = 500,
  bg = "white"
)

# Save the plotting dataframe for checking
readr::write_tsv(
  plot_df,
  paste0(out_prefix, "_plotting_table.tsv")
)