# ============================================================
# Fig. 2E-style forest plot:
# Mortality EPOCH clocks vs FinnGen / PGC disease endpoints
# LDSC genetic-correlation signals only
#
# Shows only Bonferroni-corrected significant DE-EPOCH signals:
#   P < 0.05 / 527
#
# One row = one significant mortality EPOCH clock x disease endpoint pair
# x-axis = LDSC genetic correlation rg
# horizontal line = rg +/- 1.96 * SE
#
# IMPORTANT:
#   In the current LDSC result file, the SE column is named gc_std.
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

out_prefix <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/Fig2E_DE_EPOCH_LDSC_Bonf_forest_52signals"

dat <- readr::read_tsv(infile, show_col_types = FALSE)

# ----------------------------
# 2. Helper functions
# ----------------------------
safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

get_first_col <- function(tbl, candidates, default = NA_real_) {
  hits <- intersect(candidates, names(tbl))
  if (length(hits) == 0) {
    return(rep(default, nrow(tbl)))
  }
  tbl[[hits[1]]]
}

get_first_col_name <- function(tbl, candidates) {
  hits <- intersect(candidates, names(tbl))
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[1]
}

clean_text <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_replace_all("EXALLC|INCLAVO|EXMORE", "") %>%
    str_squish()
}

format_pval <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 1e-300 ~ "<1e-300",
    p < 1e-3 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ sprintf("%.3f", p)
  )
}

shorten_endpoint <- function(x) {
  x %>%
    str_replace_all("Alzheimer's disease", "AD") %>%
    str_replace_all("Autism spectrum disorder", "ASD") %>%
    str_replace_all("Alcohol use disorder", "AUD") %>%
    str_replace_all("Bipolar disorder", "BIP") %>%
    str_replace_all("Schizophrenia", "SCZ") %>%
    str_replace_all("Attention deficit hyperactivity disorder", "ADHD") %>%
    str_replace_all("Cardiovascular disease", "CVD") %>%
    str_replace_all("Coronary heart disease", "CHD") %>%
    str_replace_all("Ischemic heart disease", "IHD") %>%
    str_replace_all("Coronary atherosclerosis", "Coronary athero.") %>%
    str_replace_all("Atrial fibrillation", "AF") %>%
    str_replace_all("Myocardial infarction", "MI") %>%
    str_replace_all("Hypertension", "HTN") %>%
    str_replace_all("Antihypertensive use", "Anti-HTN use") %>%
    str_replace_all("Hypercholesterolemia", "Hyperchol.") %>%
    str_replace_all("Hyperlipidemia", "Hyperlipid.") %>%
    str_replace_all("Lipoprotein disorder", "Lipoprotein dis.") %>%
    str_replace_all("Type 1 diabetes", "T1D") %>%
    str_replace_all("Type 2 diabetes", "T2D") %>%
    str_replace_all("Rheumatoid arthritis", "RA") %>%
    str_replace_all("Autoimmune disease", "Autoimmune dis.") %>%
    str_replace_all("Non-thyroid autoimmune disease", "Non-thyroid autoimmune dis.") %>%
    str_replace_all("Strict non-thyroid autoimmune disease", "Strict non-thyroid autoimmune dis.") %>%
    str_replace_all("Systemic connective tissue disease", "Systemic connective tissue dis.") %>%
    str_replace_all("Other systemic connective disease", "Other systemic connective dis.") %>%
    str_replace_all("Diabetic retinopathy", "DM retinopathy") %>%
    str_replace_all("Diabetes with complications", "DM w/ comp.") %>%
    str_replace_all("Non-toxic thyroid disorder", "Non-toxic thyroid dis.") %>%
    str_replace_all("Papulosquamous disorder", "Papulosquamous dis.") %>%
    str_replace_all("disease", "dis.") %>%
    str_replace_all("disorder", "dis.") %>%
    str_replace_all("complications", "comp.") %>%
    str_replace_all("treatment", "tx") %>%
    str_replace_all("therapy", "tx") %>%
    str_squish()
}

# ----------------------------
# 3. Mortality EPOCH clock labels
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
  mutate(
    clock_block = factor(clock_block, levels = organ_levels),
    clock_modality = case_when(
      str_detect(str_to_lower(mortality_clock), "mri") ~ "MRI",
      str_detect(str_to_lower(mortality_clock), "proteomics") ~ "Proteomics",
      str_detect(str_to_lower(mortality_clock), "metabolomics") ~ "Metabolomics",
      TRUE ~ "Other"
    )
  ) %>%
  arrange(clock_block, clock_suborder) %>%
  mutate(clock_rank = row_number())

# ----------------------------
# 4. Disease endpoint chapter assignment
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

chapter_levels_all <- c(
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

# ----------------------------
# 5. Extract LDSC columns robustly
# ----------------------------
p_col <- get_first_col_name(dat, c("P", "p", "pval", "p_value"))
rg_col <- get_first_col_name(dat, c("gc_mean", "rg", "rg_mean", "rg_estimate", "rg_value"))

# IMPORTANT: gc_std is the SE column in your current LDSC table.
se_col <- get_first_col_name(dat, c("gc_std", "gc_se", "rg_se", "se", "SE", "gc_sd", "rg_sd"))

if (is.na(p_col)) stop("Could not find a P-value column.")
if (is.na(rg_col)) stop("Could not find an LDSC rg column.")
if (is.na(se_col)) stop("Could not find an LDSC SE column. Expected gc_std or similar.")

message("Using columns:")
message("  P-value column: ", p_col)
message("  rg column:      ", rg_col)
message("  SE column:      ", se_col)

P_vec <- safe_num(dat[[p_col]])
rg_vec <- safe_num(dat[[rg_col]])
se_vec <- safe_num(dat[[se_col]])

# ----------------------------
# 6. Disease endpoint data
# ----------------------------
disease <- dat %>%
  mutate(
    P_num = P_vec,
    rg = rg_vec,
    rg_se = se_vec
  ) %>%
  filter(
    analysis_group == "Disease_endpoint",
    target_source %in% c("FinnGen", "PGC")
  ) %>%
  mutate(
    endpoint_id = paste0(target_source, "::", target_id),
    endpoint_label_raw = target_display,
    endpoint_label = clean_text(target_display),
    disease_chapter = assign_disease_chapter(target_source, target_display),
    disease_chapter = factor(disease_chapter, levels = chapter_levels_all),
    valid_ldsc = !is.na(P_num) & !is.na(rg)
  ) %>%
  distinct(mortality_clock, endpoint_id, .keep_all = TRUE)

if (nrow(disease) == 0) {
  stop("No FinnGen/PGC disease-endpoint rows found in the input file.")
}

# ----------------------------
# 7. Bonferroni threshold
# ----------------------------
n_disease_for_correction <- 527
p_disease_bonf <- 0.05 / n_disease_for_correction
logp_disease_bonf <- -log10(p_disease_bonf)

message("Bonferroni threshold for DE-EPOCH LDSC:")
message("  P < 0.05 / ", n_disease_for_correction, " = ", signif(p_disease_bonf, 4))
message("  -log10(P) > ", signif(logp_disease_bonf, 4))

threshold_tbl <- tibble(
  correction = "Disease-endpoint Bonferroni",
  n_disease_for_correction = n_disease_for_correction,
  p_threshold = p_disease_bonf,
  neg_log10_p_threshold = logp_disease_bonf,
  p_col = p_col,
  rg_col = rg_col,
  se_col = se_col
)

readr::write_tsv(
  threshold_tbl,
  paste0(out_prefix, "_threshold.tsv")
)

# ----------------------------
# 8. Significant DE-EPOCH signals only
# ----------------------------
disease_sig <- disease %>%
  left_join(
    clock_key %>%
      select(
        mortality_clock,
        clock_block,
        clock_label,
        clock_suborder,
        clock_modality,
        clock_rank
      ),
    by = "mortality_clock"
  ) %>%
  mutate(
    neg_log10_P = -log10(P_num),
    sig_disease_bonf = valid_ldsc & P_num < p_disease_bonf,
    
    # 95% CI for LDSC rg.
    rg_ci_low = rg - 1.96 * rg_se,
    rg_ci_high = rg + 1.96 * rg_se
  ) %>%
  filter(sig_disease_bonf) %>%
  arrange(disease_chapter, endpoint_label, clock_rank, P_num)

if (nrow(disease_sig) == 0) {
  message("Debug summary:")
  print(
    disease %>%
      summarise(
        n_rows = n(),
        n_valid_p = sum(!is.na(P_num)),
        n_valid_rg = sum(!is.na(rg)),
        n_valid_se = sum(!is.na(rg_se)),
        min_p = min(P_num, na.rm = TRUE),
        n_passing_p = sum(P_num < p_disease_bonf, na.rm = TRUE)
      )
  )
  
  stop(
    "No DE-EPOCH LDSC genetic-correlation signals passed P < 0.05/527 = ",
    signif(p_disease_bonf, 4),
    ". Check whether the input file is the expected mortality-clock LDSC result table."
  )
}

message("Number of Bonferroni-significant DE-EPOCH pairs: ", nrow(disease_sig))
message("Number of unique significant disease endpoints: ", n_distinct(disease_sig$endpoint_id))
message("Number of mortality EPOCH clocks with at least one signal: ", n_distinct(disease_sig$mortality_clock))

if (nrow(disease_sig) != 52) {
  warning(
    "Expected 52 Bonferroni-significant DE-EPOCH pairs based on prior results, but found ",
    nrow(disease_sig),
    ". This may reflect a different input file or updated LDSC results."
  )
}

readr::write_tsv(
  disease_sig,
  paste0(out_prefix, "_Bonferroni_significant_DE_EPOCH_pairs.tsv")
)

# ----------------------------
# 9. Plotting table: exactly significant rows only
# ----------------------------
plot_df <- disease_sig %>%
  mutate(
    endpoint_short = shorten_endpoint(endpoint_label),
    
    row_label = paste0(endpoint_short, " | ", clock_label),
    
    clock_modality = factor(
      clock_modality,
      levels = c("MRI", "Proteomics", "Metabolomics", "Other")
    ),
    
    assoc_sign = if_else(rg >= 0, "Positive", "Negative"),
    p_label = format_pval(P_num),
    rg_label = sprintf("rg=%+.2f", rg)
  ) %>%
  arrange(disease_chapter, endpoint_short, clock_rank, P_num) %>%
  mutate(
    row_index = row_number(),
    y = rev(row_number())
  )

message("Rows plotted: ", nrow(plot_df))

# Disease-system group labels and separators.
group_df <- plot_df %>%
  group_by(disease_chapter) %>%
  summarise(
    y_min = min(y),
    y_max = max(y),
    y_mid = mean(range(y)),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(y_mid))

separator_df <- group_df %>%
  mutate(y_sep = y_min - 0.5) %>%
  filter(y_sep > 0.5)

# ----------------------------
# 10. Plot limits
# ----------------------------
x_min <- min(plot_df$rg_ci_low, na.rm = TRUE)
x_max <- max(plot_df$rg_ci_high, na.rm = TRUE)
x_pad <- 0.08 * (x_max - x_min)

x_lim <- c(
  min(-0.05, x_min - x_pad),
  max(0.05, x_max + x_pad)
)

# Space on the left for disease-system labels.
x_group_label <- x_lim[1] - 0.10 * diff(x_lim)

# ----------------------------
# 11. Colors and shapes
# ----------------------------
chapter_palette <- c(
  "PGC brain/psychiatric" = "#3B6EA5",
  "Neurologic" = "#4C78A8",
  "Psychiatric / sleep" = "#2F5D8C",
  "Eye" = "#74A9CF",
  "Ear / sensory" = "#9ECAE1",
  "Cardiovascular" = "#2E6B45",
  "Respiratory" = "#5AA05A",
  "Digestive / oral" = "#8E63A9",
  "Endocrine / immune" = "#9BAA4F",
  "Blood" = "#B2182B",
  "Cancer / benign neoplasm" = "#8B1A1A",
  "Skin / allergy" = "#C77CFF",
  "MSK / injury" = "#E39CB1",
  "GU / reproductive" = "#00008B",
  "Pregnancy" = "#E7298A",
  "Infectious" = "#A65E2E",
  "Other" = "#999999"
)

shape_values <- c(
  "MRI" = 21,
  "Proteomics" = 24,
  "Metabolomics" = 22,
  "Other" = 23
)

# ----------------------------
# 12. Build compact forest figure
# ----------------------------
p <- ggplot(plot_df, aes(x = rg, y = y)) +
  
  # Row gridlines, exactly one per significant signal.
  geom_hline(
    aes(yintercept = y),
    color = "grey91",
    linewidth = 0.22
  ) +
  
  # Disease-system separators.
  geom_hline(
    data = separator_df,
    aes(yintercept = y_sep),
    color = "grey70",
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  
  # Null line.
  geom_vline(
    xintercept = 0,
    color = "grey55",
    linewidth = 0.42
  ) +
  
  # Error bars from SE: rg +/- 1.96*SE.
  geom_segment(
    aes(
      x = rg_ci_low,
      xend = rg_ci_high,
      yend = y,
      color = disease_chapter
    ),
    linewidth = 0.70,
    alpha = 0.78,
    lineend = "round"
  ) +
  
  # Point estimates.
  geom_point(
    aes(
      color = disease_chapter,
      fill = disease_chapter,
      shape = clock_modality
    ),
    size = 3.6,
    stroke = 0.40,
    alpha = 0.98
  ) +
  
  # Disease-system labels on the left.
  geom_text(
    data = group_df,
    aes(
      x = x_group_label,
      y = y_mid,
      label = disease_chapter
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 3.0,
    fontface = "bold",
    color = "grey25"
  ) +
  
  scale_color_manual(values = chapter_palette, drop = TRUE) +
  scale_fill_manual(values = chapter_palette, drop = TRUE) +
  scale_shape_manual(values = shape_values, drop = TRUE) +
  
  scale_x_continuous(
    limits = x_lim,
    breaks = pretty_breaks(n = 6)
  ) +
  
  scale_y_continuous(
    breaks = plot_df$y,
    labels = plot_df$row_label,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    tag = "E",
    title = "Disease-endpoint genetic correlations with mortality EPOCH clocks",
    subtitle = paste0(
      "Bonferroni-significant DE-EPOCH LDSC signals only: P < 0.05/",
      n_disease_for_correction, " = ", signif(p_disease_bonf, 3),
      ". Points show LDSC genetic correlation; horizontal bars show rg +/- 1.96 x SE."
    ),
    x = "LDSC genetic correlation",
    y = NULL,
    color = "Disease system",
    fill = "Disease system",
    shape = "Mortality EPOCH modality",
    caption = paste0(
      "Each row represents one significant mortality EPOCH clock x disease endpoint pair. ",
      "Rows plotted: ", nrow(plot_df), "."
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 4),
      title.position = "top",
      title.hjust = 0
    ),
    fill = "none",
    shape = guide_legend(
      title.position = "top",
      title.hjust = 0
    )
  ) +
  
  theme_minimal(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.30),
    
    axis.text.y = element_text(
      size = 6.2,
      color = "grey20",
      lineheight = 0.90
    ),
    axis.text.x = element_text(
      size = 8,
      color = "grey20"
    ),
    axis.title.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(t = 6)
    ),
    
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.4),
    
    plot.tag = element_text(face = "bold", size = 22),
    plot.tag.position = c(0.006, 0.995),
    
    plot.title = element_text(face = "bold", size = 14.5),
    plot.subtitle = element_text(size = 8.7, color = "grey25", lineheight = 1.05),
    plot.caption = element_text(size = 7.2, color = "grey35", hjust = 0, lineheight = 1.05),
    
    plot.margin = margin(8, 14, 8, 95)
  )

print(p)

# ----------------------------
# 13. Save outputs
# ----------------------------
fig_height <- max(9.5, 0.19 * nrow(plot_df) + 2.4)
fig_width <- 14.2

ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 500,
  bg = "white"
)

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, ".svg"),
    plot = p,
    width = fig_width,
    height = fig_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

message("Done.")
message("Saved:")
message("  ", paste0(out_prefix, ".pdf"))
message("  ", paste0(out_prefix, ".png"))
if (requireNamespace("svglite", quietly = TRUE)) {
  message("  ", paste0(out_prefix, ".svg"))
}
message("  ", paste0(out_prefix, "_threshold.tsv"))
message("  ", paste0(out_prefix, "_Bonferroni_significant_DE_EPOCH_pairs.tsv"))