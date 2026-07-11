# ============================================================
# plot_clock_BAG_GCTA_h2_22_pairs.R
#
# Purpose:
#   Compare the 22 matched mortality clock-BAG pairs using
#   GCTA/HEreg SNP heritability estimates.
#
# Key improvements:
#   1. Mortality clock and BAG bars are visually distinct:
#        - Mortality clock: solid bar
#        - BAG: striped bar
#      If ggpattern is installed, diagonal striping is used.
#      If ggpattern is not installed, a robust manual vertical-stripe
#      fallback is used.
#
#   2. Scatter plot points use the SAME organ/system colors as
#      the bar plot.
#
#   3. Reproductive_female and Reproductive_male are matched using
#      full mortality_clock ID and full bag_id, not simplified organ labels.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})

use_ggpattern <- requireNamespace("ggpattern", quietly = TRUE)

# ============================================================
# 0. Paths
# ============================================================

root_candidates <- c(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock"
)

root_dir <- root_candidates[dir.exists(root_candidates)][1]

if (is.na(root_dir) || !dir.exists(root_dir)) {
  stop(
    "Cannot find WholeBodyClock root directory. Checked:\n",
    paste(root_candidates, collapse = "\n")
  )
}

mortality_h2_dir <- file.path(
  root_dir,
  "Figure",
  "heritability_3method_panel"
)

mortality_long_candidates <- c(
  file.path(mortality_h2_dir, "h2_3method_summary_long.tsv"),
  file.path(mortality_h2_dir, "h2_3method_long_format_estimates.tsv"),
  file.path(mortality_h2_dir, "h2_3method_long_format_direct_read_used_for_plot.tsv"),
  file.path(mortality_h2_dir, "h2_3method_summary_long_direct_read.tsv")
)

mortality_wide_candidates <- c(
  file.path(mortality_h2_dir, "h2_3method_summary_wide.tsv"),
  file.path(mortality_h2_dir, "h2_3method_df_wide_used_for_plot.tsv"),
  file.path(mortality_h2_dir, "h2_3method_df_wide_direct_read_used_for_plot.tsv")
)

bag_h2_candidates <- c(
  file.path(root_dir, "Result", "BAG_23_GCTA_h2_summary.tsv"),
  file.path(root_dir, "Result", "BAG_23_GCTA_h2_summary_wide.tsv")
)

out_dir <- file.path(
  root_dir,
  "Figure",
  "clock_BAG_GCTA_h2_pairs"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pick_existing <- function(candidates, label) {
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      "Cannot find ", label, ". Checked:\n",
      paste(candidates, collapse = "\n")
    )
  }
  hits[1]
}

mortality_long_file <- mortality_long_candidates[file.exists(mortality_long_candidates)][1]
mortality_wide_file <- mortality_wide_candidates[file.exists(mortality_wide_candidates)][1]
bag_h2_file <- pick_existing(bag_h2_candidates, "BAG h2 summary file")

message("============================================================")
message("Comparing mortality clock-BAG GCTA/HEreg h2 pairs")
message("============================================================")
message("Root directory       : ", root_dir)
message("Mortality long file : ", ifelse(is.na(mortality_long_file), "not found", mortality_long_file))
message("Mortality wide file : ", ifelse(is.na(mortality_wide_file), "not found", mortality_wide_file))
message("BAG h2 file         : ", bag_h2_file)
message("Output directory    : ", out_dir)
message("ggpattern available : ", use_ggpattern)
message("============================================================")

# ============================================================
# 1. Helper functions
# ============================================================

clean_colname <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

standardize_modality <- function(x) {
  x0 <- stringr::str_to_lower(as.character(x))
  
  case_when(
    str_detect(x0, "mri") ~ "MRI",
    str_detect(x0, "prot") ~ "Proteomics",
    str_detect(x0, "met") ~ "Metabolomics",
    TRUE ~ as.character(x)
  )
}

parse_modality_from_id <- function(x) {
  x0 <- stringr::str_to_lower(as.character(x))
  
  case_when(
    str_detect(x0, "_mri|mri") ~ "MRI",
    str_detect(x0, "proteomics|prot") ~ "Proteomics",
    str_detect(x0, "metabolomics|met") ~ "Metabolomics",
    TRUE ~ "Unknown"
  )
}

canonical_organ_key <- function(x) {
  x0 <- as.character(x)
  x0 <- stringr::str_replace_all(x0, "\\s+", "_")
  x0 <- stringr::str_replace_all(x0, "-", "_")
  x0 <- stringr::str_replace_all(x0, "/", "_")
  x0 <- stringr::str_to_lower(x0)
  
  case_when(
    str_detect(x0, "reproductive_female") ~ "reproductive_female",
    str_detect(x0, "reproductive_male") ~ "reproductive_male",
    
    str_detect(x0, "^brain|_brain") ~ "brain",
    str_detect(x0, "^adipose|_adipose") ~ "adipose",
    str_detect(x0, "^heart|_heart") ~ "heart",
    str_detect(x0, "^kidney|_kidney|^renal|_renal") ~ "renal",
    str_detect(x0, "^liver|_liver|^hepatic|_hepatic") ~ "hepatic",
    str_detect(x0, "^pancreas|_pancreas") ~ "pancreas",
    str_detect(x0, "^spleen|_spleen") ~ "spleen",
    str_detect(x0, "^endocrine|_endocrine") ~ "endocrine",
    str_detect(x0, "^eye|_eye") ~ "eye",
    str_detect(x0, "^immune|_immune") ~ "immune",
    str_detect(x0, "^pulmonary|_pulmonary|^lung|_lung") ~ "pulmonary",
    str_detect(x0, "^skin|_skin") ~ "skin",
    str_detect(x0, "^digestive|_digestive") ~ "digestive",
    str_detect(x0, "^metabolic|_metabolic") ~ "metabolic",
    
    str_detect(x0, "^reproductive|_reproductive") ~ "reproductive",
    
    TRUE ~ x0
  )
}

parse_organ_from_clock <- function(clock_id) {
  key <- canonical_organ_key(clock_id)
  
  case_when(
    key == "adipose" ~ "Adipose",
    key == "brain" ~ "Brain",
    key == "digestive" ~ "Digestive",
    key == "endocrine" ~ "Endocrine",
    key == "eye" ~ "Eye",
    key == "heart" ~ "Heart",
    key == "immune" ~ "Immune",
    key == "renal" ~ "Kidney/Renal",
    key == "hepatic" ~ "Liver/Hepatic",
    key == "pancreas" ~ "Pancreas",
    key == "pulmonary" ~ "Pulmonary",
    key == "reproductive_female" ~ "Reproductive female",
    key == "reproductive_male" ~ "Reproductive male",
    key == "skin" ~ "Skin",
    key == "spleen" ~ "Spleen",
    key == "metabolic" ~ "Metabolic",
    TRUE ~ stringr::str_to_title(as.character(clock_id))
  )
}

organ_display_from_key <- function(key) {
  case_when(
    key == "adipose" ~ "Adipose",
    key == "brain" ~ "Brain",
    key == "digestive" ~ "Digestive",
    key == "endocrine" ~ "Endocrine",
    key == "eye" ~ "Eye",
    key == "heart" ~ "Heart",
    key == "immune" ~ "Immune",
    key == "renal" ~ "Kidney/Renal",
    key == "hepatic" ~ "Liver/Hepatic",
    key == "pancreas" ~ "Pancreas",
    key == "pulmonary" ~ "Pulmonary",
    key == "reproductive_female" ~ "Reproductive\nfemale",
    key == "reproductive_male" ~ "Reproductive\nmale",
    key == "skin" ~ "Skin",
    key == "spleen" ~ "Spleen",
    key == "metabolic" ~ "Metabolic",
    TRUE ~ stringr::str_to_title(as.character(key))
  )
}

organ_short_label <- function(key) {
  case_when(
    key == "adipose" ~ "Adipose",
    key == "brain" ~ "Brain",
    key == "digestive" ~ "Digestive",
    key == "endocrine" ~ "Endocrine",
    key == "eye" ~ "Eye",
    key == "heart" ~ "Heart",
    key == "immune" ~ "Immune",
    key == "renal" ~ "Renal",
    key == "hepatic" ~ "Hepatic",
    key == "pancreas" ~ "Pancreas",
    key == "pulmonary" ~ "Pulmonary",
    key == "reproductive_female" ~ "Reproductive\nfemale",
    key == "reproductive_male" ~ "Reproductive\nmale",
    key == "skin" ~ "Skin",
    key == "spleen" ~ "Spleen",
    key == "metabolic" ~ "Metabolic",
    TRUE ~ stringr::str_to_title(as.character(key))
  )
}

format_pvalue <- function(p) {
  if (!is.finite(p)) {
    return("P = NA")
  } else if (p < 2.2e-16) {
    return("P < 2.2e-16")
  } else if (p < 0.001) {
    return(paste0("P = ", formatC(p, format = "e", digits = 2)))
  } else {
    return(paste0("P = ", signif(p, 3)))
  }
}

safe_cor_test <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x_ok <- x[ok]
  y_ok <- y[ok]
  
  if (length(x_ok) < 3) {
    return(list(estimate = NA_real_, p_value = NA_real_, n = length(x_ok)))
  }
  
  if (sd(x_ok, na.rm = TRUE) == 0 || sd(y_ok, na.rm = TRUE) == 0) {
    return(list(estimate = NA_real_, p_value = NA_real_, n = length(x_ok)))
  }
  
  ct <- suppressWarnings(
    cor.test(x_ok, y_ok, method = method, exact = FALSE)
  )
  
  list(
    estimate = unname(ct$estimate),
    p_value = ct$p.value,
    n = length(x_ok)
  )
}

# ============================================================
# 2. Color palette and plotting order
# ============================================================

organ_palette <- c(
  "Adipose" = "#E0AD3A",
  "Brain" = "#315C86",
  "Digestive" = "#9A7A43",
  "Endocrine" = "#E5C35A",
  "Eye" = "#577DB4",
  "Heart" = "#B85D48",
  "Immune" = "#6B966C",
  "Kidney/Renal" = "#4E8C7E",
  "Liver/Hepatic" = "#D09B24",
  "Pancreas" = "#7866AA",
  "Pulmonary" = "#5A9DB1",
  "Reproductive\nfemale" = "#B984AA",
  "Reproductive\nmale" = "#8675B7",
  "Skin" = "#CF823F",
  "Spleen" = "#35678F",
  "Metabolic" = "#A58A4F"
)

modality_order <- c("Proteomics", "MRI", "Metabolomics")

organ_order_by_modality <- list(
  "Proteomics" = c(
    "brain",
    "endocrine",
    "eye",
    "heart",
    "hepatic",
    "immune",
    "pulmonary",
    "renal",
    "reproductive_female",
    "reproductive_male",
    "skin"
  ),
  "MRI" = c(
    "adipose",
    "brain",
    "heart",
    "renal",
    "hepatic",
    "pancreas",
    "spleen"
  ),
  "Metabolomics" = c(
    "digestive",
    "endocrine",
    "hepatic",
    "immune",
    "metabolic"
  )
)

# Filled shapes, so fill = organ color works correctly.
modality_shapes <- c(
  "Proteomics" = 24,
  "MRI" = 21,
  "Metabolomics" = 22
)

theme_clean <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(
        face = "bold",
        color = "#1F2A44"
      ),
      plot.subtitle = element_text(
        color = "#4A4A4A"
      ),
      axis.title = element_text(
        face = "bold",
        color = "#1F2A44"
      ),
      axis.text = element_text(
        color = "#1F2A44"
      ),
      panel.grid.major.y = element_line(
        color = "#E7DED1",
        linewidth = 0.45
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(
        face = "bold"
      )
    )
}

# ============================================================
# 3. Read mortality clock raw GCTA/HEreg h2 estimates
# ============================================================

read_mortality_h2 <- function(long_file, wide_file) {
  
  if (!is.na(long_file) && file.exists(long_file)) {
    
    x <- readr::read_tsv(long_file, show_col_types = FALSE, progress = FALSE)
    colnames(x) <- clean_colname(colnames(x))
    
    if (!"mortality_clock_chr" %in% colnames(x)) {
      if ("mortality_clock" %in% colnames(x)) {
        x <- x %>% rename(mortality_clock_chr = mortality_clock)
      } else {
        stop("Long mortality file has no mortality_clock_chr or mortality_clock column: ", long_file)
      }
    }
    
    if (!"method" %in% colnames(x)) {
      stop("Long mortality file has no method column: ", long_file)
    }
    
    if (!"h2_mean" %in% colnames(x)) {
      if ("h2" %in% colnames(x)) {
        x <- x %>% mutate(h2_mean = h2)
      } else if ("h2_estimate" %in% colnames(x)) {
        x <- x %>% mutate(h2_mean = h2_estimate)
      } else {
        stop("Long mortality file has no h2_mean/h2/h2_estimate column: ", long_file)
      }
    }
    
    if (!"h2_se" %in% colnames(x)) {
      if ("se" %in% colnames(x)) {
        x <- x %>% mutate(h2_se = se)
      } else {
        x <- x %>% mutate(h2_se = NA_real_)
      }
    }
    
    if (!"modality" %in% colnames(x)) {
      x <- x %>% mutate(modality = parse_modality_from_id(mortality_clock_chr))
    }
    
    if (!"organ" %in% colnames(x)) {
      x <- x %>% mutate(organ = parse_organ_from_clock(mortality_clock_chr))
    }
    
    out <- x %>%
      mutate(
        method_clean = stringr::str_to_lower(as.character(method)),
        method_clean = stringr::str_replace_all(method_clean, "[^a-z0-9]+", "_")
      ) %>%
      filter(
        stringr::str_detect(method_clean, "raw_gcta_hereg|gcta|hereg|he_reg|raw")
      ) %>%
      transmute(
        mortality_clock = as.character(mortality_clock_chr),
        modality = standardize_modality(modality),
        organ_raw = as.character(organ),
        
        # Critical fix: use full mortality clock ID.
        organ_key = canonical_organ_key(mortality_clock),
        
        organ_system = organ_display_from_key(organ_key),
        pair_id = paste(modality, organ_key, sep = "__"),
        
        clock_h2 = safe_numeric(h2_mean),
        clock_se = safe_numeric(h2_se),
        mortality_source_file = long_file
      ) %>%
      filter(is.finite(clock_h2)) %>%
      distinct(pair_id, .keep_all = TRUE)
    
    return(out)
  }
  
  if (!is.na(wide_file) && file.exists(wide_file)) {
    
    x <- readr::read_tsv(wide_file, show_col_types = FALSE, progress = FALSE)
    colnames(x) <- clean_colname(colnames(x))
    
    if (!"mortality_clock_chr" %in% colnames(x)) {
      if ("mortality_clock" %in% colnames(x)) {
        x <- x %>% rename(mortality_clock_chr = mortality_clock)
      } else {
        stop("Wide mortality file has no mortality_clock_chr or mortality_clock column: ", wide_file)
      }
    }
    
    if (!"h2_raw_gcta_hereg" %in% colnames(x)) {
      stop("Wide mortality file has no h2_raw_gcta_hereg column: ", wide_file)
    }
    
    if (!"modality" %in% colnames(x)) {
      x <- x %>% mutate(modality = parse_modality_from_id(mortality_clock_chr))
    }
    
    if (!"organ" %in% colnames(x)) {
      x <- x %>% mutate(organ = parse_organ_from_clock(mortality_clock_chr))
    }
    
    out <- x %>%
      transmute(
        mortality_clock = as.character(mortality_clock_chr),
        modality = standardize_modality(modality),
        organ_raw = as.character(organ),
        
        # Critical fix: use full mortality clock ID.
        organ_key = canonical_organ_key(mortality_clock),
        
        organ_system = organ_display_from_key(organ_key),
        pair_id = paste(modality, organ_key, sep = "__"),
        
        clock_h2 = safe_numeric(h2_raw_gcta_hereg),
        clock_se = NA_real_,
        mortality_source_file = wide_file
      ) %>%
      filter(is.finite(clock_h2)) %>%
      distinct(pair_id, .keep_all = TRUE)
    
    return(out)
  }
  
  stop(
    "No mortality h2 summary file found.\n",
    "Long candidates:\n", paste(mortality_long_candidates, collapse = "\n"), "\n",
    "Wide candidates:\n", paste(mortality_wide_candidates, collapse = "\n")
  )
}

mortality_h2 <- read_mortality_h2(
  long_file = mortality_long_file,
  wide_file = mortality_wide_file
)

message("Mortality-clock raw GCTA/HEreg h2 rows:")
print(mortality_h2 %>% count(modality, name = "n"))

message("Mortality reproductive pair IDs:")
print(
  mortality_h2 %>%
    filter(str_detect(pair_id, "reproductive")) %>%
    select(pair_id, mortality_clock, organ_raw, organ_key, clock_h2)
)

# ============================================================
# 4. Read BAG GCTA/HEreg h2 estimates
# ============================================================

bag_raw <- readr::read_tsv(
  bag_h2_file,
  show_col_types = FALSE,
  progress = FALSE
)

colnames(bag_raw) <- clean_colname(colnames(bag_raw))

if (!"bag_id" %in% colnames(bag_raw)) {
  stop("BAG h2 file missing bag_id column: ", bag_h2_file)
}

if (!"h2" %in% colnames(bag_raw)) {
  if ("h2_mean" %in% colnames(bag_raw)) {
    bag_raw <- bag_raw %>% mutate(h2 = h2_mean)
  } else {
    stop("BAG h2 file missing h2 or h2_mean column: ", bag_h2_file)
  }
}

if (!"se" %in% colnames(bag_raw)) {
  if ("h2_se" %in% colnames(bag_raw)) {
    bag_raw <- bag_raw %>% mutate(se = h2_se)
  } else {
    bag_raw <- bag_raw %>% mutate(se = NA_real_)
  }
}

if (!"modality" %in% colnames(bag_raw)) {
  bag_raw <- bag_raw %>% mutate(modality = parse_modality_from_id(bag_id))
}

if (!"organ" %in% colnames(bag_raw)) {
  bag_raw <- bag_raw %>% mutate(organ = parse_organ_from_clock(bag_id))
}

if (!"parse_status" %in% colnames(bag_raw)) {
  bag_raw <- bag_raw %>% mutate(parse_status = "ok")
}

bag_h2 <- bag_raw %>%
  mutate(
    modality = standardize_modality(modality),
    
    # Critical fix: use full BAG ID.
    organ_key = canonical_organ_key(bag_id),
    
    organ_system = organ_display_from_key(organ_key),
    pair_id = paste(modality, organ_key, sep = "__")
  ) %>%
  filter(
    parse_status == "ok",
    is.finite(safe_numeric(h2))
  ) %>%
  transmute(
    bag_id = as.character(bag_id),
    modality,
    organ_raw = as.character(organ),
    organ_key,
    organ_system,
    pair_id,
    bag_h2 = safe_numeric(h2),
    bag_se = safe_numeric(se),
    bag_pvalue = if ("pvalue" %in% colnames(bag_raw)) safe_numeric(pvalue) else NA_real_,
    bag_n = if ("n" %in% colnames(bag_raw)) safe_numeric(n) else NA_real_,
    bag_source_file = bag_h2_file
  ) %>%
  distinct(pair_id, .keep_all = TRUE)

message("BAG GCTA/HEreg h2 rows:")
print(bag_h2 %>% count(modality, name = "n"))

message("BAG reproductive pair IDs:")
print(
  bag_h2 %>%
    filter(str_detect(pair_id, "reproductive")) %>%
    select(pair_id, bag_id, organ_raw, organ_key, bag_h2)
)

# ============================================================
# 5. Pair mortality clocks and BAGs
# ============================================================

clock_bag_pairs <- mortality_h2 %>%
  inner_join(
    bag_h2 %>%
      select(
        pair_id,
        bag_id,
        bag_h2,
        bag_se,
        bag_pvalue,
        bag_n,
        bag_source_file
      ),
    by = "pair_id"
  ) %>%
  mutate(
    modality = factor(modality, levels = modality_order),
    organ_system = organ_display_from_key(organ_key),
    organ_label_short = organ_short_label(organ_key),
    h2_difference_clock_minus_BAG = clock_h2 - bag_h2,
    h2_ratio_clock_over_BAG = clock_h2 / bag_h2
  )

pair_order_tbl <- purrr::map_dfr(modality_order, function(mod) {
  tibble(
    modality = mod,
    organ_key = organ_order_by_modality[[mod]],
    within_modality_order = seq_along(organ_order_by_modality[[mod]])
  )
})

clock_bag_pairs <- clock_bag_pairs %>%
  left_join(
    pair_order_tbl,
    by = c("modality", "organ_key")
  ) %>%
  arrange(
    modality,
    within_modality_order,
    organ_key
  ) %>%
  mutate(
    pair_index = row_number(),
    pair_label = organ_short_label(organ_key),
    organ_system = factor(
      organ_system,
      levels = names(organ_palette)
    )
  )

unpaired_mortality <- mortality_h2 %>%
  anti_join(bag_h2, by = "pair_id") %>%
  arrange(modality, organ_key)

unpaired_bags <- bag_h2 %>%
  anti_join(mortality_h2, by = "pair_id") %>%
  arrange(modality, organ_key)

if (nrow(clock_bag_pairs) != 22) {
  warning(
    "Expected 22 paired clock-BAG rows, but found ",
    nrow(clock_bag_pairs),
    ". Check the unpaired output tables."
  )
}

message("Paired rows: ", nrow(clock_bag_pairs))
message("Unpaired mortality clocks: ", nrow(unpaired_mortality))
message("Unpaired BAGs: ", nrow(unpaired_bags))

message("Matched reproductive pairs:")
print(
  clock_bag_pairs %>%
    filter(str_detect(pair_id, "reproductive")) %>%
    select(
      pair_id,
      modality,
      organ_key,
      organ_system,
      mortality_clock,
      bag_id,
      clock_h2,
      bag_h2
    )
)

# ============================================================
# 6. Save paired and unpaired tables
# ============================================================

paired_file <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs.tsv")
unpaired_mortality_file <- file.path(out_dir, "clock_BAG_GCTA_h2_unpaired_mortality_clocks.tsv")
unpaired_bags_file <- file.path(out_dir, "clock_BAG_GCTA_h2_unpaired_BAGs.tsv")

readr::write_tsv(clock_bag_pairs, paired_file)
readr::write_tsv(unpaired_mortality, unpaired_mortality_file)
readr::write_tsv(unpaired_bags, unpaired_bags_file)

# ============================================================
# 7. Prepare long-format plotting table
# ============================================================

plot_long <- bind_rows(
  clock_bag_pairs %>%
    transmute(
      pair_id,
      pair_index,
      pair_label,
      mortality_clock,
      bag_id,
      modality,
      organ_key,
      organ_system,
      source = "Mortality clock",
      h2 = clock_h2,
      se = clock_se
    ),
  clock_bag_pairs %>%
    transmute(
      pair_id,
      pair_index,
      pair_label,
      mortality_clock,
      bag_id,
      modality,
      organ_key,
      organ_system,
      source = "BAG",
      h2 = bag_h2,
      se = bag_se
    )
) %>%
  mutate(
    source = factor(source, levels = c("Mortality clock", "BAG")),
    modality = factor(as.character(modality), levels = modality_order),
    organ_system = factor(as.character(organ_system), levels = names(organ_palette))
  ) %>%
  arrange(pair_index, source)

plot_long_file <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_long_for_plot.tsv")
readr::write_tsv(plot_long, plot_long_file)

modality_bounds <- clock_bag_pairs %>%
  group_by(modality) %>%
  summarise(
    x_min = min(pair_index),
    x_max = max(pair_index),
    x_mid = mean(range(pair_index)),
    .groups = "drop"
  )

separator_tbl <- modality_bounds %>%
  arrange(x_min) %>%
  mutate(separator_x = lag(x_min) - 0.5) %>%
  filter(is.finite(separator_x))

# ============================================================
# 8. Pairwise statistics
# ============================================================

pearson <- safe_cor_test(
  clock_bag_pairs$bag_h2,
  clock_bag_pairs$clock_h2,
  method = "pearson"
)

spearman <- safe_cor_test(
  clock_bag_pairs$bag_h2,
  clock_bag_pairs$clock_h2,
  method = "spearman"
)

paired_t <- suppressWarnings(
  t.test(
    clock_bag_pairs$clock_h2,
    clock_bag_pairs$bag_h2,
    paired = TRUE
  )
)

paired_w <- suppressWarnings(
  wilcox.test(
    clock_bag_pairs$clock_h2,
    clock_bag_pairs$bag_h2,
    paired = TRUE,
    exact = FALSE
  )
)

comparison_stats <- tibble(
  n_pairs = nrow(clock_bag_pairs),
  pearson_r = pearson$estimate,
  pearson_p = pearson$p_value,
  spearman_rho = spearman$estimate,
  spearman_p = spearman$p_value,
  mean_clock_h2 = mean(clock_bag_pairs$clock_h2, na.rm = TRUE),
  mean_BAG_h2 = mean(clock_bag_pairs$bag_h2, na.rm = TRUE),
  mean_difference_clock_minus_BAG = mean(clock_bag_pairs$h2_difference_clock_minus_BAG, na.rm = TRUE),
  median_difference_clock_minus_BAG = median(clock_bag_pairs$h2_difference_clock_minus_BAG, na.rm = TRUE),
  paired_t_p = paired_t$p.value,
  paired_wilcoxon_p = paired_w$p.value
)

stats_file <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_statistics.tsv")
readr::write_tsv(comparison_stats, stats_file)

scatter_label <- paste0(
  "N = ", comparison_stats$n_pairs,
  "\nPearson r = ", sprintf("%.2f", comparison_stats$pearson_r),
  "\n", format_pvalue(comparison_stats$pearson_p),
  "\nMean diff = ", sprintf("%.3f", comparison_stats$mean_difference_clock_minus_BAG)
)

# ============================================================
# 9. Panel A: paired bar plot
# ============================================================

pd <- position_dodge(width = 0.76)

h2_ymax <- max(
  plot_long$h2 + ifelse(is.finite(plot_long$se), plot_long$se, 0),
  na.rm = TRUE
)

h2_ymax <- h2_ymax * 1.22

# ------------------------------------------------------------
# Preferred version: ggpattern diagonal stripes for BAG bars.
# ------------------------------------------------------------

if (use_ggpattern) {
  
  p_bar <- ggplot(
    plot_long,
    aes(
      x = pair_index,
      y = h2,
      fill = organ_system,
      pattern = source
    )
  ) +
    geom_vline(
      data = separator_tbl,
      aes(xintercept = separator_x),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "grey35",
      linewidth = 0.55
    ) +
    ggpattern::geom_col_pattern(
      position = pd,
      width = 0.62,
      color = "black",
      linewidth = 0.28,
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_angle = 45,
      pattern_density = 0.28,
      pattern_spacing = 0.035,
      pattern_alpha = 0.45,
      pattern_key_scale_factor = 0.75
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(h2 - se, 0),
        ymax = h2 + se
      ),
      position = pd,
      width = 0.17,
      linewidth = 0.35,
      color = "black",
      na.rm = TRUE
    ) +
    geom_text(
      data = modality_bounds,
      aes(
        x = x_mid,
        y = h2_ymax * 0.985,
        label = modality
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.7,
      color = "#1F2A44"
    ) +
    scale_fill_manual(
      values = organ_palette,
      drop = FALSE,
      name = "Organ/system"
    ) +
    ggpattern::scale_pattern_manual(
      values = c(
        "Mortality clock" = "none",
        "BAG" = "stripe"
      ),
      name = "Biomarker"
    ) +
    scale_x_continuous(
      breaks = clock_bag_pairs$pair_index,
      labels = clock_bag_pairs$pair_label,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, h2_ymax),
      labels = number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.00, 0.02))
    ) +
    labs(
      title = expression("GCTA/HEreg SNP heritability (" * h^2 * "): mortality clock vs BAG"),
      subtitle = paste0(
        "Matched by modality and organ/system; ",
        nrow(clock_bag_pairs),
        " mortality clock-BAG pairs"
      ),
      x = NULL,
      y = expression(h^2 * " estimate " %+-% " SE")
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(pattern = "none"),
        ncol = 1
      ),
      pattern = guide_legend(
        override.aes = list(fill = "grey70")
      )
    )
  
} else {
  
  # ----------------------------------------------------------
  # Fallback version: manual vertical stripes for BAG bars.
  # This avoids requiring the ggpattern package.
  # ----------------------------------------------------------
  
  bar_width <- 0.34
  dodge_offset <- 0.20
  stripe_spacing <- 0.055
  
  plot_long_rect <- plot_long %>%
    mutate(
      x_pos = pair_index + if_else(source == "Mortality clock", -dodge_offset, dodge_offset),
      bar_left = x_pos - bar_width / 2,
      bar_right = x_pos + bar_width / 2
    )
  
  stripe_df <- plot_long_rect %>%
    filter(source == "BAG") %>%
    rowwise() %>%
    do({
      d <- .
      xs <- seq(
        from = d$bar_left + stripe_spacing / 2,
        to = d$bar_right - stripe_spacing / 2,
        by = stripe_spacing
      )
      tibble(
        pair_index = d$pair_index,
        x = xs,
        y = 0,
        xend = xs,
        yend = d$h2
      )
    }) %>%
    ungroup()
  
  p_bar <- ggplot() +
    geom_vline(
      data = separator_tbl,
      aes(xintercept = separator_x),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "grey35",
      linewidth = 0.55
    ) +
    geom_rect(
      data = plot_long_rect,
      aes(
        xmin = bar_left,
        xmax = bar_right,
        ymin = 0,
        ymax = h2,
        fill = organ_system,
        alpha = source
      ),
      color = "black",
      linewidth = 0.28
    ) +
    geom_segment(
      data = stripe_df,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.18,
      alpha = 0.65
    ) +
    geom_errorbar(
      data = plot_long_rect,
      aes(
        x = x_pos,
        ymin = pmax(h2 - se, 0),
        ymax = h2 + se
      ),
      width = 0.12,
      linewidth = 0.35,
      color = "black",
      na.rm = TRUE
    ) +
    geom_text(
      data = modality_bounds,
      aes(
        x = x_mid,
        y = h2_ymax * 0.985,
        label = modality
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.7,
      color = "#1F2A44"
    ) +
    scale_fill_manual(
      values = organ_palette,
      drop = FALSE,
      name = "Organ/system"
    ) +
    scale_alpha_manual(
      values = c(
        "Mortality clock" = 0.95,
        "BAG" = 0.55
      ),
      name = "Biomarker"
    ) +
    scale_x_continuous(
      breaks = clock_bag_pairs$pair_index,
      labels = clock_bag_pairs$pair_label,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = c(0, h2_ymax),
      labels = number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.00, 0.02))
    ) +
    labs(
      title = expression("GCTA/HEreg SNP heritability (" * h^2 * "): mortality clock vs BAG"),
      subtitle = paste0(
        "Matched by modality and organ/system; ",
        nrow(clock_bag_pairs),
        " mortality clock-BAG pairs"
      ),
      x = NULL,
      y = expression(h^2 * " estimate " %+-% " SE")
    )
}

p_bar <- p_bar +
  theme_clean(base_size = 11) +
  theme(
    plot.title = element_text(size = 15),
    plot.subtitle = element_text(size = 10),
    axis.title.y = element_text(size = 10.5),
    axis.text.y = element_text(size = 8.5),
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 7.4
    ),
    legend.position = "right",
    plot.margin = margin(8, 10, 8, 8)
  )

# ============================================================
# 10. Panel B: paired scatter plot
# ============================================================

scatter_x_rng <- range(clock_bag_pairs$bag_h2, na.rm = TRUE)
scatter_y_rng <- range(clock_bag_pairs$clock_h2, na.rm = TRUE)

scatter_min <- min(scatter_x_rng[1], scatter_y_rng[1], 0)
scatter_max <- max(scatter_x_rng[2], scatter_y_rng[2])
scatter_pad <- 0.06 * (scatter_max - scatter_min)

if (!is.finite(scatter_pad) || scatter_pad == 0) {
  scatter_pad <- 0.05
}

scatter_xlim <- c(scatter_min, scatter_max + scatter_pad)
scatter_ylim <- c(scatter_min, scatter_max + scatter_pad)

scatter_annot <- tibble(
  x = scatter_xlim[1] + 0.04 * diff(scatter_xlim),
  y = scatter_ylim[2] - 0.05 * diff(scatter_ylim),
  label = scatter_label
)

p_scatter <- ggplot(
  clock_bag_pairs,
  aes(
    x = bag_h2,
    y = clock_h2,
    fill = organ_system,
    shape = modality
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.55
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "#1F2A44",
    linewidth = 0.75
  ) +
  geom_point(
    size = 3.4,
    color = "black",
    stroke = 0.45,
    alpha = 0.96
  ) +
  geom_label(
    data = scatter_annot,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.2,
    fontface = "bold",
    color = "#1F2A44",
    fill = "white",
    label.size = 0.25,
    label.padding = unit(0.18, "lines")
  ) +
  scale_fill_manual(
    values = organ_palette,
    drop = FALSE,
    name = "Organ/system"
  ) +
  scale_shape_manual(
    values = modality_shapes,
    drop = FALSE,
    name = "Modality"
  ) +
  coord_cartesian(
    xlim = scatter_xlim,
    ylim = scatter_ylim,
    clip = "off"
  ) +
  scale_x_continuous(
    labels = number_format(accuracy = 0.01)
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    title = expression("Pairwise comparison of " * h^2 * " estimates"),
    x = expression("BAG GCTA/HEreg " * h^2),
    y = expression("Mortality clock GCTA/HEreg " * h^2)
  ) +
  theme_clean(base_size = 11) +
  theme(
    plot.title = element_text(size = 13.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5),
    legend.position = "right",
    plot.margin = margin(8, 10, 8, 8)
  )

# ============================================================
# 11. Combined figure
# ============================================================

combined_plot <- p_bar / p_scatter +
  plot_layout(
    heights = c(1.55, 1.00),
    guides = "collect"
  ) +
  plot_annotation(
    title = expression("Matched comparison of mortality clock and BAG SNP heritability (" * h^2 * ")"),
    subtitle = "GCTA/HEreg estimates paired by organ/system and modality",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        color = "#1F2A44"
      ),
      plot.subtitle = element_text(
        size = 11,
        color = "#4A4A4A"
      ),
      plot.tag = element_text(
        face = "bold",
        size = 20,
        color = "black"
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      legend.position = "right"
    )
  )

# ============================================================
# 12. Save figures
# ============================================================

bar_pdf <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_paired_barplot.pdf")
bar_png <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_paired_barplot.png")

scatter_pdf <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_scatter.pdf")
scatter_png <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_scatter.png")

combined_pdf <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_comparison_figure.pdf")
combined_png <- file.path(out_dir, "clock_BAG_GCTA_h2_22_pairs_comparison_figure.png")

ggsave(
  bar_pdf,
  p_bar,
  width = 14.5,
  height = 5.8,
  device = cairo_pdf
)

ggsave(
  bar_png,
  p_bar,
  width = 14.5,
  height = 5.8,
  dpi = 500
)

ggsave(
  scatter_pdf,
  p_scatter,
  width = 7.0,
  height = 5.8,
  device = cairo_pdf
)

ggsave(
  scatter_png,
  p_scatter,
  width = 7.0,
  height = 5.8,
  dpi = 500
)

ggsave(
  combined_pdf,
  combined_plot,
  width = 14.5,
  height = 10.5,
  device = cairo_pdf
)

ggsave(
  combined_png,
  combined_plot,
  width = 14.5,
  height = 10.5,
  dpi = 500
)

# ============================================================
# 13. Completion message
# ============================================================

message("============================================================")
message("Done.")
message("Saved output files:")
message("  ", paired_file)
message("  ", plot_long_file)
message("  ", stats_file)
message("  ", unpaired_mortality_file)
message("  ", unpaired_bags_file)
message("  ", bar_pdf)
message("  ", bar_png)
message("  ", scatter_pdf)
message("  ", scatter_png)
message("  ", combined_pdf)
message("  ", combined_png)
message("============================================================")

message("Paired clock-BAG h2 table preview:")
print(
  clock_bag_pairs %>%
    select(
      modality,
      organ_system,
      pair_id,
      mortality_clock,
      bag_id,
      clock_h2,
      clock_se,
      bag_h2,
      bag_se,
      h2_difference_clock_minus_BAG
    )
)