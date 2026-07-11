# ============================================================
# Fig_heritability_bar_plot_3_methods.R
#
# Three-panel figure:
#   A. h2 estimates across mortality clocks and methods
#   B. Pairwise agreement of h2 estimates across methods
#   C. Bland-Altman comparison of h2 estimates
#
# Methods:
#   1. GCTA/HEreg raw genotype
#   2. SBayesS summary
#   3. LDSC summary
#
# Inputs:
#   GCTA/HEreg:
#     ~/Reproducibile_paper/WholeBodyClock/mortality_clock/GCTA_h2/*/*.hsq
#     ~/Reproducibile_paper/WholeBodyClock/mortality_clock/GCTA_h2/*/*.HEreg
#
#   LDSC:
#     ~/Reproducibile_paper/WholeBodyClock/Result/LDSC_h2_intercept_mortality_clocks.tsv
#
#   SBayesS:
#     ~/Reproducibile_paper/WholeBodyClock/Result/GCTB_SBayesS_parameters_mortality_clocks.tsv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})

# ============================================================
# 0. Paths
# ============================================================

root_dir <- path.expand("/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock")

gcta_hereg_dir <- file.path(
  root_dir,
  "mortality_clock",
  "GCTA_h2"
)

ldsc_file <- file.path(
  root_dir,
  "Result",
  "LDSC_h2_intercept_mortality_clocks.tsv"
)

sbayess_file <- file.path(
  root_dir,
  "Result",
  "GCTB_SBayesS_parameters_mortality_clocks.tsv"
)

output_dir <- file.path(
  root_dir,
  "Figure",
  "heritability_3method_panel"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(gcta_hereg_dir)) {
  stop("Cannot find GCTA/HEreg directory: ", gcta_hereg_dir)
}

if (!file.exists(ldsc_file)) {
  stop("Cannot find LDSC file: ", ldsc_file)
}

if (!file.exists(sbayess_file)) {
  stop("Cannot find SBayesS file: ", sbayess_file)
}

message("============================================================")
message("Generating 3-panel h2 consistency figure")
message("============================================================")
message("GCTA/HEreg directory: ", gcta_hereg_dir)
message("LDSC file          : ", ldsc_file)
message("SBayesS file       : ", sbayess_file)
message("Output directory   : ", output_dir)
message("============================================================")

# ============================================================
# 1. Helper functions
# ============================================================

clean_colname <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    tolower()
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
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
    return(
      list(
        estimate = NA_real_,
        p_value = NA_real_,
        n = length(x_ok)
      )
    )
  }
  
  if (sd(x_ok, na.rm = TRUE) == 0 || sd(y_ok, na.rm = TRUE) == 0) {
    return(
      list(
        estimate = NA_real_,
        p_value = NA_real_,
        n = length(x_ok)
      )
    )
  }
  
  ct <- suppressWarnings(
    cor.test(
      x_ok,
      y_ok,
      method = method,
      exact = FALSE
    )
  )
  
  list(
    estimate = unname(ct$estimate),
    p_value = ct$p.value,
    n = length(x_ok)
  )
}

lin_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 3) {
    return(NA_real_)
  }
  
  mx <- mean(x)
  my <- mean(y)
  
  vx <- var(x)
  vy <- var(y)
  
  sxy <- cov(x, y)
  
  denom <- vx + vy + (mx - my)^2
  
  if (!is.finite(denom) || denom == 0) {
    return(NA_real_)
  }
  
  (2 * sxy) / denom
}

parse_clock_organ <- function(clock_id) {
  x <- str_to_lower(as.character(clock_id))
  
  case_when(
    str_detect(x, "^brain") ~ "Brain",
    str_detect(x, "^eye") ~ "Eye",
    str_detect(x, "^pulmonary|^lung") ~ "Pulmonary",
    str_detect(x, "^heart") ~ "Heart",
    str_detect(x, "^hepatic|^liver") ~ "Hepatic",
    str_detect(x, "^renal|^kidney") ~ "Renal",
    str_detect(x, "^pancreas") ~ "Pancreas",
    str_detect(x, "^spleen") ~ "Spleen",
    str_detect(x, "^immune") ~ "Immune",
    str_detect(x, "^endocrine") ~ "Endocrine",
    str_detect(x, "^digestive") ~ "Digestive",
    str_detect(x, "^metabolic") ~ "Metabolic",
    str_detect(x, "^adipose") ~ "Adipose",
    str_detect(x, "^skin") ~ "Skin",
    str_detect(x, "^reproductive") ~ "Reproductive",
    TRUE ~ "Other"
  )
}

parse_clock_display_organ <- function(clock_id) {
  x <- str_to_lower(as.character(clock_id))
  
  case_when(
    str_detect(x, "^reproductive_female") ~ "Reproductive female",
    str_detect(x, "^reproductive_male") ~ "Reproductive male",
    str_detect(x, "^kidney") ~ "Renal",
    str_detect(x, "^liver") ~ "Hepatic",
    TRUE ~ parse_clock_organ(clock_id)
  )
}

parse_clock_modality <- function(clock_id) {
  x <- str_to_lower(as.character(clock_id))
  
  case_when(
    str_detect(x, "_mri$|_mri_") ~ "MRI",
    str_detect(x, "proteomics") ~ "Proteomics",
    str_detect(x, "metabolomics") ~ "Metabolomics",
    TRUE ~ "Other"
  )
}

pretty_clock_label_one_line <- function(clock_id) {
  organ <- parse_clock_display_organ(clock_id)
  modality <- parse_clock_modality(clock_id)
  
  modality_short <- case_when(
    modality == "MRI" ~ "MRI",
    modality == "Proteomics" ~ "Prot",
    modality == "Metabolomics" ~ "Met",
    TRUE ~ modality
  )
  
  paste(organ, modality_short)
}

# ============================================================
# 2. Method labels, colors, and plotting controls
# ============================================================

method_order <- c(
  "raw_gcta_hereg",
  "sbayess",
  "ldsc"
)

method_labels <- c(
  raw_gcta_hereg = "GCTA/HEreg raw genotype",
  sbayess = "SBayesS summary",
  ldsc = "LDSC summary"
)

method_colors <- c(
  "GCTA/HEreg raw genotype" = "#376B9B",
  "SBayesS summary" = "#E5A11A",
  "LDSC summary" = "#5C9465"
)

modality_order <- c("MRI", "Proteomics", "Metabolomics", "Other")

modality_shapes <- c(
  MRI = 16,
  Proteomics = 17,
  Metabolomics = 15,
  Other = 4
)

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(face = "bold", color = "#1F2A44"),
      plot.subtitle = element_text(color = "#4A4A4A"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", color = "black"),
      panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1)
    )
}

# ============================================================
# 3. Read GCTA / HEreg h2 estimates
# ============================================================

read_hsq_one <- function(file, clock_folder) {
  dat <- readr::read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  colnames(dat) <- clean_colname(colnames(dat))
  
  if (!all(c("source", "variance", "se") %in% colnames(dat))) {
    return(
      tibble(
        clock_folder = clock_folder,
        mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
        modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
        organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
        method = "raw_gcta_hereg",
        h2_mean = NA_real_,
        h2_se = NA_real_,
        p_value = NA_real_,
        n = NA_real_,
        raw_method_source = "GCTA_REML_hsq",
        source_file = file,
        read_status = "missing_required_columns"
      )
    )
  }
  
  h2_row <- dat %>%
    filter(source == "V(G)/Vp") %>%
    slice_head(n = 1)
  
  p_row <- dat %>%
    filter(source == "Pval") %>%
    slice_head(n = 1)
  
  n_row <- dat %>%
    filter(source == "n") %>%
    slice_head(n = 1)
  
  tibble(
    clock_folder = clock_folder,
    mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
    modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
    organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
    method = "raw_gcta_hereg",
    h2_mean = ifelse(nrow(h2_row) == 0, NA_real_, safe_numeric(h2_row$variance[1])),
    h2_se = ifelse(nrow(h2_row) == 0, NA_real_, safe_numeric(h2_row$se[1])),
    p_value = ifelse(nrow(p_row) == 0, NA_real_, safe_numeric(p_row$variance[1])),
    n = ifelse(nrow(n_row) == 0, NA_real_, safe_numeric(n_row$variance[1])),
    raw_method_source = "GCTA_REML_hsq",
    source_file = file,
    read_status = ifelse(nrow(h2_row) == 0, "missing_h2_row", "ok")
  )
}

read_hereg_one <- function(file, clock_folder, prefer_block = "HE-CP") {
  lines <- readLines(file, warn = FALSE)
  lines_trim <- str_trim(lines)
  
  block_idx <- which(lines_trim == prefer_block)
  
  if (length(block_idx) == 0) {
    block_idx <- which(lines_trim == "HE-SD")
    prefer_block <- "HE-SD"
  }
  
  if (length(block_idx) == 0) {
    return(
      tibble(
        clock_folder = clock_folder,
        mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
        modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
        organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
        method = "raw_gcta_hereg",
        h2_mean = NA_real_,
        h2_se = NA_real_,
        p_value = NA_real_,
        n = NA_real_,
        raw_method_source = "HEreg",
        source_file = file,
        read_status = "missing_HE_block"
      )
    )
  }
  
  start <- block_idx[1]
  header_i <- start + 1
  data_start <- start + 2
  
  next_block <- which(
    seq_along(lines_trim) > data_start &
      lines_trim %in% c("HE-CP", "HE-SD")
  )
  
  if (length(next_block) > 0) {
    data_end <- next_block[1] - 1
  } else {
    data_end <- length(lines_trim)
  }
  
  block_lines <- lines[data_start:data_end]
  block_lines <- block_lines[str_trim(block_lines) != ""]
  
  if (length(block_lines) == 0) {
    return(
      tibble(
        clock_folder = clock_folder,
        mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
        modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
        organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
        method = "raw_gcta_hereg",
        h2_mean = NA_real_,
        h2_se = NA_real_,
        p_value = NA_real_,
        n = NA_real_,
        raw_method_source = paste0("HEreg_", prefer_block),
        source_file = file,
        read_status = "empty_HE_block"
      )
    )
  }
  
  header_line <- lines[header_i]
  
  dat <- tryCatch(
    read.table(
      text = paste(c(header_line, block_lines), collapse = "\n"),
      header = TRUE,
      stringsAsFactors = FALSE,
      fill = TRUE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0) {
    return(
      tibble(
        clock_folder = clock_folder,
        mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
        modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
        organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
        method = "raw_gcta_hereg",
        h2_mean = NA_real_,
        h2_se = NA_real_,
        p_value = NA_real_,
        n = NA_real_,
        raw_method_source = paste0("HEreg_", prefer_block),
        source_file = file,
        read_status = "failed_parse_HE_block"
      )
    )
  }
  
  colnames(dat) <- clean_colname(colnames(dat))
  
  h2_row <- dat %>%
    filter(coefficient == "V(G)/Vp") %>%
    slice_head(n = 1)
  
  if (nrow(h2_row) == 0) {
    return(
      tibble(
        clock_folder = clock_folder,
        mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
        modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
        organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
        method = "raw_gcta_hereg",
        h2_mean = NA_real_,
        h2_se = NA_real_,
        p_value = NA_real_,
        n = NA_real_,
        raw_method_source = paste0("HEreg_", prefer_block),
        source_file = file,
        read_status = "missing_h2_row"
      )
    )
  }
  
  h2_se_value <- if ("se_jackknife" %in% colnames(h2_row)) {
    safe_numeric(h2_row$se_jackknife[1])
  } else if ("se_ols" %in% colnames(h2_row)) {
    safe_numeric(h2_row$se_ols[1])
  } else {
    NA_real_
  }
  
  p_value <- if ("p_jackknife" %in% colnames(h2_row)) {
    safe_numeric(h2_row$p_jackknife[1])
  } else if ("p_ols" %in% colnames(h2_row)) {
    safe_numeric(h2_row$p_ols[1])
  } else {
    NA_real_
  }
  
  tibble(
    clock_folder = clock_folder,
    mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
    modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
    organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
    method = "raw_gcta_hereg",
    h2_mean = safe_numeric(h2_row$estimate[1]),
    h2_se = h2_se_value,
    p_value = p_value,
    n = NA_real_,
    raw_method_source = paste0("HEreg_", prefer_block),
    source_file = file,
    read_status = "ok"
  )
}

read_raw_one_clock <- function(clock_dir) {
  clock_folder <- basename(clock_dir)
  
  hsq_files <- list.files(
    clock_dir,
    pattern = "\\.hsq$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  hereg_files <- list.files(
    clock_dir,
    pattern = "\\.hereg$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(hsq_files) > 0) {
    return(read_hsq_one(hsq_files[1], clock_folder))
  }
  
  if (length(hereg_files) > 0) {
    return(read_hereg_one(hereg_files[1], clock_folder, prefer_block = "HE-CP"))
  }
  
  tibble(
    clock_folder = clock_folder,
    mortality_clock = str_remove(clock_folder, "_mortality_clock$"),
    modality = parse_clock_modality(str_remove(clock_folder, "_mortality_clock$")),
    organ = parse_clock_organ(str_remove(clock_folder, "_mortality_clock$")),
    method = "raw_gcta_hereg",
    h2_mean = NA_real_,
    h2_se = NA_real_,
    p_value = NA_real_,
    n = NA_real_,
    raw_method_source = NA_character_,
    source_file = NA_character_,
    read_status = "missing_hsq_or_HEreg"
  )
}

clock_dirs <- list.dirs(
  gcta_hereg_dir,
  full.names = TRUE,
  recursive = FALSE
)

raw_gcta_hereg <- map_dfr(
  clock_dirs,
  read_raw_one_clock
) %>%
  mutate(
    mortality_clock = as.character(mortality_clock),
    modality = if_else(is.na(modality) | modality == "Other",
                       parse_clock_modality(mortality_clock),
                       modality),
    organ = parse_clock_organ(mortality_clock),
    h2_mean = safe_numeric(h2_mean),
    h2_se = safe_numeric(h2_se),
    p_value = safe_numeric(p_value)
  )

write_tsv(
  raw_gcta_hereg,
  file.path(output_dir, "h2_raw_gcta_hereg_parsed.tsv")
)

message("Raw GCTA/HEreg read status:")
print(raw_gcta_hereg %>% count(read_status, raw_method_source, name = "n"))

# ============================================================
# 4. Read LDSC and SBayesS summaries
# ============================================================

ldsc_tbl <- read_tsv(
  ldsc_file,
  show_col_types = FALSE,
  progress = FALSE
)

sbayess_tbl <- read_tsv(
  sbayess_file,
  show_col_types = FALSE,
  progress = FALSE
)

colnames(ldsc_tbl) <- clean_colname(colnames(ldsc_tbl))
colnames(sbayess_tbl) <- clean_colname(colnames(sbayess_tbl))

required_ldsc_cols <- c("mortality_clock", "h2_mean")
required_sbayess_cols <- c("mortality_clock", "h2_mean")

if (!all(required_ldsc_cols %in% colnames(ldsc_tbl))) {
  stop(
    "LDSC file missing required columns: ",
    paste(setdiff(required_ldsc_cols, colnames(ldsc_tbl)), collapse = ", ")
  )
}

if (!all(required_sbayess_cols %in% colnames(sbayess_tbl))) {
  stop(
    "SBayesS file missing required columns: ",
    paste(setdiff(required_sbayess_cols, colnames(sbayess_tbl)), collapse = ", ")
  )
}

ldsc_h2 <- ldsc_tbl %>%
  transmute(
    clock_folder = if ("clock_folder" %in% colnames(ldsc_tbl)) as.character(clock_folder) else NA_character_,
    mortality_clock = as.character(mortality_clock),
    modality = if ("modality" %in% colnames(ldsc_tbl)) as.character(modality) else parse_clock_modality(mortality_clock),
    organ = parse_clock_organ(mortality_clock),
    method = "ldsc",
    h2_mean = safe_numeric(h2_mean),
    h2_se = if ("h2_std" %in% colnames(ldsc_tbl)) safe_numeric(h2_std) else NA_real_,
    p_value = NA_real_,
    raw_method_source = "LDSC",
    source_file = ldsc_file,
    read_status = if ("parse_status" %in% colnames(ldsc_tbl)) as.character(parse_status) else "ok"
  ) %>%
  filter(is.finite(h2_mean))

sbayess_h2 <- sbayess_tbl %>%
  transmute(
    clock_folder = if ("clock_folder" %in% colnames(sbayess_tbl)) as.character(clock_folder) else NA_character_,
    mortality_clock = as.character(mortality_clock),
    modality = if ("modality" %in% colnames(sbayess_tbl)) as.character(modality) else parse_clock_modality(mortality_clock),
    organ = parse_clock_organ(mortality_clock),
    method = "sbayess",
    h2_mean = safe_numeric(h2_mean),
    h2_se = if ("h2_se" %in% colnames(sbayess_tbl)) safe_numeric(h2_se) else NA_real_,
    p_value = NA_real_,
    raw_method_source = "GCTB_SBayesS",
    source_file = sbayess_file,
    read_status = if ("status" %in% colnames(sbayess_tbl)) as.character(status) else "ok"
  ) %>%
  filter(is.finite(h2_mean))

# ============================================================
# 5. Combine and build df_wide
# ============================================================

df_h2_long <- bind_rows(
  raw_gcta_hereg %>%
    filter(read_status == "ok", is.finite(h2_mean)) %>%
    select(
      clock_folder,
      mortality_clock,
      modality,
      organ,
      method,
      h2_mean,
      h2_se,
      p_value,
      raw_method_source,
      source_file,
      read_status
    ),
  sbayess_h2,
  ldsc_h2
) %>%
  mutate(
    mortality_clock_chr = as.character(mortality_clock),
    modality = factor(as.character(modality), levels = modality_order),
    organ = as.character(organ),
    method = factor(method, levels = method_order),
    h2_mean = safe_numeric(h2_mean),
    h2_se = safe_numeric(h2_se)
  ) %>%
  filter(
    !is.na(method),
    is.finite(h2_mean)
  )

write_tsv(
  df_h2_long,
  file.path(output_dir, "h2_3method_summary_long.tsv")
)

message("h2 long-format counts:")
print(df_h2_long %>% count(method, raw_method_source, name = "n"))

df_wide <- df_h2_long %>%
  group_by(
    mortality_clock_chr,
    modality,
    organ,
    method
  ) %>%
  summarise(
    h2_mean = mean(h2_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = c(mortality_clock_chr, modality, organ),
    names_from = method,
    values_from = h2_mean,
    names_prefix = "h2_"
  ) %>%
  arrange(
    modality,
    organ,
    mortality_clock_chr
  )

write_tsv(
  df_wide,
  file.path(output_dir, "h2_3method_summary_wide.tsv")
)

message("df_wide dimensions: ", nrow(df_wide), " rows x ", ncol(df_wide), " columns")
print(colnames(df_wide))

# ============================================================
# 6. Pairwise consistency statistics
# ============================================================

method_pairs <- list(
  c("raw_gcta_hereg", "ldsc"),
  c("raw_gcta_hereg", "sbayess"),
  c("sbayess", "ldsc")
)

pairwise_stats <- map_dfr(method_pairs, function(pair) {
  
  m1 <- pair[1]
  m2 <- pair[2]
  
  x_col <- paste0("h2_", m1)
  y_col <- paste0("h2_", m2)
  
  x <- df_wide[[x_col]]
  y <- df_wide[[y_col]]
  
  ok <- is.finite(x) & is.finite(y)
  
  x_ok <- x[ok]
  y_ok <- y[ok]
  diff <- x_ok - y_ok
  
  pearson <- safe_cor_test(x, y, method = "pearson")
  spearman <- safe_cor_test(x, y, method = "spearman")
  
  paired_t <- if (length(diff) >= 3) {
    suppressWarnings(t.test(x_ok, y_ok, paired = TRUE))
  } else {
    NULL
  }
  
  paired_w <- if (length(diff) >= 3) {
    suppressWarnings(wilcox.test(x_ok, y_ok, paired = TRUE, exact = FALSE))
  } else {
    NULL
  }
  
  tibble(
    method_1 = m1,
    method_2 = m2,
    method_1_label = unname(method_labels[m1]),
    method_2_label = unname(method_labels[m2]),
    comparison = paste(unname(method_labels[m1]), "vs", unname(method_labels[m2])),
    n_clocks = sum(ok),
    
    pearson_r = pearson$estimate,
    pearson_p = pearson$p_value,
    
    spearman_rho = spearman$estimate,
    spearman_p = spearman$p_value,
    
    concordance_correlation = lin_ccc(x, y),
    
    mean_h2_method_1 = mean(x_ok, na.rm = TRUE),
    mean_h2_method_2 = mean(y_ok, na.rm = TRUE),
    
    mean_difference_method1_minus_method2 = mean(diff, na.rm = TRUE),
    sd_difference = sd(diff, na.rm = TRUE),
    median_difference = median(diff, na.rm = TRUE),
    
    mean_absolute_difference = mean(abs(diff), na.rm = TRUE),
    root_mean_square_difference = sqrt(mean(diff^2, na.rm = TRUE)),
    
    paired_t_p = ifelse(is.null(paired_t), NA_real_, paired_t$p.value),
    paired_wilcoxon_p = ifelse(is.null(paired_w), NA_real_, paired_w$p.value),
    
    annotation_label = paste0(
      "Pearson r = ",
      ifelse(is.finite(pearson$estimate), sprintf("%.2f", pearson$estimate), "NA"),
      "\n",
      format_pvalue(pearson$p_value)
    )
  )
}) %>%
  mutate(
    comparison = factor(
      comparison,
      levels = c(
        "GCTA/HEreg raw genotype vs LDSC summary",
        "GCTA/HEreg raw genotype vs SBayesS summary",
        "SBayesS summary vs LDSC summary"
      )
    )
  )

write_tsv(
  pairwise_stats,
  file.path(output_dir, "h2_3method_pairwise_consistency_statistics.tsv")
)

print(pairwise_stats)

# ============================================================
# 7. Panel A: h2 estimates across clocks and methods
# ============================================================

df_plot_long <- df_h2_long %>%
  mutate(
    method_label = factor(
      unname(method_labels[as.character(method)]),
      levels = unname(method_labels[method_order])
    ),
    modality = factor(as.character(modality), levels = modality_order),
    clock_label = pretty_clock_label_one_line(mortality_clock_chr)
  )

clock_order_tbl <- df_plot_long %>%
  group_by(mortality_clock_chr, modality, organ, clock_label) %>%
  summarise(
    mean_h2_across_methods = mean(h2_mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    modality,
    organ,
    mortality_clock_chr
  ) %>%
  mutate(
    x_index = row_number()
  )

df_plot_long <- df_plot_long %>%
  left_join(
    clock_order_tbl %>%
      select(mortality_clock_chr, x_index, clock_label),
    by = "mortality_clock_chr",
    suffix = c("", "_ordered")
  ) %>%
  mutate(
    clock_label_plot = clock_label_ordered
  )

modality_bounds <- clock_order_tbl %>%
  group_by(modality) %>%
  summarise(
    x_min = min(x_index),
    x_max = max(x_index),
    x_mid = mean(range(x_index)),
    .groups = "drop"
  )

separator_tbl <- modality_bounds %>%
  arrange(x_min) %>%
  mutate(separator_x = lag(x_min) - 0.5) %>%
  filter(is.finite(separator_x))

p_line <- ggplot(
  df_plot_long,
  aes(
    x = x_index,
    y = h2_mean,
    color = method_label,
    group = method_label
  )
) +
  geom_vline(
    data = separator_tbl,
    aes(xintercept = separator_x),
    inherit.aes = FALSE,
    linetype = "dotted",
    color = "grey50",
    linewidth = 0.45
  ) +
  geom_line(
    linewidth = 0.75,
    alpha = 0.95
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(h2_mean - h2_se, 0),
      ymax = h2_mean + h2_se
    ),
    width = 0.12,
    linewidth = 0.35,
    alpha = 0.75,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.0,
    alpha = 0.98
  ) +
  geom_text(
    data = modality_bounds,
    aes(
      x = x_mid,
      y = max(df_plot_long$h2_mean + df_plot_long$h2_se, na.rm = TRUE) * 1.05,
      label = modality
    ),
    inherit.aes = FALSE,
    size = 3.1,
    fontface = "bold",
    color = "#1F2A44"
  ) +
  scale_color_manual(
    values = method_colors,
    name = "Method"
  ) +
  scale_x_continuous(
    breaks = clock_order_tbl$x_index,
    labels = clock_order_tbl$clock_label,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(df_plot_long$h2_mean + df_plot_long$h2_se, na.rm = TRUE) * 1.16
    ),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.00, 0.02))
  ) +
  labs(
    x = NULL,
    y = expression("SNP heritability estimate " %+-% " SE")
  ) +
  theme_pub(base_size = 10) +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.box.margin = margin(-6, 0, -4, 0),
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 7.2
    ),
    axis.text.y = element_text(size = 8.5),
    axis.title.y = element_text(size = 10.2),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(output_dir, "h2_3method_panel_A_line_plot.pdf"),
  p_line,
  width = 10.5,
  height = 3.2,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_panel_A_line_plot.png"),
  p_line,
  width = 10.5,
  height = 3.2,
  dpi = 500
)

# ============================================================
# 8. Panel B: pairwise agreement scatter plots
# ============================================================

df_scatter <- map_dfr(method_pairs, function(pair) {
  
  m1 <- pair[1]
  m2 <- pair[2]
  
  x_col <- paste0("h2_", m1)
  y_col <- paste0("h2_", m2)
  
  tibble(
    mortality_clock = df_wide$mortality_clock_chr,
    modality = df_wide$modality,
    organ = df_wide$organ,
    method_1 = m1,
    method_2 = m2,
    method_1_label = unname(method_labels[m1]),
    method_2_label = unname(method_labels[m2]),
    comparison = paste(unname(method_labels[m1]), "vs", unname(method_labels[m2])),
    h2_method_1 = df_wide[[x_col]],
    h2_method_2 = df_wide[[y_col]]
  )
}) %>%
  drop_na(h2_method_1, h2_method_2) %>%
  mutate(
    modality = factor(as.character(modality), levels = modality_order),
    comparison = factor(
      comparison,
      levels = levels(pairwise_stats$comparison)
    )
  )

scatter_annotation <- df_scatter %>%
  group_by(comparison) %>%
  summarise(
    x_min = min(h2_method_1, na.rm = TRUE),
    x_max = max(h2_method_1, na.rm = TRUE),
    y_min = min(h2_method_2, na.rm = TRUE),
    y_max = max(h2_method_2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_range = x_max - x_min,
    y_range = y_max - y_min,
    x_pos = x_min + 0.06 * ifelse(x_range == 0, 1, x_range),
    y_pos = y_max - 0.08 * ifelse(y_range == 0, 1, y_range)
  ) %>%
  left_join(
    pairwise_stats %>%
      select(comparison, annotation_label),
    by = "comparison"
  )

p_scatter <- ggplot(
  df_scatter,
  aes(
    x = h2_method_1,
    y = h2_method_2,
    shape = modality
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "#555555",
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.7,
    alpha = 0.95,
    color = "#376B9B"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.75,
    color = "#E5A11A"
  ) +
  geom_label(
    data = scatter_annotation,
    aes(
      x = x_pos,
      y = y_pos,
      label = annotation_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    fontface = "bold",
    color = "#1F2A44",
    fill = "#F7FAFF",
    label.size = 0.25,
    label.padding = unit(0.16, "lines")
  ) +
  facet_wrap(
    ~ comparison,
    scales = "free",
    nrow = 1
  ) +
  scale_shape_manual(
    values = modality_shapes,
    name = "Modality",
    drop = FALSE
  ) +
  labs(
    title = expression("Pairwise agreement of " * h^2 * " estimates across methods"),
    x = expression("Method 1 " * h^2),
    y = expression("Method 2 " * h^2)
  ) +
  theme_pub(base_size = 10) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 13.2),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 7.2),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(output_dir, "h2_3method_panel_B_pairwise_scatter.pdf"),
  p_scatter,
  width = 10.5,
  height = 3.3,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_panel_B_pairwise_scatter.png"),
  p_scatter,
  width = 10.5,
  height = 3.3,
  dpi = 500
)

# ============================================================
# 9. Panel C: Bland-Altman plots
# ============================================================

df_bland <- map_dfr(method_pairs, function(pair) {
  
  m1 <- pair[1]
  m2 <- pair[2]
  
  x_col <- paste0("h2_", m1)
  y_col <- paste0("h2_", m2)
  
  tibble(
    mortality_clock = df_wide$mortality_clock_chr,
    modality = df_wide$modality,
    organ = df_wide$organ,
    method_1 = m1,
    method_2 = m2,
    method_1_label = unname(method_labels[m1]),
    method_2_label = unname(method_labels[m2]),
    comparison = paste(unname(method_labels[m1]), "vs", unname(method_labels[m2])),
    h2_method_1 = df_wide[[x_col]],
    h2_method_2 = df_wide[[y_col]]
  )
}) %>%
  filter(
    is.finite(h2_method_1),
    is.finite(h2_method_2)
  ) %>%
  mutate(
    modality = factor(as.character(modality), levels = modality_order),
    comparison = factor(
      comparison,
      levels = levels(pairwise_stats$comparison)
    ),
    mean_h2 = (h2_method_1 + h2_method_2) / 2,
    diff_h2 = h2_method_1 - h2_method_2
  )

bland_summary <- df_bland %>%
  group_by(comparison) %>%
  summarise(
    mean_diff = mean(diff_h2, na.rm = TRUE),
    sd_diff = sd(diff_h2, na.rm = TRUE),
    loa_low = mean_diff - 1.96 * sd_diff,
    loa_high = mean_diff + 1.96 * sd_diff,
    .groups = "drop"
  )

write_tsv(
  bland_summary,
  file.path(output_dir, "h2_3method_bland_altman_summary.tsv")
)

p_bland <- ggplot(
  df_bland,
  aes(
    x = mean_h2,
    y = diff_h2,
    shape = modality
  )
) +
  geom_hline(
    data = bland_summary,
    aes(yintercept = mean_diff),
    color = "#376B9B",
    linewidth = 0.55
  ) +
  geom_hline(
    data = bland_summary,
    aes(yintercept = loa_low),
    color = "#C53A3A",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_hline(
    data = bland_summary,
    aes(yintercept = loa_high),
    color = "#C53A3A",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = 0,
    color = "grey35",
    linetype = "dotted",
    linewidth = 0.45
  ) +
  geom_point(
    size = 2.7,
    alpha = 0.95,
    color = "#5C9465"
  ) +
  facet_wrap(
    ~ comparison,
    scales = "free",
    nrow = 1
  ) +
  scale_shape_manual(
    values = modality_shapes,
    name = NULL,
    drop = FALSE
  ) +
  labs(
    title = expression("Bland\u2013Altman comparison of " * h^2 * " estimates"),
    x = expression("Mean " * h^2 * " between two methods"),
    y = expression("Difference in " * h^2 * ", method 1 \u2212 method 2")
  ) +
  theme_pub(base_size = 10) +
  theme(
    legend.position = "top",
    plot.title = element_text(size = 13.2),
    axis.title = element_text(size = 9.5),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 7.2),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(
  file.path(output_dir, "h2_3method_panel_C_bland_altman.pdf"),
  p_bland,
  width = 10.5,
  height = 3.3,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_panel_C_bland_altman.png"),
  p_bland,
  width = 10.5,
  height = 3.3,
  dpi = 500
)

# ============================================================
# 10. Combined 3-panel figure
# ============================================================

combined_plot <- p_line / p_scatter / p_bland +
  plot_layout(
    heights = c(1.05, 1.00, 1.00)
  ) +
  plot_annotation(
    title = expression("Consistency of " * h^2 * " estimates across three genetic-architecture methods"),
    subtitle = paste0(
      "Raw-genotype GCTA/HEreg, summary-level SBayesS, and summary-level LDSC across ",
      nrow(df_wide),
      " mortality L\u2019EPOCH clocks"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 17,
        color = "#1F2A44",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = 10.5,
        color = "#4A4A4A",
        hjust = 0
      ),
      plot.tag = element_text(
        face = "bold",
        size = 24,
        color = "black"
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )
  )

ggsave(
  file.path(output_dir, "h2_3method_consistency_3panel_figure.pdf"),
  combined_plot,
  width = 10.8,
  height = 9.8,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_consistency_3panel_figure.png"),
  combined_plot,
  width = 10.8,
  height = 9.8,
  dpi = 500
)

# Keep old output names for downstream compatibility.
ggsave(
  file.path(output_dir, "h2_3method_consistency_combined_figure.pdf"),
  combined_plot,
  width = 10.8,
  height = 9.8,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_consistency_combined_figure.png"),
  combined_plot,
  width = 10.8,
  height = 9.8,
  dpi = 500
)

# ============================================================
# 11. Completion message
# ============================================================

message("============================================================")
message("Done.")
message("Saved outputs to:")
message("  ", output_dir)
message("")
message("Main files:")
message("  ", file.path(output_dir, "h2_3method_summary_long.tsv"))
message("  ", file.path(output_dir, "h2_3method_summary_wide.tsv"))
message("  ", file.path(output_dir, "h2_3method_pairwise_consistency_statistics.tsv"))
message("  ", file.path(output_dir, "h2_3method_panel_A_line_plot.pdf"))
message("  ", file.path(output_dir, "h2_3method_panel_B_pairwise_scatter.pdf"))
message("  ", file.path(output_dir, "h2_3method_panel_C_bland_altman.pdf"))
message("  ", file.path(output_dir, "h2_3method_consistency_3panel_figure.pdf"))
message("  ", file.path(output_dir, "h2_3method_consistency_3panel_figure.png"))
message("============================================================")