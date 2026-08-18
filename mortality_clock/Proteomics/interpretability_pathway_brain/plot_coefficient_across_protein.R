#!/usr/bin/env Rscript

# ==============================================================================
# Brain proteomics mortality EPOCH:
# plot Elastic-Net coefficients + create STRING-ready input files
#
# REVISION:
#   ZERO-coefficient proteins are EXCLUDED from ALL STRING input files.
#
# INPUTS
# ------
# 1) brain_proteomics_mortality_clock_nonzero_coefficients.tsv
#    -> used for plotting the selected Elastic-Net protein coefficients
#
# 2) brain_proteomics_mortality_clock_coefficients.tsv
#    -> used to identify all brain-proteomics coefficients
#    -> coefficient == 0 proteins are removed before creating STRING files
#
# OUTPUT DIRECTORY
# ----------------
# /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/
# Brain_proteomics_mortality_clock/STRING_brain_EPOCH_coefficients
#
# STRING OUTPUTS
# --------------
# 1) STRING_values_ranks_absolute_coefficients.tsv
#    protein<TAB>|coefficient|
#    NO HEADER
#    ONLY non-zero brain-proteomics proteins
#
# 2) STRING_values_ranks_signed_coefficients.tsv
#    protein<TAB>signed coefficient
#    NO HEADER
#    ONLY non-zero brain-proteomics proteins
#
# 3) STRING_top20_proteins_multiple_input.txt
#    top 20 proteins by |coefficient|, one protein per line
#
# 4) STRING_all_candidate_brain_proteins_background.txt
#    retained for compatibility with the earlier workflow, but now contains
#    ONLY non-zero selected brain-proteomics proteins
#
# 5) STRING_all_selected_nonzero_proteins_multiple_input.txt
#    all non-zero selected proteins, one protein per line
#
# 6) STRING_excluded_zero_coefficient_proteins.tsv
#    audit file listing proteins removed from the STRING inputs
#
# STRING organism:
#    Homo sapiens (9606)
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ==============================================================================
# 1. Paths and settings
# ==============================================================================

base_dir <- paste0(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/",
  "Brain_proteomics_mortality_clock"
)

nonzero_file <- file.path(
  base_dir,
  "brain_proteomics_mortality_clock_nonzero_coefficients.tsv"
)

full_coeff_file <- file.path(
  base_dir,
  "brain_proteomics_mortality_clock_coefficients.tsv"
)

output_dir <- file.path(
  base_dir,
  "STRING_brain_EPOCH_coefficients"
)

# 0 = plot all selected proteins.
top_n <- 0L

# Number of proteins used for the compact STRING PPI network.
string_network_top_n <- 20L

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Output directory:")
message("  ", output_dir)

# ==============================================================================
# 2. Input checks
# ==============================================================================

if (!file.exists(nonzero_file)) {
  stop(
    "Non-zero coefficient file not found:\n",
    nonzero_file
  )
}

if (!file.exists(full_coeff_file)) {
  stop(
    "Full coefficient file not found:\n",
    full_coeff_file
  )
}

# ==============================================================================
# 3. Helper functions
# ==============================================================================

to_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  tolower(as.character(x)) %in% c(
    "true",
    "t",
    "1",
    "yes"
  )
}

clean_protein_name <- function(x) {
  x %>%
    str_remove("^num__") %>%
    str_trim()
}

read_coefficient_file <- function(path) {
  
  x <- read_tsv(
    path,
    show_col_types = FALSE,
    na = c("", "NA", "NaN")
  )
  
  required <- c(
    "feature",
    "coefficient"
  )
  
  missing_cols <- setdiff(
    required,
    names(x)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s) in ",
      basename(path),
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  x <- x %>%
    mutate(
      coefficient = as.numeric(coefficient),
      abs_coefficient = abs(coefficient)
    )
  
  if ("is_nonzero" %in% names(x)) {
    x <- x %>%
      mutate(
        is_nonzero = to_logical(is_nonzero)
      )
  } else {
    x <- x %>%
      mutate(
        is_nonzero =
          is.finite(coefficient) &
          coefficient != 0
      )
  }
  
  if ("is_brain_proteomics_feature" %in% names(x)) {
    
    x <- x %>%
      mutate(
        is_brain_proteomics_feature =
          to_logical(is_brain_proteomics_feature)
      )
    
  } else if ("penalty_factor" %in% names(x)) {
    
    # In this model:
    #   penalty_factor = 1 -> penalized proteomics feature
    #   penalty_factor = 0 -> unpenalized covariate
    x <- x %>%
      mutate(
        is_brain_proteomics_feature =
          suppressWarnings(
            as.numeric(penalty_factor)
          ) == 1
      )
    
  } else {
    
    stop(
      paste0(
        "Cannot identify brain-proteomics rows in ",
        basename(path),
        ". Expected either 'is_brain_proteomics_feature' ",
        "or 'penalty_factor'."
      )
    )
  }
  
  x
}

write_no_header_tsv <- function(df, path) {
  
  write.table(
    df,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

write_one_per_line <- function(x, path) {
  
  writeLines(
    as.character(x),
    con = path,
    sep = "\n",
    useBytes = TRUE
  )
}

validate_values_ranks <- function(path, expected_n) {
  
  x <- read.delim(
    path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (ncol(x) != 2) {
    stop(
      "STRING Values/Ranks file does not have exactly 2 columns:\n",
      path
    )
  }
  
  if (nrow(x) != expected_n) {
    stop(
      "Unexpected number of STRING rows in ",
      path,
      "\nObserved: ",
      nrow(x),
      "\nExpected: ",
      expected_n
    )
  }
  
  if (
    any(
      is.na(x[[1]]) |
      trimws(x[[1]]) == ""
    )
  ) {
    stop(
      "Missing protein identifier(s) in STRING file:\n",
      path
    )
  }
  
  values <- suppressWarnings(
    as.numeric(x[[2]])
  )
  
  if (any(!is.finite(values))) {
    stop(
      "Non-numeric value(s) detected in STRING file:\n",
      path
    )
  }
  
  # Explicitly enforce the requested zero-coefficient removal.
  if (any(values == 0)) {
    stop(
      "QC failure: zero values remain in STRING file:\n",
      path
    )
  }
  
  invisible(TRUE)
}

# ==============================================================================
# 4. Read selected non-zero coefficients for plotting
# ==============================================================================

dat_nonzero <- read_coefficient_file(
  nonzero_file
)

prot_selected_all <- dat_nonzero %>%
  filter(
    is_brain_proteomics_feature,
    is.finite(coefficient),
    coefficient != 0
  ) %>%
  mutate(
    protein = clean_protein_name(feature),
    abs_coefficient = abs(coefficient),
    direction = if_else(
      coefficient > 0,
      "Positive coefficient",
      "Negative coefficient"
    )
  ) %>%
  arrange(
    desc(abs_coefficient),
    protein
  )

if (nrow(prot_selected_all) == 0) {
  stop(
    "No non-zero brain-proteomics features found in ",
    basename(nonzero_file)
  )
}

if (anyDuplicated(prot_selected_all$protein)) {
  
  duplicated_proteins <- unique(
    prot_selected_all$protein[
      duplicated(prot_selected_all$protein)
    ]
  )
  
  stop(
    "Duplicated selected protein names detected: ",
    paste(
      duplicated_proteins,
      collapse = ", "
    )
  )
}

n_selected_total <- nrow(
  prot_selected_all
)

# ==============================================================================
# 5. Read full coefficients and REMOVE zero-coefficient proteins
# ==============================================================================

dat_full <- read_coefficient_file(
  full_coeff_file
)

# Audit list: these proteins are excluded from STRING.
zero_proteins <- dat_full %>%
  filter(
    is_brain_proteomics_feature,
    is.finite(coefficient),
    coefficient == 0
  ) %>%
  mutate(
    protein = clean_protein_name(feature)
  ) %>%
  arrange(protein)

zero_audit_out <- file.path(
  output_dir,
  "STRING_excluded_zero_coefficient_proteins.tsv"
)

zero_proteins %>%
  select(
    protein,
    feature,
    coefficient,
    abs_coefficient,
    everything()
  ) %>%
  write_tsv(
    zero_audit_out
  )

# IMPORTANT:
# This is the protein set used for ALL STRING inputs.
# Zero coefficients are explicitly excluded.
prot_string <- dat_full %>%
  filter(
    is_brain_proteomics_feature,
    is.finite(coefficient),
    coefficient != 0
  ) %>%
  mutate(
    protein = clean_protein_name(feature),
    abs_coefficient = abs(coefficient),
    direction = if_else(
      coefficient > 0,
      "Positive coefficient",
      "Negative coefficient"
    )
  ) %>%
  arrange(
    desc(abs_coefficient),
    protein
  )

if (nrow(prot_string) == 0) {
  stop(
    "No non-zero brain-proteomics features found in ",
    basename(full_coeff_file)
  )
}

if (any(prot_string$coefficient == 0)) {
  stop(
    "Internal QC failure: zero-coefficient proteins remain in prot_string."
  )
}

if (anyDuplicated(prot_string$protein)) {
  
  duplicated_proteins <- unique(
    prot_string$protein[
      duplicated(prot_string$protein)
    ]
  )
  
  stop(
    "Duplicated protein names detected in STRING protein set: ",
    paste(
      duplicated_proteins,
      collapse = ", "
    )
  )
}

n_string_proteins <- nrow(
  prot_string
)

n_zero_excluded <- nrow(
  zero_proteins
)

# ==============================================================================
# 6. Prepare plotting data
# ==============================================================================

prot_plot <- prot_selected_all

if (top_n > 0) {
  
  prot_plot <- prot_plot %>%
    slice_head(
      n = min(
        top_n,
        nrow(prot_plot)
      )
    )
}

prot_plot <- prot_plot %>%
  mutate(
    rank_abs = row_number(),
    protein_plot = factor(
      protein,
      levels = rev(protein)
    )
  )

# ==============================================================================
# 7. Save ranked tables with headers
# ==============================================================================

selected_ranked_out <- file.path(
  output_dir,
  paste0(
    "brain_proteomics_EPOCH_selected_proteins_",
    "ranked_by_absolute_coefficient.tsv"
  )
)

prot_selected_all %>%
  mutate(
    rank_abs = row_number()
  ) %>%
  select(
    rank_abs,
    protein,
    coefficient,
    abs_coefficient,
    direction,
    everything()
  ) %>%
  write_tsv(
    selected_ranked_out
  )

string_ranked_audit_out <- file.path(
  output_dir,
  paste0(
    "brain_proteomics_EPOCH_STRING_nonzero_proteins_",
    "ranked_by_absolute_coefficient.tsv"
  )
)

prot_string %>%
  mutate(
    rank_abs = row_number()
  ) %>%
  select(
    rank_abs,
    protein,
    coefficient,
    abs_coefficient,
    direction,
    everything()
  ) %>%
  write_tsv(
    string_ranked_audit_out
  )

# ==============================================================================
# 8. Create STRING-ready files: NON-ZERO proteins only
# ==============================================================================

# ------------------------------------------------------------------------------
# 8A. PRIMARY ranked analysis:
#     protein<TAB>|Elastic-Net coefficient|
#     no header
# ------------------------------------------------------------------------------

string_abs_out <- file.path(
  output_dir,
  "STRING_values_ranks_absolute_coefficients.tsv"
)

string_abs <- prot_string %>%
  select(
    protein,
    value = abs_coefficient
  )

write_no_header_tsv(
  string_abs,
  string_abs_out
)

# ------------------------------------------------------------------------------
# 8B. SECONDARY direction-aware ranked analysis:
#     protein<TAB>signed Elastic-Net coefficient
#     no header
# ------------------------------------------------------------------------------

string_signed_out <- file.path(
  output_dir,
  "STRING_values_ranks_signed_coefficients.tsv"
)

string_signed <- prot_string %>%
  select(
    protein,
    value = coefficient
  )

write_no_header_tsv(
  string_signed,
  string_signed_out
)

# ------------------------------------------------------------------------------
# 8C. Top-20 proteins for a compact STRING PPI network
# ------------------------------------------------------------------------------

top_network_n <- min(
  string_network_top_n,
  nrow(prot_string)
)

top20 <- prot_string %>%
  slice_head(
    n = top_network_n
  )

string_top20_out <- file.path(
  output_dir,
  "STRING_top20_proteins_multiple_input.txt"
)

write_one_per_line(
  top20$protein,
  string_top20_out
)

string_top20_audit_out <- file.path(
  output_dir,
  "STRING_top20_proteins_with_coefficients.tsv"
)

top20 %>%
  mutate(
    rank_abs = row_number()
  ) %>%
  select(
    rank_abs,
    protein,
    coefficient,
    abs_coefficient,
    direction
  ) %>%
  write_tsv(
    string_top20_audit_out
  )

# ------------------------------------------------------------------------------
# 8D. Background/input list:
#     for compatibility, keep the original filename, but it now contains
#     ONLY non-zero selected proteins
# ------------------------------------------------------------------------------

string_background_out <- file.path(
  output_dir,
  "STRING_all_candidate_brain_proteins_background.txt"
)

write_one_per_line(
  prot_string$protein,
  string_background_out
)

# More explicit duplicate filename for clarity.
string_nonzero_background_out <- file.path(
  output_dir,
  "STRING_all_nonzero_selected_brain_proteins_background.txt"
)

write_one_per_line(
  prot_string$protein,
  string_nonzero_background_out
)

# ------------------------------------------------------------------------------
# 8E. All selected non-zero proteins for STRING "Multiple proteins"
# ------------------------------------------------------------------------------

string_selected_out <- file.path(
  output_dir,
  "STRING_all_selected_nonzero_proteins_multiple_input.txt"
)

write_one_per_line(
  prot_string$protein,
  string_selected_out
)

# ==============================================================================
# 9. Validate STRING files
# ==============================================================================

validate_values_ranks(
  string_abs_out,
  n_string_proteins
)

validate_values_ranks(
  string_signed_out,
  n_string_proteins
)

if (any(string_abs$value == 0)) {
  stop(
    "QC failure: zero values remain in absolute STRING input."
  )
}

if (any(string_signed$value == 0)) {
  stop(
    "QC failure: zero values remain in signed STRING input."
  )
}

# ==============================================================================
# 10. Plot signed coefficients ranked by absolute magnitude
# ==============================================================================

subtitle_text <- if (
  top_n > 0 &&
  top_n < n_selected_total
) {
  
  paste0(
    "Top ",
    nrow(prot_plot),
    " of ",
    n_selected_total,
    " selected proteins, ranked by absolute coefficient"
  )
  
} else {
  
  paste0(
    "All ",
    n_selected_total,
    " non-zero Elastic-Net-selected brain proteins"
  )
}

p_signed <- ggplot(
  prot_plot,
  aes(
    x = coefficient,
    y = protein_plot,
    fill = direction
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  labs(
    title =
      "Brain proteomics mortality EPOCH coefficients",
    subtitle =
      subtitle_text,
    x =
      "Elastic-Net coefficient",
    y =
      NULL,
    fill =
      "Direction"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0,
      margin = margin(
        b = 8
      )
    ),
    axis.text.y = element_text(
      size = if (
        nrow(prot_plot) > 35
      ) {
        7.5
      } else {
        9
      }
    ),
    axis.text.x = element_text(
      size = 9
    ),
    axis.title.x = element_text(
      size = 10
    ),
    legend.position = "bottom",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

# ==============================================================================
# 11. Plot absolute coefficients
# ==============================================================================

p_abs <- ggplot(
  prot_plot,
  aes(
    x = abs_coefficient,
    y = protein_plot
  )
) +
  geom_col(
    width = 0.72
  ) +
  labs(
    title =
      "Brain proteomics mortality EPOCH",
    subtitle =
      subtitle_text,
    x =
      "|Elastic-Net coefficient|",
    y =
      NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0,
      margin = margin(
        b = 8
      )
    ),
    axis.text.y = element_text(
      size = if (
        nrow(prot_plot) > 35
      ) {
        7.5
      } else {
        9
      }
    ),
    axis.text.x = element_text(
      size = 9
    ),
    axis.title.x = element_text(
      size = 10
    ),
    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

# ==============================================================================
# 12. Save figures
# ==============================================================================

fig_height <- max(
  7,
  min(
    15,
    2.2 + 0.24 * nrow(prot_plot)
  )
)

fig_width <- 7.5

signed_pdf <- file.path(
  output_dir,
  paste0(
    "brain_proteomics_EPOCH_signed_coefficients_",
    "ranked_by_absolute_value.pdf"
  )
)

signed_png <- file.path(
  output_dir,
  paste0(
    "brain_proteomics_EPOCH_signed_coefficients_",
    "ranked_by_absolute_value.png"
  )
)

absolute_pdf <- file.path(
  output_dir,
  "brain_proteomics_EPOCH_absolute_coefficients_ranked.pdf"
)

absolute_png <- file.path(
  output_dir,
  "brain_proteomics_EPOCH_absolute_coefficients_ranked.png"
)

ggsave(
  filename = signed_pdf,
  plot = p_signed,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = signed_png,
  plot = p_signed,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = absolute_pdf,
  plot = p_abs,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = absolute_png,
  plot = p_abs,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 400,
  bg = "white"
)

if (
  requireNamespace(
    "svglite",
    quietly = TRUE
  )
) {
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "brain_proteomics_EPOCH_signed_coefficients_",
        "ranked_by_absolute_value.svg"
      )
    ),
    plot = p_signed,
    width = fig_width,
    height = fig_height,
    units = "in",
    device = svglite::svglite
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      "brain_proteomics_EPOCH_absolute_coefficients_ranked.svg"
    ),
    plot = p_abs,
    width = fig_width,
    height = fig_height,
    units = "in",
    device = svglite::svglite
  )
}

# ==============================================================================
# 13. README
# ==============================================================================

readme_out <- file.path(
  output_dir,
  "README_STRING_inputs.txt"
)

readme_lines <- c(
  "Brain proteomics mortality EPOCH - STRING input files",
  "=====================================================",
  "",
  paste0(
    "Non-zero proteins included in STRING: N = ",
    n_string_proteins
  ),
  paste0(
    "Zero-coefficient proteins excluded from STRING: N = ",
    n_zero_excluded
  ),
  "",
  "IMPORTANT",
  "---------",
  "All zero-coefficient proteins have been removed from every STRING input.",
  "",
  "PRIMARY RANKED ENRICHMENT",
  "-------------------------",
  "File: STRING_values_ranks_absolute_coefficients.tsv",
  "STRING input: Proteins with Values/Ranks",
  "Organism: Homo sapiens (9606)",
  "Column 1: protein/gene symbol",
  "Column 2: absolute Elastic-Net coefficient",
  "Header: NONE",
  "Contains only non-zero Elastic-Net-selected brain proteins.",
  "",
  "SECONDARY DIRECTION-AWARE RANKED ENRICHMENT",
  "-------------------------------------------",
  "File: STRING_values_ranks_signed_coefficients.tsv",
  "STRING input: Proteins with Values/Ranks",
  "Organism: Homo sapiens (9606)",
  "Column 1: protein/gene symbol",
  "Column 2: signed Elastic-Net coefficient",
  "Header: NONE",
  "Contains only non-zero Elastic-Net-selected brain proteins.",
  "",
  "TOP-20 PPI NETWORK",
  "------------------",
  "File: STRING_top20_proteins_multiple_input.txt",
  "STRING input: Multiple proteins",
  "Organism: Homo sapiens (9606)",
  "Contains the 20 non-zero proteins with the largest absolute coefficients.",
  "",
  "ALL NON-ZERO SELECTED PROTEINS",
  "------------------------------",
  "File: STRING_all_selected_nonzero_proteins_multiple_input.txt",
  "Contains all non-zero Elastic-Net-selected brain proteins.",
  "",
  "ZERO-PROTEIN AUDIT",
  "------------------",
  "File: STRING_excluded_zero_coefficient_proteins.tsv",
  "Lists proteins removed because their Elastic-Net coefficient was exactly zero."
)

writeLines(
  readme_lines,
  readme_out
)

# ==============================================================================
# 14. Console summary
# ==============================================================================

message("")
message("============================================================")
message("Brain proteomics mortality EPOCH analysis complete")
message("============================================================")

message(
  "Output directory: ",
  output_dir
)

message(
  "Non-zero proteins included in STRING: ",
  n_string_proteins
)

message(
  "Zero-coefficient proteins excluded from STRING: ",
  n_zero_excluded
)

message("")
message("STRING files:")

message(
  "  PRIMARY absolute Values/Ranks: ",
  string_abs_out
)

message(
  "  Signed Values/Ranks:           ",
  string_signed_out
)

message(
  "  Top-20 PPI input:              ",
  string_top20_out
)

message(
  "  All non-zero selected:         ",
  string_selected_out
)

message(
  "  Non-zero background:           ",
  string_nonzero_background_out
)

message(
  "  Excluded zero proteins:        ",
  zero_audit_out
)

message("")
message(
  "Top ",
  top_network_n,
  " proteins by |coefficient|:"
)

print(
  top20 %>%
    mutate(
      rank_abs = row_number()
    ) %>%
    select(
      rank_abs,
      protein,
      coefficient,
      abs_coefficient
    ),
  n = nrow(top20)
)

if (n_zero_excluded > 0) {
  message("")
  message("Excluded zero-coefficient proteins:")
  print(
    zero_proteins %>%
      select(
        protein,
        coefficient
      ),
    n = n_zero_excluded
  )
}

# Display the signed coefficient plot when sourced interactively.
if (interactive()) {
  print(p_signed)
}