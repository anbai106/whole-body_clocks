#!/usr/bin/env Rscript

# =============================================================================
# FUMA GENE2FUNC enrichment extraction and domain-ordered visualization
# Systemic EPOCH-BAG resilience-associated expression genes
#
# INPUT:
#   GS.txt      FUMA "Enrichment of input genes in Gene Sets"
#   summary.txt optional; used to recover the number of input genes
#
# OUTPUTS:
#   1) enrichment_significant_by_biological_domain.tsv
#   2) enrichment_domain_summary.tsv
#   3) enrichment_gene_contributions.tsv
#   4) enrichment_functional_pathways_only.tsv
#   5) enrichment_by_biological_domain.{png,pdf,svg}
#   6) enrichment_functional_pathways_only.{png,pdf,svg}
#
# Biological-domain grouping is intentionally transparent and editable below.
#
# Important interpretation:
#   - GO/KEGG terms are treated as functional pathways/processes.
#   - GWAS Catalog terms are prior trait-association annotations, not
#     mechanistic pathways.
#   - Positional gene sets indicate genomic clustering, not pathway biology.
#   - Terms dominated only by APOE/APOC1 are explicitly flagged.
#
# Usage on your local Mac:
#   Rscript plot_FUMA_GENE2FUNC_resilience_enrichment_local_mac.R
#
# The default local input directory is:
#   /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/
#   fuma/Brain_proteomics_mortality_clock/EPOCH_BAG_residual/gene2func
#
# The corresponding cluster directory is:
#   /gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/
#   mortality_clock/fuma/Brain_proteomics_mortality_clock/
#   EPOCH_BAG_residual/gene2func
#
# Optional custom input directory:
#   Rscript plot_FUMA_GENE2FUNC_resilience_enrichment_local_mac.R \
#       /path/to/gene2func
#
# Optional custom input and output directories:
#   Rscript plot_FUMA_GENE2FUNC_resilience_enrichment_local_mac.R \
#       /path/to/gene2func /path/to/output_directory
#
# Environment alternatives:
#   GENE2FUNC_DIR=/path/to/gene2func
#   OUT_DIR=/path/to/output
#
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(svglite)
})

# -----------------------------------------------------------------------------
# 1. INPUT / OUTPUT
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------------------------------------------------------
# Default FUMA GENE2FUNC directories
# -----------------------------------------------------------------------------
#
# Local Mac path:
#   /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/
#   fuma/Brain_proteomics_mortality_clock/EPOCH_BAG_residual/gene2func
#
# CUBIC/CBICA cluster path:
#   /gpfs/fs001/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/
#   mortality_clock/fuma/Brain_proteomics_mortality_clock/
#   EPOCH_BAG_residual/gene2func
#
# The script is intended to run directly on your Mac. It first checks the
# Mac-mounted cubic-home path and then falls back to the cluster path.
#
# You can always override the directory by:
#   Rscript plot_FUMA_GENE2FUNC_resilience_enrichment_local_mac.R \
#       /custom/path/to/gene2func
#
# or:
#   GENE2FUNC_DIR=/custom/path/to/gene2func \
#       Rscript plot_FUMA_GENE2FUNC_resilience_enrichment_local_mac.R
# -----------------------------------------------------------------------------

DEFAULT_MAC_GENE2FUNC_DIR <- file.path(
  "/Users/hao/cubic-home",
  "Reproducibile_paper",
  "WholeBodyClock",
  "mortality_clock",
  "fuma",
  "Brain_proteomics_mortality_clock",
  "EPOCH_BAG_residual",
  "gene2func"
)

DEFAULT_CLUSTER_GENE2FUNC_DIR <- file.path(
  "/gpfs/fs001/cbica/home/wenju",
  "Reproducibile_paper",
  "WholeBodyClock",
  "mortality_clock",
  "fuma",
  "Brain_proteomics_mortality_clock",
  "EPOCH_BAG_residual",
  "gene2func"
)

resolve_gene2func_dir <- function() {

  # Highest priority: first command-line argument.
  if (length(args) >= 1 && nzchar(args[[1]])) {
    return(
      normalizePath(
        path.expand(args[[1]]),
        mustWork = FALSE
      )
    )
  }

  # Second priority: environment variable.
  env_dir <- Sys.getenv(
    "GENE2FUNC_DIR",
    unset = ""
  )

  if (nzchar(env_dir)) {
    return(
      normalizePath(
        path.expand(env_dir),
        mustWork = FALSE
      )
    )
  }

  # Preferred local-Mac location.
  if (dir.exists(DEFAULT_MAC_GENE2FUNC_DIR)) {
    return(
      normalizePath(
        DEFAULT_MAC_GENE2FUNC_DIR,
        mustWork = TRUE
      )
    )
  }

  # Cluster fallback.
  if (dir.exists(DEFAULT_CLUSTER_GENE2FUNC_DIR)) {
    return(
      normalizePath(
        DEFAULT_CLUSTER_GENE2FUNC_DIR,
        mustWork = TRUE
      )
    )
  }

  stop(
    paste0(
      "Could not find the FUMA GENE2FUNC directory.\n\n",
      "Checked local Mac path:\n  ",
      DEFAULT_MAC_GENE2FUNC_DIR,
      "\n\nChecked cluster path:\n  ",
      DEFAULT_CLUSTER_GENE2FUNC_DIR,
      "\n\nEither mount cubic-home on your Mac or run with:\n",
      "  Rscript plot_FUMA_GENE2FUNC_resilience_enrichment_local_mac.R ",
      "/path/to/gene2func\n"
    )
  )
}

GENE2FUNC_DIR <- resolve_gene2func_dir()

GS_PATH <- file.path(
  GENE2FUNC_DIR,
  "GS.txt"
)

SUMMARY_PATH <- file.path(
  GENE2FUNC_DIR,
  "summary.txt"
)

GENE_IDS_PATH <- file.path(
  GENE2FUNC_DIR,
  "geneIDs.txt"
)

GENE_TABLE_PATH <- file.path(
  GENE2FUNC_DIR,
  "geneTable.txt"
)

GTEX_GENERAL_LOG2_PATH <- file.path(
  GENE2FUNC_DIR,
  "gtex_v8_ts_general_avg_log2TPM_exp.txt"
)

GTEX_GENERAL_NORM_PATH <- file.path(
  GENE2FUNC_DIR,
  "gtex_v8_ts_general_avg_normTPM_exp.txt"
)

GTEX_GENERAL_DEG_PATH <- file.path(
  GENE2FUNC_DIR,
  "gtex_v8_ts_general_DEG.txt"
)

# By default, place all custom outputs in a clean subdirectory so the
# original FUMA files remain untouched.
OUT_DIR <- if (length(args) >= 2 && nzchar(args[[2]])) {
  normalizePath(
    path.expand(args[[2]]),
    mustWork = FALSE
  )
} else {
  Sys.getenv(
    "OUT_DIR",
    unset = file.path(
      GENE2FUNC_DIR,
      "resilience_enrichment_plots"
    )
  )
}

FDR_THRESHOLD <- 0.05

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

message("")
message("============================================================")
message("FUMA GENE2FUNC resilience enrichment")
message("============================================================")
message("GENE2FUNC input directory:")
message("  ", GENE2FUNC_DIR)
message("GS input:")
message("  ", GS_PATH)
message("Output directory:")
message("  ", OUT_DIR)
message("")

if (!file.exists(GS_PATH)) {
  stop(
    "GS.txt was not found:\n  ",
    GS_PATH,
    "\n\nResolved GENE2FUNC directory:\n  ",
    GENE2FUNC_DIR,
    "\n\nOn the local Mac, confirm that /Users/hao/cubic-home is mounted ",
    "and contains the FUMA output, or pass the directory explicitly."
  )
}

# -----------------------------------------------------------------------------
# 2. READ FUMA ENRICHMENT RESULTS
# -----------------------------------------------------------------------------

gs <- fread(
  GS_PATH,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  na.strings = c("", "NA", "NaN")
)

required_cols <- c(
  "Category",
  "GeneSet",
  "N_genes",
  "N_overlap",
  "p",
  "adjP",
  "genes"
)

missing_cols <- setdiff(required_cols, names(gs))

if (length(missing_cols) > 0) {
  stop(
    "GS.txt is missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# Recover the FUMA input-gene count from summary.txt if available.
N_INPUT <- NA_integer_

if (file.exists(SUMMARY_PATH)) {
  summary_df <- fread(
    SUMMARY_PATH,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    data.table = FALSE
  )

  hit <- summary_df[
    summary_df[[1]] == "Number of input genes",
    ,
    drop = FALSE
  ]

  if (nrow(hit) >= 1) {
    N_INPUT <- suppressWarnings(
      as.integer(hit[[2]][1])
    )
  }
}

# Current analysis contains 40 systemic resilience-associated genes.
# This fallback is only used if summary.txt is unavailable.
if (is.na(N_INPUT)) {
  N_INPUT <- 40L
  warning(
    "Could not read Number of input genes from summary.txt; ",
    "using N_INPUT = 40."
  )
}

# -----------------------------------------------------------------------------
# 3. CLEAN LABELS
# -----------------------------------------------------------------------------

clean_gene_set_label <- function(category, gene_set) {

  out <- gene_set

  # Remove MSigDB-style source prefixes for readable plotting.
  out <- str_remove(
    out,
    regex("^GOBP_", ignore_case = TRUE)
  )

  out <- str_remove(
    out,
    regex("^GOMF_", ignore_case = TRUE)
  )

  out <- str_remove(
    out,
    regex("^GOCC_", ignore_case = TRUE)
  )

  out <- str_remove(
    out,
    regex("^KEGG_", ignore_case = TRUE)
  )

  out <- str_replace_all(
    out,
    "_",
    " "
  )

  # Preserve natural capitalization for GWAS Catalog and positional terms.
  if (!category %in% c(
    "GWAScatalog",
    "Positional_gene_sets"
  )) {
    out <- str_to_sentence(
      str_to_lower(out)
    )
  }

  out
}

# -----------------------------------------------------------------------------
# 4. BIOLOGICAL DOMAIN ASSIGNMENT
# -----------------------------------------------------------------------------
#
# These domains are intended for biological organization of the figure,
# NOT as a new statistical analysis.
#
# Edit this case_when() if you want to regroup terms differently.
# -----------------------------------------------------------------------------

assign_domain <- function(category, gene_set) {

  case_when(

    # Genomic clustering should be separated from biological pathways.
    category == "Positional_gene_sets" ~
      "Genomic positional loci",

    # Lipid / lipoprotein / membrane transport.
    str_detect(
      gene_set,
      regex(
        paste(
          "HIGH_DENSITY_LIPOPROTEIN",
          "LIPID_TRANSPORTER",
          "APOLIPOPROTEIN",
          "CHOLESTEROL",
          "ABC_TRANSPORTERS",
          sep = "|"
        ),
        ignore_case = TRUE
      )
    ) ~
      "Lipid & membrane transport",

    # Brain aging / AD / vascular-brain phenotypes.
    str_detect(
      gene_set,
      regex(
        paste(
          "AMYLOID",
          "COGNITIVE",
          "MICROBLEED",
          "ALZHEIMER",
          "P-TAU",
          "TAU181",
          sep = "|"
        ),
        ignore_case = TRUE
      )
    ) ~
      "Neurodegeneration & cerebrovascular",

    # Cardiovascular and lifespan-related traits.
    str_detect(
      gene_set,
      regex(
        paste(
          "CARDIOVASCULAR",
          "LONGEVITY",
          "PARENTAL",
          sep = "|"
        ),
        ignore_case = TRUE
      )
    ) ~
      "Cardiovascular & longevity",

    # Diet / smoking / exposure phenotypes.
    str_detect(
      gene_set,
      regex(
        paste(
          "MEAT",
          "LAMB",
          "SMOKER",
          "METHYLATION",
          sep = "|"
        ),
        ignore_case = TRUE
      )
    ) ~
      "Diet & exposure",

    TRUE ~
      "Other"
  )
}

domain_order <- c(
  "Lipid & membrane transport",
  "Neurodegeneration & cerebrovascular",
  "Cardiovascular & longevity",
  "Diet & exposure",
  "Genomic positional loci",
  "Other"
)

assign_evidence_class <- function(category) {
  case_when(
    category %in% c(
      "GO_bp",
      "GO_mf",
      "GO_cc",
      "KEGG"
    ) ~ "Functional pathway/process",

    category == "GWAScatalog" ~
      "GWAS trait annotation",

    category == "Positional_gene_sets" ~
      "Positional gene set",

    TRUE ~ "Other"
  )
}

source_label <- function(category) {
  case_when(
    category == "GO_bp" ~ "GO Biological Process",
    category == "GO_mf" ~ "GO Molecular Function",
    category == "GO_cc" ~ "GO Cellular Component",
    category == "KEGG" ~ "KEGG",
    category == "GWAScatalog" ~ "GWAS Catalog",
    category == "Positional_gene_sets" ~ "Positional gene set",
    TRUE ~ category
  )
}

# -----------------------------------------------------------------------------
# 5. BUILD SIGNIFICANT ENRICHMENT TABLE
# -----------------------------------------------------------------------------

sig <- gs %>%
  mutate(
    p = as.numeric(p),
    adjP = as.numeric(adjP),
    N_genes = as.integer(N_genes),
    N_overlap = as.integer(N_overlap)
  ) %>%
  filter(
    is.finite(adjP),
    adjP < FDR_THRESHOLD
  ) %>%
  mutate(
    Biological_domain = mapply(
      assign_domain,
      Category,
      GeneSet,
      USE.NAMES = FALSE
    ),
    Evidence_class = vapply(
      Category,
      assign_evidence_class,
      FUN.VALUE = character(1)
    ),
    Source = vapply(
      Category,
      source_label,
      FUN.VALUE = character(1)
    ),
    Display_term = mapply(
      clean_gene_set_label,
      Category,
      GeneSet,
      USE.NAMES = FALSE
    ),
    minus_log10_FDR = -log10(adjP),
    Input_gene_fraction = N_overlap / N_INPUT
  )

if (nrow(sig) == 0) {
  stop(
    "No gene-set enrichments passed adjP < ",
    FDR_THRESHOLD,
    "."
  )
}

# -----------------------------------------------------------------------------
# 6. FLAG APOE/APOC1-DOMINATED TERMS
# -----------------------------------------------------------------------------

split_gene_string <- function(x) {

  if (
    is.na(x) ||
    !nzchar(x)
  ) {
    return(character(0))
  }

  unique(
    str_split(
      x,
      ":",
      simplify = FALSE
    )[[1]]
  )
}

driver_stats <- lapply(
  sig$genes,
  function(x) {

    genes_vec <- split_gene_string(x)

    n_core <- sum(
      genes_vec %in% c(
        "APOE",
        "APOC1"
      )
    )

    n_other <- sum(
      !genes_vec %in% c(
        "APOE",
        "APOC1"
      )
    )

    data.frame(
      Contains_APOE = "APOE" %in% genes_vec,
      Contains_APOC1 = "APOC1" %in% genes_vec,
      N_APOE_APOC1 = n_core,
      N_non_APOE_APOC1 = n_other,
      APOE_APOC1_only = (
        length(genes_vec) > 0 &&
          n_other == 0
      ),
      stringsAsFactors = FALSE
    )
  }
)

driver_stats <- bind_rows(
  driver_stats
)

sig <- bind_cols(
  sig,
  driver_stats
)

sig <- sig %>%
  mutate(
    Driver_annotation = case_when(
      APOE_APOC1_only ~
        "Overlap only APOE/APOC1",
      Contains_APOE | Contains_APOC1 ~
        "Includes APOE/APOC1 + other genes",
      TRUE ~
        "No APOE/APOC1"
    )
  )

# -----------------------------------------------------------------------------
# 7. ORDER BY BIOLOGICAL DOMAIN, THEN FDR
# -----------------------------------------------------------------------------

sig <- sig %>%
  mutate(
    Biological_domain = factor(
      Biological_domain,
      levels = domain_order
    )
  ) %>%
  arrange(
    Biological_domain,
    adjP,
    desc(N_overlap),
    Display_term
  ) %>%
  mutate(
    Rank_within_domain = ave(
      seq_len(n()),
      Biological_domain,
      FUN = seq_along
    )
  )

# -----------------------------------------------------------------------------
# 8. OUTPUT DETAILED ENRICHMENT TABLE
# -----------------------------------------------------------------------------

output_table <- sig %>%
  select(
    Biological_domain,
    Evidence_class,
    Source,
    Category,
    GeneSet,
    Display_term,
    N_genes,
    N_overlap,
    Input_gene_fraction,
    p,
    adjP,
    minus_log10_FDR,
    genes,
    Contains_APOE,
    Contains_APOC1,
    N_APOE_APOC1,
    N_non_APOE_APOC1,
    APOE_APOC1_only,
    Driver_annotation,
    everything()
  )

fwrite(
  output_table,
  file.path(
    OUT_DIR,
    "enrichment_significant_by_biological_domain.tsv"
  ),
  sep = "\t",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 9. DOMAIN SUMMARY
# -----------------------------------------------------------------------------

domain_summary <- sig %>%
  group_by(
    Biological_domain
  ) %>%
  summarise(
    N_significant_terms = n(),
    Best_adjP = min(
      adjP,
      na.rm = TRUE
    ),
    Max_overlap = max(
      N_overlap,
      na.rm = TRUE
    ),
    N_terms_APOE_APOC1_only = sum(
      APOE_APOC1_only
    ),
    N_terms_containing_APOE = sum(
      Contains_APOE
    ),
    N_terms_containing_APOC1 = sum(
      Contains_APOC1
    ),
    .groups = "drop"
  ) %>%
  arrange(
    Biological_domain
  )

fwrite(
  domain_summary,
  file.path(
    OUT_DIR,
    "enrichment_domain_summary.tsv"
  ),
  sep = "\t",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 10. GENE CONTRIBUTION TABLE
# -----------------------------------------------------------------------------

gene_contributions <- sig %>%
  select(
    Biological_domain,
    Evidence_class,
    Source,
    Display_term,
    adjP,
    genes
  ) %>%
  separate_rows(
    genes,
    sep = ":"
  ) %>%
  filter(
    !is.na(genes),
    nzchar(genes)
  ) %>%
  rename(
    Gene = genes
  ) %>%
  group_by(
    Gene
  ) %>%
  summarise(
    N_significant_terms = n(),
    N_biological_domains = n_distinct(
      Biological_domain
    ),
    Best_adjP = min(
      adjP,
      na.rm = TRUE
    ),
    Domains = paste(
      unique(
        as.character(
          Biological_domain
        )
      ),
      collapse = "; "
    ),
    Terms = paste(
      unique(
        Display_term
      ),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(N_significant_terms),
    Best_adjP,
    Gene
  )

fwrite(
  gene_contributions,
  file.path(
    OUT_DIR,
    "enrichment_gene_contributions.tsv"
  ),
  sep = "\t",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 11. FUNCTIONAL-PATHWAY-ONLY TABLE
# -----------------------------------------------------------------------------

functional_only <- sig %>%
  filter(
    Evidence_class ==
      "Functional pathway/process"
  )

fwrite(
  functional_only,
  file.path(
    OUT_DIR,
    "enrichment_functional_pathways_only.tsv"
  ),
  sep = "\t",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 12. PLOT HELPER
# -----------------------------------------------------------------------------

make_enrichment_plot <- function(
    dat,
    title_text,
    subtitle_text = NULL
) {

  if (nrow(dat) == 0) {
    return(NULL)
  }

  plot_dat <- dat %>%
    arrange(
      Biological_domain,
      adjP,
      desc(N_overlap)
    ) %>%
    mutate(
      Plot_label = str_wrap(
        Display_term,
        width = 50
      )
    )

  # Unique labels are required for factor ordering.
  # Add source only if an accidental duplicate label exists.
  if (
    anyDuplicated(
      plot_dat$Plot_label
    )
  ) {
    plot_dat <- plot_dat %>%
      mutate(
        Plot_label = paste0(
          Plot_label,
          " [",
          Source,
          "]"
        )
      )
  }

  # Reverse global factor levels so the most significant item within each
  # domain appears nearest the top of its facet.
  plot_dat$Plot_label <- factor(
    plot_dat$Plot_label,
    levels = rev(
      plot_dat$Plot_label
    )
  )

  ggplot(
    plot_dat,
    aes(
      x = minus_log10_FDR,
      y = Plot_label
    )
  ) +

    geom_segment(
      aes(
        x = 0,
        xend = minus_log10_FDR,
        yend = Plot_label,
        color = Biological_domain
      ),
      linewidth = 0.6,
      alpha = 0.72,
      show.legend = FALSE
    ) +

    geom_point(
      aes(
        size = N_overlap,
        color = Biological_domain,
        shape = Evidence_class
      ),
      stroke = 0.55
    ) +

    geom_vline(
      xintercept = -log10(
        FDR_THRESHOLD
      ),
      linetype = "dashed",
      linewidth = 0.45
    ) +

    facet_grid(
      Biological_domain ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +

    scale_size_continuous(
      name = "Input genes\nin gene set",
      range = c(
        2.6,
        6.6
      ),
      breaks = sort(
        unique(
          plot_dat$N_overlap
        )
      )
    ) +

    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = expression(
        -log[10](
          "FDR-adjusted P"
        )
      ),
      y = NULL,
      shape = "Evidence class"
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
        size = 9
      ),
      axis.title.x = element_text(
        margin = margin(
          t = 8
        )
      ),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold",
        hjust = 1,
        size = 10
      ),
      panel.spacing.y = grid::unit(
        0.55,
        "lines"
      ),
      legend.position = "right",
      plot.margin = margin(
        8,
        12,
        8,
        8
      )
    )
}

# -----------------------------------------------------------------------------
# 13. ALL SIGNIFICANT TERMS PLOT
# -----------------------------------------------------------------------------

p_all <- make_enrichment_plot(
  sig,
  title_text =
    "Gene-set enrichment of systemic resilience-associated expression genes",
  subtitle_text =
    paste0(
      "FUMA GENE2FUNC; ",
      N_INPUT,
      " input genes; terms ordered by biological domain"
    )
)

all_png <- file.path(
  OUT_DIR,
  "enrichment_by_biological_domain.png"
)

all_pdf <- file.path(
  OUT_DIR,
  "enrichment_by_biological_domain.pdf"
)

all_svg <- file.path(
  OUT_DIR,
  "enrichment_by_biological_domain.svg"
)

ggsave(
  all_png,
  p_all,
  width = 11.5,
  height = 11.5,
  units = "in",
  dpi = 500,
  bg = "white"
)

ggsave(
  all_pdf,
  p_all,
  width = 11.5,
  height = 11.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  all_svg,
  p_all,
  width = 11.5,
  height = 11.5,
  units = "in",
  device = svglite::svglite,
  bg = "white"
)

# -----------------------------------------------------------------------------
# 14. FUNCTIONAL-PATHWAY-ONLY PLOT
# -----------------------------------------------------------------------------

if (nrow(functional_only) > 0) {

  p_functional <- make_enrichment_plot(
    functional_only,
    title_text =
      "Functional enrichment of systemic resilience-associated expression genes",
    subtitle_text =
      "GO and KEGG terms only; GWAS trait annotations and positional loci excluded"
  )

  ggsave(
    file.path(
      OUT_DIR,
      "enrichment_functional_pathways_only.png"
    ),
    p_functional,
    width = 10,
    height = 5.5,
    units = "in",
    dpi = 500,
    bg = "white"
  )

  ggsave(
    file.path(
      OUT_DIR,
      "enrichment_functional_pathways_only.pdf"
    ),
    p_functional,
    width = 10,
    height = 5.5,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    file.path(
      OUT_DIR,
      "enrichment_functional_pathways_only.svg"
    ),
    p_functional,
    width = 10,
    height = 5.5,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

# -----------------------------------------------------------------------------
# 15. CONSOLE SUMMARY
# -----------------------------------------------------------------------------

message("")
message("============================================================")
message("FUMA GENE2FUNC enrichment extraction complete")
message("============================================================")
message("Input genes: ", N_INPUT)
message(
  "Significant enriched terms (adjP < ",
  FDR_THRESHOLD,
  "): ",
  nrow(sig)
)
message("")

message("Significant terms by biological domain:")

print(
  domain_summary,
  row.names = FALSE
)

message("")
message(
  "Terms containing APOE: ",
  sum(sig$Contains_APOE),
  " / ",
  nrow(sig)
)

message(
  "Terms containing APOC1: ",
  sum(sig$Contains_APOC1),
  " / ",
  nrow(sig)
)

message(
  "Terms whose overlap consists only of APOE/APOC1: ",
  sum(sig$APOE_APOC1_only),
  " / ",
  nrow(sig)
)

message("")
message("Outputs written to:")
message("  ", normalizePath(OUT_DIR, mustWork = FALSE))
message("")
message(
  "Main table: enrichment_significant_by_biological_domain.tsv"
)
message(
  "Main plot:  enrichment_by_biological_domain.{png,pdf,svg}"
)
message(
  "Functional-only plot: enrichment_functional_pathways_only.{png,pdf,svg}"
)