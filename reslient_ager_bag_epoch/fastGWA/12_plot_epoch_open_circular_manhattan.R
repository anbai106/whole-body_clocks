#!/usr/bin/env Rscript

# =============================================================================
# Brain-proteomics EPOCH-BAG residual GWAS
# Directional open circular Manhattan plot -- lightweight vector version
#
# Residual phenotype:
#   EPOCH-BAG residual = observed mortality EPOCH
#                        - EPOCH expected from BAG
#
# Direction of the fastGWA effect allele:
#   BETA < 0  -> lower EPOCH-BAG residual -> resilience-like direction
#   BETA > 0  -> higher EPOCH-BAG residual -> vulnerability-like direction
#
# This figure is intentionally optimized for publication/vector editing:
#
#   1) NO outer chromosome circle/ring.
#   2) Chromosome labels are placed INSIDE the open circular plot.
#   3) NO GWAS Catalog word clouds. The open upper-right sector contains only
#      the directional legend.
#   4) Individual sub-threshold/non-significant GWAS points are NOT plotted.
#      Instead, each chromosome is represented by a filled annular sector
#      ("curved rectangle") at the base of the Manhattan plot.
#   5) Genome-wide significant variants are summarized into genomic bins
#      (default 250 kb) and drawn as filled curved signal blocks rather than
#      thousands of individual points. FUMA top lead SNPs remain individually
#      displayed and labeled.
#
# Thus, the full GWAS is still used for:
#   - chromosome/genome layout
#   - lead-SNP BETA matching
#   - genome-wide significant signal summarization
# but the vector output contains far fewer graphical objects.
#
# DEFAULT FULL GWAS INPUT:
# /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/
# fastGWA/output/Brain_proteomics_mortality_clock/EPOCH_BAG_residual/
# organ_pheno_normalized_residualized.fastGWA
#
# DEFAULT FUMA DIRECTORY:
# /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/
# fuma/Brain_proteomics_mortality_clock/EPOCH_BAG_residual
#
# Required R packages:
#   data.table, stringr, ggplot2, patchwork, svglite
#
# =============================================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

###############################################################################
# 1. USER SETTINGS
###############################################################################

root_dir <- Sys.getenv(
  "ROOT_DIR",
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock"
)

clock_name <- "Brain_proteomics_mortality_clock"
analysis_name <- "EPOCH_BAG_residual"

default_gwas_path <- file.path(
  root_dir,
  "fastGWA",
  "output",
  clock_name,
  analysis_name,
  "organ_pheno_normalized_residualized.fastGWA"
)

gwas_path <- Sys.getenv("GWAS_PATH", default_gwas_path)

fuma_dir <- Sys.getenv(
  "FUMA_DIR",
  file.path(root_dir, "fuma", clock_name, analysis_name)
)

loci_path <- file.path(fuma_dir, "GenomicRiskLoci.txt")

# Output prefix.
out_stem <- "EPOCH_BAG_residual_directional_open_circular_manhattan"
out_prefix <- file.path(fuma_dir, out_stem)

# Overwrite existing figure outputs.
overwrite_existing_outputs <- TRUE

# Stop if any FUMA top lead SNP cannot be assigned a finite non-zero beta.
strict_lead_beta <- TRUE

# Genome-wide significance threshold.
p_genome_wide <- 5e-8

# Open-circle geometry.
#
# 296 degrees gives a broad opening in the upper-right part of the figure,
# similar to the reference layout supplied by the user.
theta_start_deg <- 82
sweep_deg <- 296
gap_fraction <- 0.004

# ---------------------------------------------------------------------------
# Lightweight signal representation
# ---------------------------------------------------------------------------
# NO individual P > 5e-8 SNP points are drawn.
#
# Genome-wide significant SNPs are condensed into genomic bins. For each bin,
# the most significant SNP determines the radial height of one filled curved
# block. Smaller bins retain more local detail but create more vector objects.
signal_bin_bp <- 250000L

# If TRUE, write the aggregated genome-wide-significant bins to a TSV.
write_signal_bin_table <- TRUE

# Cap displayed -log10(P), if desired. NULL = determine automatically.
p_cap <- NULL

# ---------------------------------------------------------------------------
# Radial layout
# ---------------------------------------------------------------------------
# Chromosome curved ribbon that replaces the dense background point cloud.
chromosome_track_inner <- 0.62
chromosome_track_outer <- 0.70

# Chromosome labels are INSIDE that ribbon.
chromosome_label_radius <- 0.555

# Genome-wide significance dashed arc.
threshold_radius <- 0.755

# Genome-wide signal blocks extend from threshold_radius outward.
signal_outer_radius <- 1.045

# Lead-SNP labels may extend somewhat farther.
lead_label_max_radius <- 1.205

# Plot limit.
plot_outer_radius <- 1.29

# ---------------------------------------------------------------------------
# Visual settings
# ---------------------------------------------------------------------------
figure_width <- 13
figure_height <- 9
png_dpi <- 500

# Lead-locus colors.
resilience_color <- "#2C8AC4"
vulnerability_color <- "#E04B5A"
unknown_color <- "#777777"

# Genome-wide threshold.
threshold_color <- "#ED1C24"

# Chromosome colors.
chr_order <- c(as.character(1:22), "X")
chr_index <- setNames(seq_along(chr_order), chr_order)

chr_colors <- c(
  "#59B6D9", "#F5A134", "#59C68A", "#F16073", "#A88AD0",
  "#B98B8B", "#EDA4C7", "#A9ADB1", "#C7C75A", "#6FB5D4",
  "#EAB56D", "#75C993", "#B49AD3", "#9E87C8", "#EAA5C7",
  "#A9AFB5", "#C4C64F", "#B29ACB", "#E29AB9", "#78B7D9",
  "#94C995", "#B4BBC0", "#92989E"
)
names(chr_colors) <- chr_order

# Slight transparency for the chromosome ribbon and significant blocks.
chromosome_track_alpha <- 0.86
signal_block_alpha <- 0.88

# Lead marker/label sizes.
lead_point_size_resilience <- 2.25
lead_point_size_vulnerability <- 2.35
lead_label_size <- 2.35

# Optional title. FALSE is recommended if this will be inserted as a panel.
show_title <- FALSE
figure_title <- "Brain proteomics EPOCH-BAG discordance residual GWAS"

# Optional panel tag such as "C"; leave blank for none.
panel_tag <- ""

###############################################################################
# 2. DEPENDENCY CHECK
###############################################################################

required_packages <- c(
  "data.table", "stringr", "ggplot2", "patchwork", "svglite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall with:\ninstall.packages(c(\"",
    paste(missing_packages, collapse = "\", \""),
    "\"))"
  )
}

###############################################################################
# 3. GENERAL HELPERS
###############################################################################

normalize_chr <- function(x) {
  x <- toupper(
    str_remove(
      str_remove(str_trim(as.character(x)), regex("^chr", ignore_case = TRUE)),
      "\\.0$"
    )
  )
  
  x[x == "23"] <- "X"
  x[!x %in% chr_order] <- NA_character_
  x
}

detect_column <- function(dt, candidates, required = TRUE) {
  lookup <- setNames(names(dt), tolower(str_trim(names(dt))))
  
  for (candidate in candidates) {
    key <- tolower(candidate)
    if (key %in% names(lookup)) {
      return(unname(lookup[[key]]))
    }
  }
  
  if (required) {
    stop(
      "Could not identify required column. Accepted names: ",
      paste(candidates, collapse = ", "),
      "\nObserved columns: ",
      paste(names(dt), collapse = ", ")
    )
  }
  
  NULL
}

safe_file_size <- function(path) {
  x <- suppressWarnings(file.info(path)$size)
  if (length(x) == 0L || is.na(x)) return(0)
  x
}

format_n <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

###############################################################################
# 4. READ FULL fastGWA RESULTS
###############################################################################

read_fastgwa <- function(path) {
  if (!file.exists(path) || safe_file_size(path) == 0L) {
    stop("Missing or empty residual fastGWA file: ", path)
  }
  
  message("Reading full fastGWA file: ", path)
  
  dt <- fread(
    path,
    sep = "auto",
    header = TRUE,
    showProgress = TRUE
  )
  
  chr_col <- detect_column(
    dt,
    c("CHR", "CHROM", "CHROMOSOME", "#CHROM")
  )
  
  pos_col <- detect_column(
    dt,
    c("POS", "BP", "POSITION", "BASE_PAIR_LOCATION")
  )
  
  p_col <- detect_column(
    dt,
    c("P", "PVAL", "PVALUE", "P_VALUE")
  )
  
  beta_col <- detect_column(
    dt,
    c("BETA", "B", "EFFECT", "ESTIMATE", "EFFECT_SIZE")
  )
  
  snp_col <- detect_column(
    dt,
    c("SNP", "RSID", "ID", "MARKERNAME"),
    required = FALSE
  )
  
  a1_col <- detect_column(
    dt,
    c("A1", "EA", "EFFECT_ALLELE", "ALLELE1"),
    required = FALSE
  )
  
  a2_col <- detect_column(
    dt,
    c("A2", "NEA", "OTHER_ALLELE", "NON_EFFECT_ALLELE", "ALLELE2"),
    required = FALSE
  )
  
  se_col <- detect_column(
    dt,
    c("SE", "STDERR", "STANDARD_ERROR"),
    required = FALSE
  )
  
  af_col <- detect_column(
    dt,
    c("AF1", "EAF", "AF", "A1FREQ"),
    required = FALSE
  )
  
  keep_cols <- unique(
    c(
      chr_col, pos_col, p_col, beta_col,
      snp_col, a1_col, a2_col, se_col, af_col
    )
  )
  
  keep_cols <- keep_cols[
    !is.na(keep_cols) & nzchar(keep_cols)
  ]
  
  out <- copy(dt[, ..keep_cols])
  
  setnames(out, chr_col, "CHR")
  setnames(out, pos_col, "POS")
  setnames(out, p_col, "P")
  setnames(out, beta_col, "BETA")
  
  if (!is.null(snp_col)) setnames(out, snp_col, "SNP")
  if (!is.null(a1_col)) setnames(out, a1_col, "A1")
  if (!is.null(a2_col)) setnames(out, a2_col, "A2")
  if (!is.null(se_col)) setnames(out, se_col, "SE")
  if (!is.null(af_col)) setnames(out, af_col, "AF1")
  
  out[, CHR := normalize_chr(CHR)]
  out[, POS := suppressWarnings(as.numeric(POS))]
  out[, P := suppressWarnings(as.numeric(P))]
  out[, BETA := suppressWarnings(as.numeric(BETA))]
  
  if ("SE" %in% names(out)) {
    out[, SE := suppressWarnings(as.numeric(SE))]
  }
  
  if ("AF1" %in% names(out)) {
    out[, AF1 := suppressWarnings(as.numeric(AF1))]
  }
  
  out <- out[
    !is.na(CHR) &
      is.finite(POS) &
      POS > 0 &
      is.finite(P) &
      P > 0 &
      P <= 1
  ]
  
  out[, POS := as.integer(POS)]
  
  if (!"SNP" %in% names(out)) {
    out[, SNP := paste0(CHR, ":", POS)]
  } else {
    out[, SNP := str_trim(as.character(SNP))]
  }
  
  if (!"A1" %in% names(out)) out[, A1 := NA_character_]
  if (!"A2" %in% names(out)) out[, A2 := NA_character_]
  if (!"SE" %in% names(out)) out[, SE := NA_real_]
  if (!"AF1" %in% names(out)) out[, AF1 := NA_real_]
  
  out[, CHR_INDEX := unname(chr_index[CHR])]
  setorder(out, CHR_INDEX, POS, P)
  out[, CHR_INDEX := NULL]
  
  out[]
}

###############################################################################
# 5. READ FUMA GENOMIC RISK LOCI
###############################################################################

empty_loci_table <- function() {
  data.table(
    GenomicLocus = integer(),
    LabelSNP = character(),
    CHR = character(),
    POS = integer(),
    P = numeric()
  )
}

read_genomic_risk_loci <- function(path) {
  if (!file.exists(path) || safe_file_size(path) == 0L) {
    stop("Missing or empty FUMA GenomicRiskLoci.txt: ", path)
  }
  
  loci <- fread(
    path,
    sep = "\t",
    colClasses = "character",
    fill = TRUE,
    showProgress = FALSE
  )
  
  if (nrow(loci) == 0L) {
    return(empty_loci_table())
  }
  
  required <- c("GenomicLocus", "rsID", "chr", "pos", "p")
  missing <- setdiff(required, names(loci))
  
  if (length(missing) > 0L) {
    stop(
      "GenomicRiskLoci.txt is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  
  loci[, GenomicLocus := suppressWarnings(as.integer(GenomicLocus))]
  loci[, CHR := normalize_chr(chr)]
  loci[, POS := suppressWarnings(as.integer(as.numeric(pos)))]
  loci[, P := suppressWarnings(as.numeric(p))]
  loci[, LabelSNP := str_trim(as.character(rsID))]
  
  loci <- loci[
    !is.na(GenomicLocus) &
      !is.na(CHR) &
      is.finite(POS) &
      POS > 0 &
      is.finite(P) &
      P > 0 &
      P <= p_genome_wide &
      !is.na(LabelSNP) &
      LabelSNP != "" &
      tolower(LabelSNP) != "nan"
  ]
  
  if (nrow(loci) == 0L) {
    return(empty_loci_table())
  }
  
  # If FUMA includes multiple rows per locus, retain the most significant row.
  setorder(loci, GenomicLocus, P)
  
  unique(loci, by = "GenomicLocus")[
    , .(
      GenomicLocus,
      LabelSNP,
      CHR,
      POS,
      P
    )
  ]
}

###############################################################################
# 6. MATCH FUMA TOP LEAD SNPs TO fastGWA AND ASSIGN DIRECTION
###############################################################################

match_lead_snps_to_gwas <- function(gwas, loci) {
  gwas_snp <- unique(
    gwas[order(P)],
    by = "SNP"
  )
  
  gwas_pos <- unique(
    gwas[order(P)],
    by = c("CHR", "POS")
  )
  
  out <- rbindlist(
    lapply(seq_len(nrow(loci)), function(i) {
      x <- loci[i]
      
      hit <- gwas_snp[SNP == x$LabelSNP]
      source <- "unmatched"
      
      if (nrow(hit) > 0L) {
        source <- "fastGWA SNP ID"
      } else {
        hit <- gwas_pos[
          CHR == x$CHR &
            POS == x$POS
        ]
        
        if (nrow(hit) > 0L) {
          source <- "fastGWA chromosome-position"
        }
      }
      
      if (nrow(hit) > 0L) {
        data.table(
          GenomicLocus = x$GenomicLocus,
          LabelSNP = x$LabelSNP,
          CHR = hit$CHR[1],
          POS = as.integer(hit$POS[1]),
          P = as.numeric(hit$P[1]),
          BETA = as.numeric(hit$BETA[1]),
          SE = as.numeric(hit$SE[1]),
          A1 = as.character(hit$A1[1]),
          A2 = as.character(hit$A2[1]),
          AF1 = as.numeric(hit$AF1[1]),
          MatchSource = source
        )
      } else {
        data.table(
          GenomicLocus = x$GenomicLocus,
          LabelSNP = x$LabelSNP,
          CHR = x$CHR,
          POS = as.integer(x$POS),
          P = as.numeric(x$P),
          BETA = NA_real_,
          SE = NA_real_,
          A1 = NA_character_,
          A2 = NA_character_,
          AF1 = NA_real_,
          MatchSource = source
        )
      }
    })
  )
  
  out[, NEG_LOG10_P := -log10(P)]
  
  out[, Direction := fifelse(
    is.finite(BETA) & BETA < 0,
    "Resilience",
    fifelse(
      is.finite(BETA) & BETA > 0,
      "Vulnerability",
      "Unknown"
    )
  )]
  
  out[, DirectionDefinition := fifelse(
    Direction == "Resilience",
    "lead-SNP beta < 0: lower EPOCH-BAG residual",
    fifelse(
      Direction == "Vulnerability",
      "lead-SNP beta > 0: higher EPOCH-BAG residual",
      "lead-SNP beta unavailable or zero"
    )
  )]
  
  out[, CHR_INDEX := unname(chr_index[CHR])]
  setorder(out, CHR_INDEX, POS, P)
  out[, CHR_INDEX := NULL]
  
  out[]
}

validate_lead_directions <- function(lead_snps) {
  unknown <- lead_snps[
    Direction == "Unknown" |
      !is.finite(BETA) |
      BETA == 0
  ]
  
  if (nrow(unknown) == 0L) {
    return(invisible(TRUE))
  }
  
  msg <- paste0(
    "Could not assign a directional beta to ",
    nrow(unknown),
    " FUMA locus/loci:\n",
    paste(
      paste0(
        "  locus ",
        unknown$GenomicLocus,
        ": ",
        unknown$LabelSNP,
        " (Chr",
        unknown$CHR,
        ":",
        unknown$POS,
        "; match=",
        unknown$MatchSource,
        ")"
      ),
      collapse = "\n"
    )
  )
  
  if (strict_lead_beta) {
    stop(
      msg,
      "\nBecause strict_lead_beta=TRUE, plotting stops rather than ",
      "misclassifying lead-locus direction."
    )
  } else {
    warning(msg)
  }
  
  invisible(FALSE)
}

###############################################################################
# 7. GENOME LAYOUT
###############################################################################

build_genome_layout <- function(gwas) {
  chr_lengths <- gwas[
    ,
    .(
      CHR_LENGTH = max(POS, na.rm = TRUE)
    ),
    by = CHR
  ]
  
  chr_lengths[, CHR_INDEX := unname(chr_index[CHR])]
  setorder(chr_lengths, CHR_INDEX)
  
  total_bp <- sum(chr_lengths$CHR_LENGTH)
  gap_bp <- total_bp * gap_fraction
  
  cursor <- 0
  starts <- numeric(nrow(chr_lengths))
  ends <- numeric(nrow(chr_lengths))
  
  for (i in seq_len(nrow(chr_lengths))) {
    starts[i] <- cursor
    ends[i] <- cursor + chr_lengths$CHR_LENGTH[i]
    cursor <- ends[i] + gap_bp
  }
  
  chr_lengths[, START := starts]
  chr_lengths[, END := ends]
  chr_lengths[, MID := (START + END) / 2]
  
  genome_span <- cursor - gap_bp
  
  chr_lengths[, THETA_START := theta_start_deg +
                sweep_deg * START / genome_span]
  
  chr_lengths[, THETA_END := theta_start_deg +
                sweep_deg * END / genome_span]
  
  chr_lengths[, THETA_MID := theta_start_deg +
                sweep_deg * MID / genome_span]
  
  list(
    chromosomes = chr_lengths,
    genome_span = genome_span
  )
}

add_angles <- function(dt, layout) {
  out <- merge(
    copy(dt),
    layout$chromosomes[
      ,
      .(
        CHR,
        START
      )
    ],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )
  
  out[, CUM_POS := START + POS]
  
  out[, THETA_DEG := theta_start_deg +
        sweep_deg * CUM_POS / layout$genome_span]
  
  out[]
}

###############################################################################
# 8. AGGREGATE GENOME-WIDE SIGNIFICANT SIGNALS INTO CURVED BLOCKS
###############################################################################

build_significant_signal_bins <- function(gwas, layout) {
  sig <- gwas[
    P <= p_genome_wide
  ]
  
  if (nrow(sig) == 0L) {
    return(data.table())
  }
  
  # Bin only genome-wide significant variants.
  sig[, BIN_ID := floor((POS - 1L) / signal_bin_bp)]
  
  bins <- sig[
    ,
    .(
      BIN_START_POS = min(POS),
      BIN_END_POS = max(POS),
      N_SIGNIFICANT_SNPS = .N,
      MIN_P = min(P, na.rm = TRUE),
      TOP_SNP = SNP[which.min(P)][1]
    ),
    by = .(
      CHR,
      BIN_ID
    )
  ]
  
  # Give singleton bins some visible genomic width while respecting chromosome
  # boundaries. The angular width is still based on a genomic interval rather
  # than an individual point.
  bins[, NOMINAL_START := as.integer(BIN_ID * signal_bin_bp + 1L)]
  bins[, NOMINAL_END := as.integer((BIN_ID + 1L) * signal_bin_bp)]
  
  bins <- merge(
    bins,
    layout$chromosomes[
      ,
      .(
        CHR,
        START,
        CHR_LENGTH
      )
    ],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )
  
  bins[, DRAW_START_POS := pmax(1L, NOMINAL_START)]
  bins[, DRAW_END_POS := pmin(CHR_LENGTH, NOMINAL_END)]
  
  bins[, THETA_START := theta_start_deg +
         sweep_deg * (START + DRAW_START_POS) / layout$genome_span]
  
  bins[, THETA_END := theta_start_deg +
         sweep_deg * (START + DRAW_END_POS) / layout$genome_span]
  
  bins[, NEG_LOG10_P := -log10(MIN_P)]
  
  bins[]
}

###############################################################################
# 9. RADIAL TRANSFORM
###############################################################################

determine_p_cap <- function(gwas, lead_snps, signal_bins) {
  candidates <- c(
    -log10(p_genome_wide),
    lead_snps$NEG_LOG10_P
  )
  
  if (nrow(signal_bins) > 0L) {
    candidates <- c(
      candidates,
      signal_bins$NEG_LOG10_P
    )
  }
  
  observed_max <- max(
    candidates[is.finite(candidates)],
    na.rm = TRUE
  )
  
  if (is.null(p_cap)) {
    max(
      10,
      ceiling(observed_max + 1)
    )
  } else {
    as.numeric(p_cap)
  }
}

p_to_radius <- function(neg_log10_p, radial_cap) {
  # The threshold itself maps to threshold_radius.
  threshold_logp <- -log10(p_genome_wide)
  
  denom <- max(
    radial_cap - threshold_logp,
    1e-8
  )
  
  scaled <- (
    pmin(neg_log10_p, radial_cap) -
      threshold_logp
  ) / denom
  
  scaled <- pmax(
    0,
    pmin(1, scaled)
  )
  
  threshold_radius +
    scaled *
    (signal_outer_radius - threshold_radius)
}

###############################################################################
# 10. LABEL PLACEMENT
###############################################################################

spread_labels <- function(
    leads,
    angular_window = 3.6,
    radial_step = 0.038
) {
  out <- copy(leads)
  
  if (nrow(out) == 0L) {
    return(out)
  }
  
  setorder(out, THETA_DEG)
  
  recent_angle <- numeric()
  recent_level <- integer()
  label_levels <- integer(nrow(out))
  
  for (i in seq_len(nrow(out))) {
    theta <- out$THETA_DEG[i]
    
    keep <- (theta - recent_angle) < angular_window
    recent_angle <- recent_angle[keep]
    recent_level <- recent_level[keep]
    
    level <- 0L
    
    while (level %in% recent_level) {
      level <- level + 1L
    }
    
    label_levels[i] <- level
    
    recent_angle <- c(recent_angle, theta)
    recent_level <- c(recent_level, level)
  }
  
  out[, LABEL_LEVEL := label_levels]
  out[, LABEL_OFFSET := 0.045 + LABEL_LEVEL * radial_step]
  
  out[, LABEL_RADIUS := pmin(
    RADIUS + LABEL_OFFSET,
    lead_label_max_radius
  )]
  
  out[]
}

###############################################################################
# 11. CUSTOM LEGEND FOR THE OPEN SECTOR
###############################################################################

make_open_sector_legend <- function(
    n_resilience,
    n_vulnerability,
    n_signal_bins
) {
  ggplot() +
    
    # Resilience locus.
    annotate(
      "point",
      x = 0.07,
      y = 0.79,
      shape = 21,
      size = 4.0,
      fill = resilience_color,
      color = "black",
      stroke = 0.45
    ) +
    annotate(
      "text",
      x = 0.15,
      y = 0.79,
      label = paste0(
        "Resilience loci: lead-SNP beta < 0 (n=",
        n_resilience,
        ")"
      ),
      hjust = 0,
      vjust = 0.5,
      size = 3.15,
      color = resilience_color
    ) +
    
    # Vulnerability locus.
    annotate(
      "point",
      x = 0.07,
      y = 0.59,
      shape = 24,
      size = 4.1,
      fill = vulnerability_color,
      color = "black",
      stroke = 0.45
    ) +
    annotate(
      "text",
      x = 0.15,
      y = 0.59,
      label = paste0(
        "Vulnerability loci: lead-SNP beta > 0 (n=",
        n_vulnerability,
        ")"
      ),
      hjust = 0,
      vjust = 0.5,
      size = 3.15,
      color = vulnerability_color
    ) +
    
    # Genome-wide threshold.
    annotate(
      "segment",
      x = 0.035,
      xend = 0.105,
      y = 0.39,
      yend = 0.39,
      color = threshold_color,
      linewidth = 1.0,
      linetype = "22"
    ) +
    annotate(
      "text",
      x = 0.15,
      y = 0.39,
      label = "Genome-wide: P < 5 x 10^-8",
      hjust = 0,
      vjust = 0.5,
      size = 3.05,
      color = "black"
    ) +
    
    # Aggregated signal block.
    annotate(
      "rect",
      xmin = 0.035,
      xmax = 0.105,
      ymin = 0.17,
      ymax = 0.25,
      fill = "#8FAFC2",
      color = NA
    ) +
    annotate(
      "text",
      x = 0.15,
      y = 0.21,
      label = paste0(
        "Significant SNPs summarized as ",
        signal_bin_bp / 1000,
        "-kb curved blocks (n=",
        n_signal_bins,
        ")"
      ),
      hjust = 0,
      vjust = 0.5,
      size = 2.75,
      color = "black"
    ) +
    
    xlim(0, 1) +
    ylim(0, 1) +
    
    theme_void() +
    
    theme(
      plot.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.55
      ),
      plot.margin = margin(
        7,
        8,
        7,
        8
      )
    )
}

###############################################################################
# 12. BUILD THE OPEN CIRCULAR MANHATTAN PLOT
###############################################################################

make_open_circular_manhattan <- function(
    gwas,
    lead_snps,
    signal_bins
) {
  layout <- build_genome_layout(gwas)
  
  leads <- add_angles(
    lead_snps[
      CHR %in% layout$chromosomes$CHR
    ],
    layout
  )
  
  radial_cap <- determine_p_cap(
    gwas = gwas,
    lead_snps = leads,
    signal_bins = signal_bins
  )
  
  leads[, NEG_LOG10_P := -log10(P)]
  leads[, RADIUS := p_to_radius(
    NEG_LOG10_P,
    radial_cap
  )]
  
  # Ensure very significant lead SNPs are clearly outside the chromosome ribbon.
  leads[, RADIUS := pmax(
    RADIUS,
    threshold_radius
  )]
  
  # Recompute signal-bin geometry with the current layout.
  sig_blocks <- copy(signal_bins)
  
  if (nrow(sig_blocks) > 0L) {
    sig_blocks[, RADIUS := p_to_radius(
      NEG_LOG10_P,
      radial_cap
    )]
    
    sig_blocks[, RADIUS := pmax(
      RADIUS,
      threshold_radius + 0.008
    )]
  }
  
  # ---------------------------------------------------------------------------
  # Chromosome track: one filled curved rectangle per chromosome.
  # This replaces all individual sub-threshold SNP dots.
  # ---------------------------------------------------------------------------
  chr_track <- copy(
    layout$chromosomes
  )
  
  chr_track[, YMIN := chromosome_track_inner]
  chr_track[, YMAX := chromosome_track_outer]
  
  # Chromosome labels INSIDE the inner circle.
  chr_labels <- copy(
    layout$chromosomes
  )
  
  chr_labels[, LABEL_RADIUS := chromosome_label_radius]
  chr_labels[, LABEL := paste0("Chr", CHR)]
  
  # Genome-wide threshold arc.
  threshold_arc <- data.table(
    THETA_DEG = seq(
      theta_start_deg,
      theta_start_deg + sweep_deg,
      length.out = 1800
    ),
    RADIUS = threshold_radius
  )
  
  # Lead label placement.
  label_data <- spread_labels(leads)
  
  resilience_leads <- leads[
    Direction == "Resilience"
  ]
  
  vulnerability_leads <- leads[
    Direction == "Vulnerability"
  ]
  
  unknown_leads <- leads[
    Direction == "Unknown"
  ]
  
  resilience_labels <- label_data[
    Direction == "Resilience"
  ]
  
  vulnerability_labels <- label_data[
    Direction == "Vulnerability"
  ]
  
  unknown_labels <- label_data[
    Direction == "Unknown"
  ]
  
  # ---------------------------------------------------------------------------
  # Start plot.
  # ---------------------------------------------------------------------------
  p <- ggplot()
  
  # ---------------------------------------------------------------------------
  # A. Filled chromosome ribbon (curved rectangles).
  # No individual non-significant SNP points are plotted.
  # ---------------------------------------------------------------------------
  p <- p +
    geom_rect(
      data = chr_track,
      aes(
        xmin = THETA_START,
        xmax = THETA_END,
        ymin = YMIN,
        ymax = YMAX,
        fill = CHR
      ),
      alpha = chromosome_track_alpha,
      color = NA,
      inherit.aes = FALSE
    )
  
  # ---------------------------------------------------------------------------
  # B. Genome-wide significant signal blocks.
  # Each block is a genomic bin and uses the minimum P-value in that bin.
  # ---------------------------------------------------------------------------
  if (nrow(sig_blocks) > 0L) {
    p <- p +
      geom_rect(
        data = sig_blocks,
        aes(
          xmin = THETA_START,
          xmax = THETA_END,
          ymin = threshold_radius,
          ymax = RADIUS,
          fill = CHR
        ),
        alpha = signal_block_alpha,
        color = NA,
        inherit.aes = FALSE
      )
  }
  
  # ---------------------------------------------------------------------------
  # C. Genome-wide threshold.
  # ---------------------------------------------------------------------------
  p <- p +
    geom_path(
      data = threshold_arc,
      aes(
        x = THETA_DEG,
        y = RADIUS
      ),
      color = threshold_color,
      linewidth = 0.90,
      linetype = "22",
      inherit.aes = FALSE
    )
  
  # ---------------------------------------------------------------------------
  # D. Chromosome labels INSIDE the inner circle.
  # ---------------------------------------------------------------------------
  p <- p +
    geom_text(
      data = chr_labels,
      aes(
        x = THETA_MID,
        y = LABEL_RADIUS,
        label = LABEL
      ),
      size = 3.45,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    )
  
  # ---------------------------------------------------------------------------
  # E. Direction-specific FUMA top lead SNPs.
  # ---------------------------------------------------------------------------
  if (nrow(resilience_leads) > 0L) {
    p <- p +
      
      geom_segment(
        data = resilience_leads,
        aes(
          x = THETA_DEG,
          xend = THETA_DEG,
          y = chromosome_track_outer,
          yend = RADIUS
        ),
        linewidth = 1.05,
        color = resilience_color,
        alpha = 0.96,
        lineend = "round",
        inherit.aes = FALSE
      ) +
      
      geom_point(
        data = resilience_leads,
        aes(
          x = THETA_DEG,
          y = RADIUS
        ),
        shape = 21,
        size = lead_point_size_resilience,
        fill = resilience_color,
        color = "black",
        stroke = 0.45,
        inherit.aes = FALSE
      )
  }
  
  if (nrow(vulnerability_leads) > 0L) {
    p <- p +
      
      geom_segment(
        data = vulnerability_leads,
        aes(
          x = THETA_DEG,
          xend = THETA_DEG,
          y = chromosome_track_outer,
          yend = RADIUS
        ),
        linewidth = 1.05,
        color = vulnerability_color,
        alpha = 0.96,
        lineend = "round",
        inherit.aes = FALSE
      ) +
      
      geom_point(
        data = vulnerability_leads,
        aes(
          x = THETA_DEG,
          y = RADIUS
        ),
        shape = 24,
        size = lead_point_size_vulnerability,
        fill = vulnerability_color,
        color = "black",
        stroke = 0.45,
        inherit.aes = FALSE
      )
  }
  
  if (nrow(unknown_leads) > 0L) {
    p <- p +
      
      geom_segment(
        data = unknown_leads,
        aes(
          x = THETA_DEG,
          xend = THETA_DEG,
          y = chromosome_track_outer,
          yend = RADIUS
        ),
        linewidth = 0.9,
        color = unknown_color,
        alpha = 0.9,
        inherit.aes = FALSE
      ) +
      
      geom_point(
        data = unknown_leads,
        aes(
          x = THETA_DEG,
          y = RADIUS
        ),
        shape = 22,
        size = 2.1,
        fill = "white",
        color = unknown_color,
        stroke = 0.55,
        inherit.aes = FALSE
      )
  }
  
  # ---------------------------------------------------------------------------
  # F. Lead-SNP connector lines and labels.
  # ---------------------------------------------------------------------------
  if (nrow(resilience_labels) > 0L) {
    p <- p +
      
      geom_segment(
        data = resilience_labels,
        aes(
          x = THETA_DEG,
          xend = THETA_DEG,
          y = RADIUS,
          yend = LABEL_RADIUS
        ),
        linewidth = 0.34,
        color = resilience_color,
        alpha = 0.92,
        inherit.aes = FALSE
      ) +
      
      geom_text(
        data = resilience_labels,
        aes(
          x = THETA_DEG,
          y = LABEL_RADIUS,
          label = LabelSNP
        ),
        size = lead_label_size,
        fontface = "bold",
        color = resilience_color,
        check_overlap = TRUE,
        inherit.aes = FALSE
      )
  }
  
  if (nrow(vulnerability_labels) > 0L) {
    p <- p +
      
      geom_segment(
        data = vulnerability_labels,
        aes(
          x = THETA_DEG,
          xend = THETA_DEG,
          y = RADIUS,
          yend = LABEL_RADIUS
        ),
        linewidth = 0.34,
        color = vulnerability_color,
        alpha = 0.92,
        inherit.aes = FALSE
      ) +
      
      geom_text(
        data = vulnerability_labels,
        aes(
          x = THETA_DEG,
          y = LABEL_RADIUS,
          label = LabelSNP
        ),
        size = lead_label_size,
        fontface = "bold",
        color = vulnerability_color,
        check_overlap = TRUE,
        inherit.aes = FALSE
      )
  }
  
  if (nrow(unknown_labels) > 0L) {
    p <- p +
      
      geom_segment(
        data = unknown_labels,
        aes(
          x = THETA_DEG,
          xend = THETA_DEG,
          y = RADIUS,
          yend = LABEL_RADIUS
        ),
        linewidth = 0.28,
        color = unknown_color,
        inherit.aes = FALSE
      ) +
      
      geom_text(
        data = unknown_labels,
        aes(
          x = THETA_DEG,
          y = LABEL_RADIUS,
          label = LabelSNP
        ),
        size = 2.1,
        fontface = "bold",
        color = unknown_color,
        check_overlap = TRUE,
        inherit.aes = FALSE
      )
  }
  
  # ---------------------------------------------------------------------------
  # G. Scales / polar conversion.
  #
  # There is deliberately NO outer chromosome circle.
  # ---------------------------------------------------------------------------
  p <- p +
    
    scale_fill_manual(
      values = chr_colors,
      breaks = chr_order,
      drop = FALSE,
      guide = "none"
    ) +
    
    scale_x_continuous(
      limits = c(
        theta_start_deg,
        theta_start_deg + 360
      ),
      expand = expansion(mult = 0)
    ) +
    
    scale_y_continuous(
      limits = c(
        0,
        plot_outer_radius
      ),
      expand = expansion(mult = 0)
    ) +
    
    coord_polar(
      theta = "x",
      start = 0,
      direction = -1,
      clip = "off"
    ) +
    
    theme_void(base_size = 11) +
    
    theme(
      plot.margin = margin(
        14,
        16,
        14,
        14
      ),
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      )
    )
  
  list(
    plot = p,
    layout = layout,
    leads = leads,
    signal_blocks = sig_blocks,
    radial_cap = radial_cap
  )
}

###############################################################################
# 13. MAIN
###############################################################################

message("============================================================")
message("Brain proteomics EPOCH-BAG residual circular Manhattan plot")
message("Lightweight vector / directional lead-locus version")
message("============================================================")
message("GWAS file: ", gwas_path)
message("FUMA directory: ", fuma_dir)
message("FUMA loci file: ", loci_path)

if (!dir.exists(fuma_dir)) {
  stop(
    "FUMA directory not found: ",
    fuma_dir
  )
}

if (!file.exists(gwas_path) || safe_file_size(gwas_path) == 0L) {
  stop(
    "Missing or empty residual fastGWA file: ",
    gwas_path
  )
}

if (!file.exists(loci_path) || safe_file_size(loci_path) == 0L) {
  stop(
    "Missing or empty FUMA loci file: ",
    loci_path
  )
}

# Image paths.
png_path <- paste0(out_prefix, ".png")
pdf_path <- paste0(out_prefix, ".pdf")
svg_path <- paste0(out_prefix, ".svg")

if (
  !overwrite_existing_outputs &&
  any(
    file.exists(
      c(
        png_path,
        pdf_path,
        svg_path
      )
    )
  )
) {
  stop(
    "At least one output figure already exists. ",
    "Set overwrite_existing_outputs <- TRUE to regenerate."
  )
}

# ---------------------------------------------------------------------------
# Load data.
# ---------------------------------------------------------------------------
gwas <- read_fastgwa(gwas_path)

message(
  "Valid GWAS variants retained: ",
  format_n(nrow(gwas))
)

if (nrow(gwas) == 0L) {
  stop(
    "No valid GWAS variants were found."
  )
}

loci <- read_genomic_risk_loci(loci_path)

message(
  "Genome-wide significant FUMA loci retained: ",
  nrow(loci)
)

if (nrow(loci) == 0L) {
  stop(
    "No valid FUMA loci pass P < ",
    format(p_genome_wide, scientific = TRUE)
  )
}

lead_snps <- match_lead_snps_to_gwas(
  gwas,
  loci
)

validate_lead_directions(
  lead_snps
)

n_resilience <- nrow(
  lead_snps[
    Direction == "Resilience"
  ]
)

n_vulnerability <- nrow(
  lead_snps[
    Direction == "Vulnerability"
  ]
)

n_unknown <- nrow(
  lead_snps[
    Direction == "Unknown"
  ]
)

message("Lead-locus direction counts:")
message("  Resilience (beta < 0): ", n_resilience)
message("  Vulnerability (beta > 0): ", n_vulnerability)

if (n_unknown > 0L) {
  message("  Unknown: ", n_unknown)
}

# ---------------------------------------------------------------------------
# Build genome layout and significant signal bins.
# ---------------------------------------------------------------------------
layout_for_bins <- build_genome_layout(
  gwas
)

signal_bins <- build_significant_signal_bins(
  gwas,
  layout_for_bins
)

message(
  "Genome-wide significant variants summarized into ",
  nrow(signal_bins),
  " curved ",
  signal_bin_bp / 1000,
  "-kb signal blocks."
)

# ---------------------------------------------------------------------------
# Write supporting tables.
# ---------------------------------------------------------------------------
lead_all_path <- paste0(
  out_prefix,
  "_lead_loci_all.tsv"
)

lead_resilience_path <- paste0(
  out_prefix,
  "_lead_loci_resilience_beta_negative.tsv"
)

lead_vulnerability_path <- paste0(
  out_prefix,
  "_lead_loci_vulnerability_beta_positive.tsv"
)

fwrite(
  lead_snps,
  lead_all_path,
  sep = "\t"
)

fwrite(
  lead_snps[
    Direction == "Resilience"
  ],
  lead_resilience_path,
  sep = "\t"
)

fwrite(
  lead_snps[
    Direction == "Vulnerability"
  ],
  lead_vulnerability_path,
  sep = "\t"
)

if (
  write_signal_bin_table &&
  nrow(signal_bins) > 0L
) {
  signal_bin_path <- paste0(
    out_prefix,
    "_genomewide_significant_signal_bins.tsv"
  )
  
  fwrite(
    signal_bins,
    signal_bin_path,
    sep = "\t"
  )
}

direction_summary <- data.table(
  Direction = c(
    "Resilience",
    "Vulnerability",
    "Unknown"
  ),
  Lead_beta_definition = c(
    "BETA < 0: effect allele lowers EPOCH-BAG residual",
    "BETA > 0: effect allele raises EPOCH-BAG residual",
    "BETA unavailable or zero"
  ),
  N_FUMA_loci = c(
    n_resilience,
    n_vulnerability,
    n_unknown
  )
)

fwrite(
  direction_summary,
  paste0(
    out_prefix,
    "_direction_summary.tsv"
  ),
  sep = "\t"
)

# ---------------------------------------------------------------------------
# Build main circular plot.
# ---------------------------------------------------------------------------
plot_result <- make_open_circular_manhattan(
  gwas = gwas,
  lead_snps = lead_snps,
  signal_bins = signal_bins
)

manhattan_plot <- plot_result$plot

# ---------------------------------------------------------------------------
# Put only the legend in the open upper-right sector.
# ---------------------------------------------------------------------------
legend_plot <- make_open_sector_legend(
  n_resilience = n_resilience,
  n_vulnerability = n_vulnerability,
  n_signal_bins = nrow(signal_bins)
)

legend_inset <- wrap_elements(
  full = legend_plot,
  clip = FALSE,
  ignore_tag = TRUE
)

final_plot <- manhattan_plot +
  
  inset_element(
    legend_inset,
    left = 0.565,
    bottom = 0.68,
    right = 0.955,
    top = 0.955,
    align_to = "full",
    on_top = TRUE,
    clip = FALSE
  )

# Optional title / panel tag.
annotation_theme <- theme(
  plot.background = element_rect(
    fill = "white",
    color = NA
  )
)

if (show_title) {
  annotation_theme <- annotation_theme +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0,
        margin = margin(
          l = 32,
          b = 2
        )
      )
    )
}

if (nzchar(panel_tag)) {
  annotation_theme <- annotation_theme +
    theme(
      plot.tag = element_text(
        size = 20,
        face = "bold"
      ),
      plot.tag.position = c(
        0.018,
        0.985
      )
    )
}

if (show_title && nzchar(panel_tag)) {
  
  final_plot <- final_plot +
    plot_annotation(
      title = figure_title,
      tag_levels = list(panel_tag),
      theme = annotation_theme
    )
  
} else if (show_title) {
  
  final_plot <- final_plot +
    plot_annotation(
      title = figure_title,
      theme = annotation_theme
    )
  
} else if (nzchar(panel_tag)) {
  
  final_plot <- final_plot +
    plot_annotation(
      tag_levels = list(panel_tag),
      theme = annotation_theme
    )
  
} else {
  
  final_plot <- final_plot +
    plot_annotation(
      theme = annotation_theme
    )
}

# ---------------------------------------------------------------------------
# Save.
# ---------------------------------------------------------------------------
print(final_plot)

ggsave(
  filename = png_path,
  plot = final_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = png_dpi,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = pdf_path,
  plot = final_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = cairo_pdf,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = svg_path,
  plot = final_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = svglite::svglite,
  bg = "white",
  limitsize = FALSE
)

message("")
message("Finished.")
message("  PNG: ", normalizePath(png_path, mustWork = FALSE))
message("  PDF: ", normalizePath(pdf_path, mustWork = FALSE))
message("  SVG: ", normalizePath(svg_path, mustWork = FALSE))
message("  Lead loci: ", normalizePath(lead_all_path, mustWork = FALSE))
message(
  "  Resilience loci: ",
  normalizePath(
    lead_resilience_path,
    mustWork = FALSE
  )
)
message(
  "  Vulnerability loci: ",
  normalizePath(
    lead_vulnerability_path,
    mustWork = FALSE
  )
)

if (
  write_signal_bin_table &&
  exists("signal_bin_path")
) {
  message(
    "  Signal bins: ",
    normalizePath(
      signal_bin_path,
      mustWork = FALSE
    )
  )
}

message("")
message("Plotting note:")
message(
  "  Individual non-significant SNPs are intentionally omitted. ",
  "The chromosome ribbon replaces the dense background cloud."
)
message(
  "  Genome-wide significant SNPs are aggregated into ",
  signal_bin_bp / 1000,
  "-kb curved signal blocks; FUMA top lead SNPs remain individually plotted."
)