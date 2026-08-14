#!/usr/bin/env Rscript

# Generalized open circular Manhattan plots for all EPOCH mortality clocks.
#
# Key changes from the single-clock script:
#   1. Automatically discovers all clock directories with the required files.
#   2. Processes clocks sequentially to limit memory use.
#   3. Removes the manually curated category_patterns mapping.
#   4. Uses normalized GWAS Catalog trait names directly, with optional
#      data-driven canonicalization of near-duplicate labels.
#   5. Shows only the genome-wide significance threshold (P = 5e-8).
#   6. Saves PNG, PDF, SVG, lead-SNP tables, and word-cloud tables per clock.
#
# Required packages:
# install.packages(c(
#   "data.table", "stringr", "ggplot2", "ggwordcloud",
#   "patchwork", "svglite", "scales", "markdown"
# ))

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(ggwordcloud)
  library(patchwork)
  library(svglite)
  library(scales)
  library(markdown)
})

###############################################################################
# User settings
###############################################################################

root_dir <- file.path(
  "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock",
  "mortality_clock"
)

# "all" discovers and runs every eligible clock.
# "selected" runs only the clock names listed in selected_clocks.
run_mode <- "all"
selected_clocks <- c(
  "Brain_proteomics_mortality_clock"
)

# Input and output conventions.
gwas_filename <- "organ_pheno_normalized_residualized.fastGWA.zip"
loci_filename <- "GenomicRiskLoci.txt"
gwascatalog_filename <- "gwascatalog.txt"
out_stem <- "EPOCH_open_circular_R"

# Existing-output behavior. When FALSE, a clock is skipped if any one of
# PNG, PDF, or SVG already exists. Set TRUE to force regeneration.
overwrite_existing_outputs <- FALSE

# Plot settings.
p_genome_wide <- 5e-8
theta_start_deg <- 82
sweep_deg <- 296
max_background_points <- 350000L
point_size <- 0.45
figure_size <- 13
seed <- 2026L
p_cap <- NULL

# GWAS Catalog word-cloud settings.
#
# "auto_trait": recommended. Uses data-driven normalized trait labels and
#               collapses obvious reporting variants without organ-specific
#               hand-written category rules.
# "trait":      uses normalized GWAS Catalog Trait strings directly.
# "none":       omits the word cloud.
wordcloud_mode <- "auto_trait"

# Frequency counting:
# "locus" counts each term at most once per genomic locus (recommended).
# "association" counts unique locus + PMID + trait combinations.
# "row" counts all rows in gwascatalog.txt.
wordcloud_count <- "locus"
wordcloud_min_frequency <- 1L
wordcloud_max_terms <- 45L

# Keep the original trait labels in the output table for traceability.
write_trait_mapping_table <- TRUE

###############################################################################
# Dependency check
###############################################################################

required_packages <- c(
  "data.table", "stringr", "ggplot2", "ggwordcloud",
  "patchwork", "svglite", "scales", "markdown"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(\"",
    paste(missing_packages, collapse = "\", \""),
    "\"), type = \"binary\")"
  )
}

###############################################################################
# Constants
###############################################################################

chr_order <- c(as.character(1:22), "X")
chr_index <- setNames(seq_along(chr_order), chr_order)

chr_colors <- c(
  "#64B5D9", "#F39A38", "#62C58A", "#F05B68", "#A98BC7",
  "#BA8F8F", "#E8A6C4", "#A6A9AD", "#C7C681", "#72AFCF",
  "#E6B976", "#80C995", "#C1A5D5", "#A98BC7", "#EAA9C5",
  "#AAB0B5", "#C1C44C", "#B39ACB", "#E798B5", "#83B6D5",
  "#9DCA9A", "#B9BEC3", "#92979C"
)
names(chr_colors) <- chr_order

###############################################################################
# General helpers
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
    if (key %in% names(lookup)) return(unname(lookup[[key]]))
  }
  
  if (required) {
    stop(
      "Could not identify required column. Accepted names: ",
      paste(candidates, collapse = ", "),
      ". Observed columns: ",
      paste(names(dt), collapse = ", ")
    )
  }
  
  NULL
}

read_table_robust <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  
  if (endsWith(tolower(path), ".zip")) {
    members <- unzip(path, list = TRUE)$Name
    members <- members[
      !endsWith(members, "/") & !startsWith(members, "__MACOSX")
    ]
    
    if (length(members) == 0L) stop("No readable file in ZIP: ", path)
    
    if (length(members) > 1L) {
      warning("Reading first ZIP member: ", members[1])
    }
    
    cmd <- paste("unzip -p", shQuote(path), shQuote(members[1]))
    return(
      fread(
        cmd = cmd,
        sep = "auto",
        header = TRUE,
        showProgress = TRUE
      )
    )
  }
  
  fread(path, sep = "auto", header = TRUE, showProgress = TRUE)
}

clock_title_from_name <- function(clock_name) {
  x <- clock_name
  x <- str_remove(x, regex("_mortality_clock$", ignore_case = TRUE))
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  
  # Preserve common modality capitalization.
  x <- str_replace_all(x, regex("\\bmri\\b", ignore_case = TRUE), "MRI")
  x <- str_replace_all(x, regex("\\bpet\\b", ignore_case = TRUE), "PET")
  
  # Sentence-style title while preserving MRI/PET.
  words <- str_split(x, " ", simplify = TRUE)
  words <- as.character(words[1, ])
  words <- words[words != ""]
  
  words <- vapply(words, function(w) {
    if (w %in% c("MRI", "PET")) return(w)
    str_to_lower(w)
  }, character(1))
  
  if (length(words) > 0L && !words[1] %in% c("MRI", "PET")) {
    words[1] <- str_to_title(words[1])
  }
  
  paste(c(words, "mortality EPOCH"), collapse = " ")
}

###############################################################################
# Input readers
###############################################################################

read_fastgwa <- function(path) {
  dt <- read_table_robust(path)
  
  chr_col <- detect_column(dt, c("CHR", "CHROM", "CHROMOSOME", "#CHROM"))
  pos_col <- detect_column(dt, c("POS", "BP", "POSITION", "BASE_PAIR_LOCATION"))
  p_col <- detect_column(dt, c("P", "PVAL", "PVALUE", "P_VALUE"))
  snp_col <- detect_column(
    dt,
    c("SNP", "RSID", "ID", "MARKERNAME"),
    required = FALSE
  )
  
  keep <- c(chr_col, pos_col, p_col, snp_col)
  keep <- keep[!vapply(keep, is.null, logical(1))]
  out <- copy(dt[, ..keep])
  
  setnames(out, c(chr_col, pos_col, p_col), c("CHR", "POS", "P"))
  if (!is.null(snp_col)) setnames(out, snp_col, "SNP")
  
  out[, CHR := normalize_chr(CHR)]
  out[, POS := suppressWarnings(as.numeric(POS))]
  out[, P := suppressWarnings(as.numeric(P))]
  
  out <- out[
    !is.na(CHR) &
      is.finite(POS) & POS > 0 &
      is.finite(P) & P > 0 & P <= 1
  ]
  
  out[, POS := as.integer(POS)]
  
  if (!"SNP" %in% names(out)) {
    out[, SNP := paste0(CHR, ":", POS)]
  } else {
    out[, SNP := str_trim(as.character(SNP))]
  }
  
  out[, CHR_INDEX := unname(chr_index[CHR])]
  setorder(out, CHR_INDEX, POS, P)
  out[, CHR_INDEX := NULL]
  out[]
}

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
  if (!file.exists(path)) {
    return(empty_loci_table())
  }
  
  file_size <- suppressWarnings(file.info(path)$size)
  if (is.na(file_size) || file_size == 0L) {
    return(empty_loci_table())
  }
  
  loci <- tryCatch(
    fread(
      path,
      sep = "\t",
      colClasses = "character",
      fill = TRUE,
      showProgress = FALSE
    ),
    error = function(e) {
      warning(
        "Could not read FUMA GenomicRiskLoci.txt: ",
        conditionMessage(e)
      )
      data.table()
    }
  )
  
  if (nrow(loci) == 0L) {
    return(empty_loci_table())
  }
  
  required <- c("GenomicLocus", "rsID", "chr", "pos", "p")
  missing <- setdiff(required, names(loci))
  
  if (length(missing) > 0L) {
    warning(
      "GenomicRiskLoci.txt is missing required columns: ",
      paste(missing, collapse = ", ")
    )
    return(empty_loci_table())
  }
  
  loci[, GenomicLocus := suppressWarnings(as.integer(GenomicLocus))]
  loci[, CHR := normalize_chr(chr)]
  loci[, POS := suppressWarnings(as.integer(as.numeric(pos)))]
  loci[, P := suppressWarnings(as.numeric(p))]
  loci[, LabelSNP := str_trim(as.character(rsID))]
  
  # Retain only valid genome-wide significant loci.
  loci <- loci[
    !is.na(GenomicLocus) &
      !is.na(CHR) &
      is.finite(POS) & POS > 0 &
      is.finite(P) & P > 0 & P <= p_genome_wide &
      !is.na(LabelSNP) & LabelSNP != "" &
      tolower(LabelSNP) != "nan"
  ]
  
  if (nrow(loci) == 0L) {
    return(empty_loci_table())
  }
  
  setorder(loci, GenomicLocus, P)
  
  unique(loci, by = "GenomicLocus")[
    , .(GenomicLocus, LabelSNP, CHR, POS, P)
  ]
}

read_gwas_catalog <- function(path) {
  catalog <- fread(
    path,
    sep = "\t",
    colClasses = "character",
    quote = "",
    fill = TRUE
  )
  
  required <- c("GenomicLocus", "Trait")
  missing <- setdiff(required, names(catalog))
  
  if (length(missing) > 0L) {
    stop(
      "Missing columns in gwascatalog.txt: ",
      paste(missing, collapse = ", ")
    )
  }
  
  if (!"PMID" %in% names(catalog)) catalog[, PMID := ""]
  if (!"IndSigSNP" %in% names(catalog)) catalog[, IndSigSNP := ""]
  
  catalog[, GenomicLocus := suppressWarnings(as.integer(GenomicLocus))]
  catalog[, Trait := str_squish(as.character(Trait))]
  catalog[, PMID := str_trim(as.character(PMID))]
  
  catalog[
    !is.na(GenomicLocus) &
      !is.na(Trait) & Trait != "" &
      tolower(Trait) != "nan"
  ]
}

###############################################################################
# Lead-SNP matching
###############################################################################

match_lead_snps_to_gwas <- function(gwas, loci) {
  gwas_snp <- unique(gwas[order(P)], by = "SNP")
  gwas_pos <- unique(gwas[order(P)], by = c("CHR", "POS"))
  
  out <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
    x <- loci[i]
    hit <- gwas_snp[SNP == x$LabelSNP]
    source <- "FUMA locus file"
    
    if (nrow(hit) > 0L) {
      chr <- hit$CHR[1]
      pos <- hit$POS[1]
      p <- hit$P[1]
      source <- "fastGWA SNP ID"
    } else {
      hit <- gwas_pos[CHR == x$CHR & POS == x$POS]
      
      if (nrow(hit) > 0L) {
        chr <- hit$CHR[1]
        pos <- hit$POS[1]
        p <- hit$P[1]
        source <- "fastGWA chromosome-position"
      } else {
        chr <- x$CHR
        pos <- x$POS
        p <- x$P
      }
    }
    
    data.table(
      GenomicLocus = x$GenomicLocus,
      LabelSNP = x$LabelSNP,
      CHR = chr,
      POS = as.integer(pos),
      P = as.numeric(p),
      MatchSource = source
    )
  }))
  
  out[, NEG_LOG10_P := -log10(P)]
  out[, CHR_INDEX := unname(chr_index[CHR])]
  setorder(out, CHR_INDEX, POS, P)
  out[, CHR_INDEX := NULL]
  out[]
}

###############################################################################
# Generalizable GWAS Catalog trait processing
###############################################################################

normalize_trait_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[_/]", " ")
  x <- str_replace_all(x, "[–—]", "-")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  
  # Remove reporting details that generally do not change phenotype meaning.
  x <- str_remove_all(
    x,
    regex(
      paste0(
        "\\s*\\((?:",
        "adjusted[^)]*|",
        "unadjusted[^)]*|",
        "combined[^)]*|",
        "meta-analysis[^)]*|",
        "discovery[^)]*|",
        "replication[^)]*|",
        "female[^)]*|male[^)]*|",
        "women[^)]*|men[^)]*",
        ")\\)"
      ),
      ignore_case = TRUE
    )
  )
  
  x <- str_replace_all(x, "\\s*;\\s*", "; ")
  x <- str_replace_all(x, "\\s*,\\s*", ", ")
  x <- str_squish(x)
  str_trim(x, side = "both")
}

canonicalize_trait_text <- function(x) {
  x <- normalize_trait_text(x)
  key <- tolower(x)
  
  # Remove cohort/analysis qualifiers while preserving the biological trait.
  key <- str_remove_all(
    key,
    regex(
      paste0(
        "\\b(?:",
        "uk biobank|ukbb|biobank japan|bbj|finngen|",
        "european ancestry|african ancestry|asian ancestry|",
        "east asian ancestry|mixed ancestry|trans-ancestry|",
        "genome-wide association study|gwas|",
        "meta analysis|meta-analysis|",
        "quantitative trait|case control|case-control",
        ")\\b"
      ),
      ignore_case = TRUE
    )
  )
  
  # Remove common analysis suffixes but do not replace phenotype words.
  key <- str_remove_all(
    key,
    regex(
      paste0(
        "\\b(?:",
        "adjusted for [^,;]+|",
        "controlling for [^,;]+|",
        "excluding [^,;]+|",
        "interaction with [^,;]+",
        ")\\b"
      ),
      ignore_case = TRUE
    )
  )
  
  key <- str_replace_all(key, "[^a-z0-9%+.-]+", " ")
  key <- str_squish(key)
  key <- str_trim(key)
  
  # Preserve a readable display term based on the cleaned original string.
  display <- x
  display <- str_remove_all(
    display,
    regex(
      "\\b(?:UK Biobank|UKBB|Biobank Japan|BBJ|FinnGen)\\b",
      ignore_case = TRUE
    )
  )
  display <- str_replace_all(display, "\\s*[-,:;]+\\s*$", "")
  display <- str_squish(display)
  
  data.table(
    TraitNormalized = x,
    TraitKey = key,
    TraitDisplayCandidate = display
  )
}

build_wordcloud_frequencies <- function(
    catalog,
    valid_loci,
    mode = wordcloud_mode,
    count_method = wordcloud_count,
    minimum_frequency = wordcloud_min_frequency,
    maximum_terms = wordcloud_max_terms
) {
  x <- catalog[GenomicLocus %in% as.integer(valid_loci)]
  
  if (nrow(x) == 0L) {
    return(list(
      frequencies = data.table(),
      mapping = data.table()
    ))
  }
  
  normalized <- canonicalize_trait_text(x$Trait)
  x <- cbind(x, normalized)
  
  if (mode == "trait") {
    x[, WordCloudTerm := TraitNormalized]
    x[, TraitKey := tolower(TraitNormalized)]
  } else if (mode == "auto_trait") {
    # Choose the shortest readable label for each canonical key.
    representatives <- x[
      TraitKey != "" & !is.na(TraitKey),
      .SD[which.min(nchar(TraitDisplayCandidate))],
      by = TraitKey,
      .SDcols = c("TraitDisplayCandidate")
    ]
    
    setnames(
      representatives,
      "TraitDisplayCandidate",
      "WordCloudTerm"
    )
    
    x <- merge(
      x,
      representatives,
      by = "TraitKey",
      all.x = TRUE,
      sort = FALSE
    )
  } else if (mode == "none") {
    return(list(
      frequencies = data.table(),
      mapping = data.table()
    ))
  } else {
    stop("wordcloud_mode must be 'auto_trait', 'trait', or 'none'.")
  }
  
  x <- x[
    !is.na(WordCloudTerm) & WordCloudTerm != "" &
      !is.na(TraitKey) & TraitKey != ""
  ]
  
  if (nrow(x) == 0L) {
    return(list(
      frequencies = data.table(),
      mapping = data.table()
    ))
  }
  
  counted <- switch(
    count_method,
    locus = unique(x[, .(GenomicLocus, TraitKey, WordCloudTerm)]),
    association = unique(
      x[, .(GenomicLocus, PMID, TraitKey, WordCloudTerm)]
    ),
    row = x[, .(GenomicLocus, PMID, TraitKey, WordCloudTerm)],
    stop("wordcloud_count must be 'locus', 'association', or 'row'.")
  )
  
  freq <- counted[, .(Frequency = .N), by = .(TraitKey, WordCloudTerm)]
  freq <- freq[Frequency >= minimum_frequency]
  setorder(freq, -Frequency, WordCloudTerm)
  freq <- head(freq, maximum_terms)
  
  mapping <- unique(
    x[, .(
      GenomicLocus,
      PMID,
      OriginalTrait = Trait,
      TraitNormalized,
      TraitKey,
      WordCloudTerm
    )]
  )
  
  list(
    frequencies = freq[],
    mapping = mapping[]
  )
}

###############################################################################
# Genome layout
###############################################################################

build_genome_layout <- function(gwas, gap_fraction = 0.004) {
  chr_lengths <- gwas[, .(CHR_LENGTH = max(POS)), by = CHR]
  chr_lengths[, CHR_INDEX := unname(chr_index[CHR])]
  setorder(chr_lengths, CHR_INDEX)
  
  total_bp <- sum(chr_lengths$CHR_LENGTH)
  gap_bp <- total_bp * gap_fraction
  cursor <- 0
  starts <- ends <- numeric(nrow(chr_lengths))
  
  for (i in seq_len(nrow(chr_lengths))) {
    starts[i] <- cursor
    ends[i] <- cursor + chr_lengths$CHR_LENGTH[i]
    cursor <- ends[i] + gap_bp
  }
  
  chr_lengths[, START := starts]
  chr_lengths[, END := ends]
  chr_lengths[, MID := (START + END) / 2]
  
  list(
    chromosomes = chr_lengths,
    genome_span = cursor - gap_bp
  )
}

add_angles <- function(dt, layout) {
  out <- merge(
    copy(dt),
    layout$chromosomes[, .(CHR, START)],
    by = "CHR",
    all.x = TRUE,
    sort = FALSE
  )
  
  out[, CUM_POS := START + POS]
  out[, THETA_DEG := theta_start_deg +
        sweep_deg * CUM_POS / layout$genome_span]
  out[]
}

choose_background_variants <- function(gwas, lead_snps) {
  lead_keys <- paste(lead_snps$CHR, lead_snps$POS, sep = ":")
  x <- copy(gwas)
  x[, IS_LEAD := paste(CHR, POS, sep = ":") %in% lead_keys]
  
  kept <- x[P <= p_genome_wide | IS_LEAD]
  background <- x[!(P <= p_genome_wide | IS_LEAD)]
  
  if (nrow(background) > max_background_points) {
    set.seed(seed)
    background <- background[
      sample.int(nrow(background), max_background_points)
    ]
  }
  
  rbindlist(list(kept, background))[, IS_LEAD := NULL][]
}

spread_labels <- function(leads, window = 3.2, step = 0.033) {
  out <- copy(leads)
  setorder(out, THETA_DEG)
  
  recent_angle <- numeric()
  recent_level <- integer()
  levels <- integer(nrow(out))
  
  for (i in seq_len(nrow(out))) {
    theta <- out$THETA_DEG[i]
    keep <- (theta - recent_angle) < window
    recent_angle <- recent_angle[keep]
    recent_level <- recent_level[keep]
    
    level <- 0L
    while (level %in% recent_level) level <- level + 1L
    levels[i] <- level
    
    recent_angle <- c(recent_angle, theta)
    recent_level <- c(recent_level, level)
  }
  
  out[, LABEL_LEVEL := levels]
  out[, LABEL_OFFSET := 0.035 + LABEL_LEVEL * step]
  out[]
}

###############################################################################
# Plotting
###############################################################################

make_manhattan_plot <- function(gwas, lead_snps) {
  selected <- choose_background_variants(gwas, lead_snps)
  layout <- build_genome_layout(gwas)
  
  selected <- add_angles(selected, layout)
  leads <- add_angles(
    lead_snps[CHR %in% layout$chromosomes$CHR],
    layout
  )
  
  selected[, NEG_LOG10_P := -log10(P)]
  leads[, NEG_LOG10_P := -log10(P)]
  
  observed_max <- max(
    selected$NEG_LOG10_P,
    leads$NEG_LOG10_P,
    -log10(p_genome_wide),
    na.rm = TRUE
  )
  
  radial_cap <- if (is.null(p_cap)) {
    max(10, ceiling(observed_max + 1))
  } else {
    as.numeric(p_cap)
  }
  
  radial_inner <- 0.47
  radial_outer <- 1.04
  chr_radius <- 1.19
  chr_label_radius <- 1.285
  
  p_to_radius <- function(x) {
    radial_inner +
      (pmin(x, radial_cap) / radial_cap) *
      (radial_outer - radial_inner)
  }
  
  selected[, RADIUS := p_to_radius(NEG_LOG10_P)]
  leads[, RADIUS := p_to_radius(NEG_LOG10_P)]
  gws_radius <- p_to_radius(-log10(p_genome_wide))
  
  gws_segments <- selected[
    P <= p_genome_wide,
    .(
      THETA_DEG,
      RADIUS,
      RADIUS_START = radial_inner,
      CHR
    )
  ]
  
  lead_segments <- leads[
    ,
    .(
      THETA_DEG,
      RADIUS,
      RADIUS_START = radial_inner,
      CHR
    )
  ]
  
  threshold_arc <- data.table(
    THETA_DEG = seq(
      theta_start_deg,
      theta_start_deg + sweep_deg,
      length.out = 1800
    ),
    RADIUS = gws_radius
  )
  
  chr_arcs <- rbindlist(
    lapply(seq_len(nrow(layout$chromosomes)), function(i) {
      x <- layout$chromosomes[i]
      
      t1 <- theta_start_deg +
        sweep_deg * x$START / layout$genome_span
      
      t2 <- theta_start_deg +
        sweep_deg * x$END / layout$genome_span
      
      data.table(
        CHR = x$CHR,
        THETA_DEG = seq(t1, t2, length.out = 220),
        RADIUS = chr_radius
      )
    })
  )
  
  chr_labels <- copy(layout$chromosomes)
  chr_labels[, THETA_DEG := theta_start_deg +
               sweep_deg * MID / layout$genome_span]
  chr_labels[, RADIUS := chr_label_radius]
  chr_labels[, LABEL := paste0("Chr", CHR)]
  
  label_data <- spread_labels(leads)
  label_data[
    ,
    LABEL_RADIUS := pmin(
      RADIUS + LABEL_OFFSET,
      chr_radius - 0.035
    )
  ]
  
  ggplot() +
    geom_segment(
      data = gws_segments,
      aes(
        x = THETA_DEG,
        xend = THETA_DEG,
        y = RADIUS_START,
        yend = RADIUS,
        color = CHR
      ),
      linewidth = 0.28,
      alpha = 0.30,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = selected,
      aes(
        x = THETA_DEG,
        y = RADIUS,
        color = CHR
      ),
      shape = 1,
      size = point_size,
      stroke = 0.28,
      alpha = 0.66,
      show.legend = FALSE
    ) +
    geom_path(
      data = threshold_arc,
      aes(x = THETA_DEG, y = RADIUS),
      color = "#E31A1C",
      linewidth = 0.8,
      linetype = "22",
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = lead_segments,
      aes(
        x = THETA_DEG,
        xend = THETA_DEG,
        y = RADIUS_START,
        yend = RADIUS,
        color = CHR
      ),
      linewidth = 0.75,
      alpha = 0.95,
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = leads,
      aes(
        x = THETA_DEG,
        y = RADIUS,
        color = CHR
      ),
      shape = 21,
      fill = "white",
      size = 1.65,
      stroke = 0.75,
      show.legend = FALSE
    ) +
    geom_path(
      data = chr_arcs,
      aes(
        x = THETA_DEG,
        y = RADIUS,
        group = CHR
      ),
      color = "#A6A9AD",
      linewidth = 3.3,
      lineend = "butt",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = chr_labels,
      aes(
        x = THETA_DEG,
        y = RADIUS,
        label = LABEL
      ),
      size = 3.25,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = label_data,
      aes(
        x = THETA_DEG,
        xend = THETA_DEG,
        y = RADIUS,
        yend = LABEL_RADIUS
      ),
      linewidth = 0.18,
      color = "#4C4C4C",
      alpha = 0.75,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_data,
      aes(
        x = THETA_DEG,
        y = LABEL_RADIUS,
        label = LabelSNP
      ),
      size = 2.10,
      fontface = "bold",
      color = "black",
      check_overlap = TRUE,
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      values = chr_colors,
      breaks = chr_order,
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(theta_start_deg, theta_start_deg + 360),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      limits = c(0, chr_label_radius + 0.08),
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
      plot.margin = margin(18, 18, 18, 18),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

make_wordcloud_plot <- function(freq) {
  set.seed(seed)
  
  cols <- hue_pal(
    h = c(15, 375),
    c = 80,
    l = 55
  )(nrow(freq))
  names(cols) <- freq$WordCloudTerm
  
  ggplot(
    freq,
    aes(
      label = WordCloudTerm,
      size = Frequency,
      color = WordCloudTerm
    )
  ) +
    geom_text_wordcloud(
      rm_outside = TRUE,
      eccentricity = 0.60,
      area_corr = TRUE,
      grid_size = 12,
      seed = seed
    ) +
    scale_size_area(max_size = 12) +
    scale_color_manual(values = cols, guide = "none") +
    labs(
      title = "GWAS Catalog traits associated with lead loci"
    ) +
    theme_void(base_size = 9) +
    theme(
      plot.tag = element_blank(),
      plot.title = element_text(
        hjust = 0.5,
        size = 8.5,
        margin = margin(b = 4)
      ),
      plot.background = element_rect(
        fill = "white",
        color = "#B6B6B6",
        linewidth = 0.6
      ),
      plot.margin = margin(6, 6, 6, 6)
    )
}

make_threshold_legend <- function() {
  ggplot() +
    annotate(
      "segment",
      x = 0.05,
      xend = 0.22,
      y = 0.50,
      yend = 0.50,
      color = "#E31A1C",
      linewidth = 0.9,
      linetype = "22"
    ) +
    annotate(
      "text",
      x = 0.27,
      y = 0.50,
      label = "Genome-wide: P < 5 x 10^-8",
      hjust = 0,
      vjust = 0.5,
      size = 3.0
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )
}

###############################################################################
# Clock discovery and processing
###############################################################################

discover_clocks <- function(root_dir) {
  gwas_root <- file.path(root_dir, "fastGWA", "output")
  
  if (!dir.exists(gwas_root)) {
    stop("GWAS root not found: ", gwas_root)
  }
  
  clocks <- basename(
    list.dirs(
      gwas_root,
      recursive = FALSE,
      full.names = TRUE
    )
  )
  
  sort(unique(clocks[nzchar(clocks)]))
}

clock_result <- function(clock_name, status, message = "") {
  data.table(
    Clock = clock_name,
    Status = status,
    Message = message
  )
}

process_clock <- function(clock_name) {
  figure_title <- clock_title_from_name(clock_name)
  
  gwas_path <- file.path(
    root_dir,
    "fastGWA",
    "output",
    clock_name,
    gwas_filename
  )
  
  fuma_dir <- file.path(
    root_dir,
    "fuma",
    clock_name
  )
  
  loci_path <- file.path(fuma_dir, loci_filename)
  gwascatalog_path <- file.path(fuma_dir, gwascatalog_filename)
  out_prefix <- file.path(fuma_dir, out_stem)
  
  png_path <- paste0(out_prefix, ".png")
  pdf_path <- paste0(out_prefix, ".pdf")
  svg_path <- paste0(out_prefix, ".svg")
  image_outputs <- c(png_path, pdf_path, svg_path)
  
  message("\n============================================================")
  message("Clock: ", clock_name)
  message("Title: ", figure_title)
  message("============================================================")
  
  # Skip before reading millions of GWAS rows when any final image exists.
  existing_outputs <- image_outputs[file.exists(image_outputs)]
  if (!overwrite_existing_outputs && length(existing_outputs) > 0L) {
    msg <- paste0(
      "Existing figure output found; skipped: ",
      paste(basename(existing_outputs), collapse = ", ")
    )
    message(msg)
    return(clock_result(clock_name, "skipped_existing_output", msg))
  }
  
  # Robust input checks. Missing GWAS or FUMA inputs are skipped, not failed.
  if (!file.exists(gwas_path)) {
    msg <- paste0("Missing GWAS file: ", gwas_path)
    message(msg)
    return(clock_result(clock_name, "skipped_missing_gwas", msg))
  }
  
  if (!dir.exists(fuma_dir)) {
    msg <- paste0("Missing FUMA directory: ", fuma_dir)
    message(msg)
    return(clock_result(clock_name, "skipped_missing_fuma", msg))
  }
  
  missing_fuma_files <- c(loci_path, gwascatalog_path)[
    !file.exists(c(loci_path, gwascatalog_path))
  ]
  
  if (length(missing_fuma_files) > 0L) {
    msg <- paste0(
      "Missing FUMA file(s): ",
      paste(basename(missing_fuma_files), collapse = ", ")
    )
    message(msg)
    return(clock_result(clock_name, "skipped_missing_fuma", msg))
  }
  
  empty_fuma_files <- c(loci_path, gwascatalog_path)[
    vapply(
      c(loci_path, gwascatalog_path),
      function(path) {
        size <- suppressWarnings(file.info(path)$size)
        is.na(size) || size == 0L
      },
      logical(1)
    )
  ]
  
  if (length(empty_fuma_files) > 0L) {
    msg <- paste0(
      "Empty FUMA file(s): ",
      paste(basename(empty_fuma_files), collapse = ", ")
    )
    message(msg)
    return(clock_result(clock_name, "skipped_empty_fuma", msg))
  }
  
  message("Reading FUMA loci: ", loci_path)
  loci <- read_genomic_risk_loci(loci_path)
  message(
    "Retained ", nrow(loci),
    " loci at P < ", format(p_genome_wide, scientific = TRUE), "."
  )
  
  if (nrow(loci) == 0L) {
    msg <- paste0(
      "No valid FUMA loci passed P < ",
      format(p_genome_wide, scientific = TRUE),
      "; figure not generated."
    )
    message(msg)
    rm(loci)
    gc(verbose = FALSE)
    return(clock_result(clock_name, "skipped_no_genome_wide_loci", msg))
  }
  
  # Confirm that the GWAS Catalog file is readable before loading the GWAS.
  catalog <- tryCatch(
    read_gwas_catalog(gwascatalog_path),
    error = function(e) e
  )
  
  if (inherits(catalog, "error")) {
    msg <- paste0(
      "Could not read gwascatalog.txt: ",
      conditionMessage(catalog)
    )
    message(msg)
    rm(loci, catalog)
    gc(verbose = FALSE)
    return(clock_result(clock_name, "skipped_invalid_fuma", msg))
  }
  
  if (nrow(catalog) == 0L && wordcloud_mode != "none") {
    msg <- "gwascatalog.txt contains no valid annotations; figure not generated."
    message(msg)
    rm(loci, catalog)
    gc(verbose = FALSE)
    return(clock_result(clock_name, "skipped_no_gwascatalog_annotations", msg))
  }
  
  message("Reading fastGWA results: ", gwas_path)
  gwas <- read_fastgwa(gwas_path)
  message("Retained ", format(nrow(gwas), big.mark = ","), " variants.")
  
  if (nrow(gwas) == 0L) {
    msg <- "GWAS file contains no valid variants; figure not generated."
    message(msg)
    rm(gwas, loci, catalog)
    gc(verbose = FALSE)
    return(clock_result(clock_name, "skipped_no_valid_gwas_variants", msg))
  }
  
  dir.create(fuma_dir, recursive = TRUE, showWarnings = FALSE)
  
  lead_snps <- match_lead_snps_to_gwas(gwas, loci)
  
  lead_path <- paste0(out_prefix, "_lead_snps.tsv")
  fwrite(lead_snps, lead_path, sep = "\t")
  
  wordcloud_result <- build_wordcloud_frequencies(
    catalog,
    loci$GenomicLocus
  )
  
  wordcloud_freq <- wordcloud_result$frequencies
  trait_mapping <- wordcloud_result$mapping
  
  if (nrow(wordcloud_freq) > 0L) {
    fwrite(
      wordcloud_freq,
      paste0(out_prefix, "_wordcloud_frequencies.tsv"),
      sep = "\t"
    )
  }
  
  if (write_trait_mapping_table && nrow(trait_mapping) > 0L) {
    fwrite(
      trait_mapping,
      paste0(out_prefix, "_gwascatalog_trait_mapping.tsv"),
      sep = "\t"
    )
  }
  
  manhattan_plot <- make_manhattan_plot(gwas, lead_snps)
  
  threshold_plot <- make_threshold_legend()
  threshold_inset <- wrap_elements(
    full = threshold_plot,
    clip = FALSE,
    ignore_tag = TRUE
  )
  
  final_plot <- manhattan_plot +
    inset_element(
      threshold_inset,
      left = 0.59,
      bottom = 0.47,
      right = 0.89,
      top = 0.535,
      align_to = "full",
      on_top = TRUE,
      clip = FALSE
    )
  
  if (wordcloud_mode != "none" && nrow(wordcloud_freq) > 0L) {
    wordcloud_plot <- make_wordcloud_plot(wordcloud_freq)
    
    wordcloud_inset <- wrap_elements(
      full = wordcloud_plot,
      clip = FALSE,
      ignore_tag = TRUE
    )
    
    final_plot <- final_plot +
      inset_element(
        wordcloud_inset,
        left = 0.585,
        bottom = 0.57,
        right = 0.925,
        top = 0.855,
        align_to = "full",
        on_top = TRUE,
        clip = FALSE
      )
  }
  
  final_plot <- final_plot +
    plot_annotation(
      title = figure_title,
      tag_levels = "a",
      theme = theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0,
          margin = margin(l = 42, b = 3)
        ),
        plot.tag = element_text(size = 20, face = "bold"),
        plot.tag.position = c(0.025, 0.985),
        plot.background = element_rect(fill = "white", color = NA)
      )
    )
  
  print(final_plot)
  
  ggsave(
    png_path,
    final_plot,
    width = figure_size,
    height = figure_size,
    units = "in",
    dpi = 500,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    pdf_path,
    final_plot,
    width = figure_size,
    height = figure_size,
    units = "in",
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    svg_path,
    final_plot,
    width = figure_size,
    height = figure_size,
    units = "in",
    device = svglite::svglite,
    bg = "white",
    limitsize = FALSE
  )
  
  message("Wrote PNG: ", normalizePath(png_path, mustWork = FALSE))
  message("Wrote PDF: ", normalizePath(pdf_path, mustWork = FALSE))
  message("Wrote SVG: ", normalizePath(svg_path, mustWork = FALSE))
  
  objects_to_remove <- intersect(
    c(
      "gwas", "loci", "lead_snps", "catalog",
      "wordcloud_result", "wordcloud_freq", "trait_mapping",
      "manhattan_plot", "threshold_plot", "threshold_inset",
      "wordcloud_plot", "wordcloud_inset", "final_plot"
    ),
    ls()
  )
  
  rm(list = objects_to_remove)
  gc(verbose = FALSE)
  
  clock_result(clock_name, "success", "")
}

###############################################################################
# Main
###############################################################################

all_discovered_clocks <- discover_clocks(root_dir)

if (run_mode == "all") {
  clocks_to_run <- all_discovered_clocks
} else if (run_mode == "selected") {
  clocks_to_run <- unique(selected_clocks)
} else {
  stop("run_mode must be 'all' or 'selected'.")
}

if (length(clocks_to_run) == 0L) {
  stop("No clock directories were found or selected.")
}

message("Clock directories discovered: ", length(all_discovered_clocks))
message("Clocks to evaluate: ", length(clocks_to_run))
message(paste0("  - ", clocks_to_run, collapse = "\n"))

results <- lapply(clocks_to_run, function(clock_name) {
  tryCatch(
    process_clock(clock_name),
    error = function(e) {
      msg <- conditionMessage(e)
      message("ERROR for ", clock_name, ": ", msg)
      gc(verbose = FALSE)
      clock_result(clock_name, "failed", msg)
    }
  )
})

run_summary <- rbindlist(results, fill = TRUE)

summary_dir <- file.path(root_dir, "fuma")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(
  summary_dir,
  "EPOCH_open_circular_R_run_summary.tsv"
)

fwrite(run_summary, summary_path, sep = "\t")

message("\nRun summary:")
print(run_summary)
message("Wrote run summary: ", normalizePath(summary_path, mustWork = FALSE))