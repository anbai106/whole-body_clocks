#!/usr/bin/env Rscript

# =============================================================================
# Head-to-head cumulative C-index plot:
# EPOCH vs BAG, proteomics + metabolomics clocks
#
# This script compares cumulative survival prediction curves for:
#   1) 11 proteomics + 4 metabolomics mortality EPOCH clocks
#   2) matched 11 ProtBAG + 4 MetBAG clocks
#
# Key features:
#   - Default metric: 5-fold cross-validated C-index
#   - X-axis labels follow the EPOCH mortality-clock cumulative order
#   - BAG curve is plotted as a dotted line
#   - Step-wise significance is shown using symbols only:
#       * P<0.05, ** P<0.01, *** P<0.001
#
# Example:
#   DISEASE_CODE=I10 DISEASE_NAME="Hypertension" Rscript plot_epoch_vs_bag_cumulative_cindex.R
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# User-editable settings
# -----------------------------------------------------------------------------

DISEASE_CODE <- Sys.getenv("DISEASE_CODE", unset = "I500")

DISEASE_NAME <- Sys.getenv(
  "DISEASE_NAME",
  unset = "CHF"
)

EPOCH_INPUT_DIR <- Sys.getenv(
  "EPOCH_INPUT_DIR",
  unset = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_cumulative_EPOCH_PM/disease_free_cv"
)

BAG_INPUT_DIR <- Sys.getenv(
  "BAG_INPUT_DIR",
  unset = "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/output_cumulative_BAG_PM/disease_free_cv"
)

EPOCH_INPUT_FILE <- Sys.getenv(
  "EPOCH_INPUT_FILE",
  unset = file.path(
    EPOCH_INPUT_DIR,
    paste0("cox_cumulative_EPOCH_PM_", DISEASE_CODE, ".tsv")
  )
)

BAG_INPUT_FILE <- Sys.getenv(
  "BAG_INPUT_FILE",
  unset = file.path(
    BAG_INPUT_DIR,
    paste0("cox_cumulative_BAG_PM_", DISEASE_CODE, ".tsv")
  )
)

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(
    dirname(EPOCH_INPUT_DIR),
    "figures_EPOCH_vs_BAG_PM"
  )
)

# Options:
#   "apparent" = plot apparent/in-sample Harrell C-index
#   "cv"       = plot 5-fold cross-validated Harrell C-index
#   "both"     = facet apparent and 5-fold CV curves
PLOT_METRIC <- tolower(Sys.getenv("PLOT_METRIC", unset = "cv"))

# Whether to label the final points.
SHOW_FINAL_LABELS <- TRUE

# Whether to label each point with C-index values.
SHOW_POINT_CINDEX_LABELS <- FALSE

# Whether to show step-wise significance symbols.
# Symbols are based on sequential likelihood-ratio tests comparing each
# cumulative model with the immediately preceding cumulative model.
SHOW_STEPWISE_SIGNIFICANCE <- TRUE

# Output dimensions.
PDF_WIDTH <- 13.5
PDF_HEIGHT <- 7.2
PNG_DPI <- 320

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "ggplot2",
  "scales",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Please install these R packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(tibble)
})

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

normalize_icd <- function(x) {
  x <- toupper(trimws(as.character(x)))
  gsub("\\.", "", x)
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

to_logical_safe <- function(x) {
  x_chr <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("true", "t", "1", "yes", "y") ~ TRUE,
    x_chr %in% c("false", "f", "0", "no", "n") ~ FALSE,
    TRUE ~ NA
  )
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

significance_symbol <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

pretty_score_label <- function(x) {
  dplyr::case_when(
    x == "EPOCH" ~ "Mortality EPOCH",
    x == "BAG" ~ "Matched BAG",
    TRUE ~ x
  )
}

pretty_epoch_clock_label <- function(x) {
  x <- as.character(x)

  out <- dplyr::case_when(
    x == "BASE" ~ "Clinical\nbaseline",
    TRUE ~ x
  )

  out <- gsub("_proteomics$", "\nProt", out, ignore.case = TRUE)
  out <- gsub("_metabolomics$", "\nMet", out, ignore.case = TRUE)

  out <- gsub("Reproductive_female", "Reprod. female", out, ignore.case = TRUE)
  out <- gsub("Reproductive_male", "Reprod. male", out, ignore.case = TRUE)

  out <- gsub("_", " ", out)
  out
}

read_cumulative_file <- function(path, score_type, disease_code_clean) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  dat <- readr::read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    na = c("", "NA", "NaN", "nan")
  )

  required_core <- c(
    "disease_id",
    "N",
    "N_case",
    "N_noncase",
    "event_rate",
    "median_followup_years",
    "cumulative_step",
    "added_clock",
    "n_clocks",
    "c_index",
    "base_c_index",
    "delta_c_index_vs_base",
    "delta_c_index_vs_previous",
    "sequential_lr_p_vs_previous",
    "status"
  )

  missing_core <- setdiff(required_core, names(dat))
  if (length(missing_core) > 0L) {
    stop(
      score_type,
      " input file is missing required columns: ",
      paste(missing_core, collapse = ", "),
      call. = FALSE
    )
  }

  dat <- dat %>%
    mutate(
      disease_id_clean = normalize_icd(disease_id),
      score_type = score_type,
      score_label = pretty_score_label(score_type),
      cumulative_step = as.integer(cumulative_step),
      added_clock = as.character(added_clock),
      n_clocks = as.integer(n_clocks),
      N = as.integer(N),
      N_case = as.integer(N_case),
      N_noncase = as.integer(N_noncase),
      event_rate = to_numeric_safe(event_rate),
      median_followup_years = to_numeric_safe(median_followup_years),
      c_index = to_numeric_safe(c_index),
      base_c_index = to_numeric_safe(base_c_index),
      delta_c_index_vs_base = to_numeric_safe(delta_c_index_vs_base),
      delta_c_index_vs_previous = to_numeric_safe(delta_c_index_vs_previous),
      sequential_lr_p_vs_previous = to_numeric_safe(sequential_lr_p_vs_previous),
      status = as.character(status)
    ) %>%
    filter(disease_id_clean == disease_code_clean) %>%
    arrange(cumulative_step)

  if (nrow(dat) == 0L) {
    stop(
      "Disease code ",
      disease_code_clean,
      " was not found in ",
      score_type,
      " file: ",
      path,
      call. = FALSE
    )
  }

  if (anyDuplicated(dat$cumulative_step)) {
    stop(
      "More than one row exists for at least one cumulative step in ",
      score_type,
      " file.",
      call. = FALSE
    )
  }

  dat
}

make_metric_dataset <- function(dat, metric) {
  metric <- tolower(metric)

  if (metric == "apparent") {
    out <- dat %>%
      transmute(
        disease_id = disease_id_clean,
        score_type,
        score_label,
        metric = "apparent",
        metric_label = "Apparent C-index",
        cumulative_step,
        added_clock,
        n_clocks,
        N,
        N_case,
        N_noncase,
        event_rate,
        median_followup_years,
        c_index_plot = c_index,
        base_c_index_plot = base_c_index,
        delta_vs_base_plot = delta_c_index_vs_base,
        delta_vs_previous_plot = delta_c_index_vs_previous,
        sequential_p = sequential_lr_p_vs_previous,
        status
      )
    return(out)
  }

  if (metric == "cv") {
    cv_required <- c(
      "cv_c_index",
      "cv_base_c_index",
      "delta_cv_c_index_vs_base",
      "delta_cv_c_index_vs_previous",
      "cv_status"
    )

    missing_cv <- setdiff(cv_required, names(dat))
    if (length(missing_cv) > 0L) {
      warning(
        "Skipping CV metric for ",
        unique(dat$score_type),
        " because these columns are missing: ",
        paste(missing_cv, collapse = ", "),
        call. = FALSE
      )
      return(tibble())
    }

    out <- dat %>%
      mutate(
        cv_c_index = to_numeric_safe(cv_c_index),
        cv_base_c_index = to_numeric_safe(cv_base_c_index),
        delta_cv_c_index_vs_base = to_numeric_safe(delta_cv_c_index_vs_base),
        delta_cv_c_index_vs_previous = to_numeric_safe(delta_cv_c_index_vs_previous),
        cv_status = as.character(cv_status)
      ) %>%
      transmute(
        disease_id = disease_id_clean,
        score_type,
        score_label,
        metric = "cv",
        metric_label = "5-fold CV C-index",
        cumulative_step,
        added_clock,
        n_clocks,
        N,
        N_case,
        N_noncase,
        event_rate,
        median_followup_years,
        c_index_plot = cv_c_index,
        base_c_index_plot = cv_base_c_index,
        delta_vs_base_plot = delta_cv_c_index_vs_base,
        delta_vs_previous_plot = delta_cv_c_index_vs_previous,
        sequential_p = sequential_lr_p_vs_previous,
        status = cv_status
      )
    return(out)
  }

  stop("Unknown metric: ", metric, call. = FALSE)
}

make_epoch_x_labels <- function(epoch_raw) {
  epoch_labels <- epoch_raw %>%
    arrange(cumulative_step) %>%
    select(cumulative_step, added_clock) %>%
    mutate(
      x_label = pretty_epoch_clock_label(added_clock)
    )

  label_vec <- epoch_labels$x_label
  names(label_vec) <- as.character(epoch_labels$cumulative_step)

  label_vec
}

# -----------------------------------------------------------------------------
# Read EPOCH and BAG files
# -----------------------------------------------------------------------------

if (!PLOT_METRIC %in% c("apparent", "cv", "both")) {
  stop(
    "PLOT_METRIC must be one of: apparent, cv, both.",
    call. = FALSE
  )
}

disease_code_clean <- normalize_icd(DISEASE_CODE)

epoch_raw <- read_cumulative_file(
  path = EPOCH_INPUT_FILE,
  score_type = "EPOCH",
  disease_code_clean = disease_code_clean
)

bag_raw <- read_cumulative_file(
  path = BAG_INPUT_FILE,
  score_type = "BAG",
  disease_code_clean = disease_code_clean
)

epoch_x_labels <- make_epoch_x_labels(epoch_raw)

metrics_to_plot <- if (PLOT_METRIC == "both") {
  c("apparent", "cv")
} else {
  PLOT_METRIC
}

plot_dat <- bind_rows(
  lapply(metrics_to_plot, function(m) make_metric_dataset(epoch_raw, m)),
  lapply(metrics_to_plot, function(m) make_metric_dataset(bag_raw, m))
) %>%
  filter(
    is.finite(c_index_plot),
    tolower(status) %in% c("ok", "success")
  ) %>%
  mutate(
    score_label = factor(
      score_label,
      levels = c("Mortality EPOCH", "Matched BAG")
    ),
    metric_label = factor(
      metric_label,
      levels = c("Apparent C-index", "5-fold CV C-index")
    ),
    stepwise_sig = significance_symbol(sequential_p),
    cindex_label = sprintf("%.3f", c_index_plot),
    delta_label = sprintf("%+.4f", delta_vs_base_plot)
  )

if (nrow(plot_dat) == 0L) {
  stop("No valid rows are available for plotting.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Check whether EPOCH and BAG used the same sample
# -----------------------------------------------------------------------------

sample_info <- plot_dat %>%
  distinct(
    score_label,
    metric_label,
    N,
    N_case,
    N_noncase,
    event_rate,
    median_followup_years
  )

sample_info_simple <- sample_info %>%
  distinct(score_label, N, N_case, N_noncase, event_rate, median_followup_years)

same_sample <- nrow(
  sample_info_simple %>%
    distinct(N, N_case, N_noncase, event_rate, median_followup_years)
) == 1L

if (!same_sample) {
  warning(
    "EPOCH and BAG do not appear to have identical sample summaries. ",
    "The plot is still generated, but head-to-head interpretation should note this.",
    call. = FALSE
  )
}

sample_caption <- if (same_sample) {
  s <- sample_info_simple %>% slice(1)
  paste0(
    "N = ", scales::comma(s$N),
    "; cases/noncases = ", scales::comma(s$N_case), "/",
    scales::comma(s$N_noncase),
    "; event rate = ", scales::percent(s$event_rate, accuracy = 0.1),
    "; median follow-up = ", sprintf("%.2f", s$median_followup_years),
    " years."
  )
} else {
  paste(
    apply(
      sample_info_simple,
      1,
      function(x) {
        paste0(
          x[["score_label"]],
          ": N=",
          x[["N"]],
          ", cases/noncases=",
          x[["N_case"]],
          "/",
          x[["N_noncase"]]
        )
      }
    ),
    collapse = "; "
  )
}

# -----------------------------------------------------------------------------
# Summary statistics
# -----------------------------------------------------------------------------

summary_by_curve <- plot_dat %>%
  group_by(score_type, score_label, metric, metric_label) %>%
  summarise(
    baseline_c_index = c_index_plot[cumulative_step == 0][1],
    final_step = max(cumulative_step, na.rm = TRUE),
    final_c_index = c_index_plot[which.max(cumulative_step)],
    final_delta_vs_base = final_c_index - baseline_c_index,
    best_step = cumulative_step[which.max(c_index_plot)],
    best_c_index = max(c_index_plot, na.rm = TRUE),
    best_delta_vs_base = best_c_index - baseline_c_index,
    .groups = "drop"
  )

comparison_summary <- summary_by_curve %>%
  select(
    score_type,
    metric,
    metric_label,
    final_c_index,
    final_delta_vs_base,
    best_c_index,
    best_delta_vs_base
  ) %>%
  {
    epoch <- filter(., score_type == "EPOCH") %>%
      rename(
        final_c_index_EPOCH = final_c_index,
        final_delta_vs_base_EPOCH = final_delta_vs_base,
        best_c_index_EPOCH = best_c_index,
        best_delta_vs_base_EPOCH = best_delta_vs_base
      ) %>%
      select(-score_type)

    bag <- filter(., score_type == "BAG") %>%
      rename(
        final_c_index_BAG = final_c_index,
        final_delta_vs_base_BAG = final_delta_vs_base,
        best_c_index_BAG = best_c_index,
        best_delta_vs_base_BAG = best_delta_vs_base
      ) %>%
      select(-score_type)

    inner_join(epoch, bag, by = c("metric", "metric_label"))
  } %>%
  mutate(
    final_c_index_difference_EPOCH_minus_BAG =
      final_c_index_EPOCH - final_c_index_BAG,
    final_delta_difference_EPOCH_minus_BAG =
      final_delta_vs_base_EPOCH - final_delta_vs_base_BAG,
    best_c_index_difference_EPOCH_minus_BAG =
      best_c_index_EPOCH - best_c_index_BAG,
    best_delta_difference_EPOCH_minus_BAG =
      best_delta_vs_base_EPOCH - best_delta_vs_base_BAG
  )

# -----------------------------------------------------------------------------
# Plot setup
# -----------------------------------------------------------------------------

max_step <- max(plot_dat$cumulative_step, na.rm = TRUE)

baseline_df <- plot_dat %>%
  group_by(metric_label) %>%
  summarise(
    baseline_c = mean(c_index_plot[cumulative_step == 0], na.rm = TRUE),
    .groups = "drop"
  )

final_label_df <- plot_dat %>%
  group_by(score_label, metric_label) %>%
  filter(cumulative_step == max(cumulative_step, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    final_label = paste0(
      as.character(score_label),
      "\nC = ",
      sprintf("%.3f", c_index_plot),
      "\nΔ = ",
      sprintf("%+.4f", delta_vs_base_plot)
    )
  )

c_min <- min(plot_dat$c_index_plot, na.rm = TRUE)
c_max <- max(plot_dat$c_index_plot, na.rm = TRUE)
c_range <- max(c_max - c_min, 0.01)

sig_label_df <- plot_dat %>%
  filter(
    SHOW_STEPWISE_SIGNIFICANCE,
    cumulative_step > 0,
    stepwise_sig != ""
  ) %>%
  mutate(
    sig_y = case_when(
      score_label == "Mortality EPOCH" ~ c_index_plot + 0.12 * c_range,
      score_label == "Matched BAG" ~ c_index_plot - 0.12 * c_range,
      TRUE ~ c_index_plot + 0.12 * c_range
    )
  )

y_lower <- max(
  0,
  min(
    c_min,
    ifelse(nrow(sig_label_df) > 0L, min(sig_label_df$sig_y, na.rm = TRUE), c_min)
  ) - 0.20 * c_range
)

y_upper <- min(
  1,
  max(
    c_max,
    ifelse(nrow(sig_label_df) > 0L, max(sig_label_df$sig_y, na.rm = TRUE), c_max)
  ) + 0.35 * c_range
)

score_palette <- c(
  "Mortality EPOCH" = "#2F5D8C",
  "Matched BAG" = "#C98245"
)

score_linetypes <- c(
  "Mortality EPOCH" = "solid",
  "Matched BAG" = "dotted"
)

score_shapes <- c(
  "Mortality EPOCH" = 16,
  "Matched BAG" = 17
)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

p <- ggplot(
  plot_dat,
  aes(
    x = cumulative_step,
    y = c_index_plot,
    color = score_label,
    shape = score_label,
    linetype = score_label,
    group = score_label
  )
) +
  geom_hline(
    data = baseline_df,
    aes(yintercept = baseline_c),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.65,
    color = "grey35"
  ) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3.5, stroke = 0.8) +
  scale_color_manual(
    values = score_palette,
    name = ""
  ) +
  scale_shape_manual(
    values = score_shapes,
    name = ""
  ) +
  scale_linetype_manual(
    values = score_linetypes,
    name = ""
  ) +
  scale_x_continuous(
    breaks = seq(0, max_step, by = 1),
    labels = function(x) {
      out <- epoch_x_labels[as.character(x)]
      out[is.na(out)] <- as.character(x[is.na(out)])
      out
    },
    limits = c(-0.2, max_step + ifelse(SHOW_FINAL_LABELS, 2.6, 0.4)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(y_lower, y_upper),
    breaks = scales::pretty_breaks(n = 6),
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Cumulative model step, labeled by EPOCH mortality-clock order",
    y = "Harrell C-index",
    title = paste0(
      DISEASE_NAME,
      " (",
      disease_code_clean,
      "): cumulative prediction by EPOCH versus BAG"
    ),
    subtitle = paste0(
      "Clinical baseline = age, sex, smoking, BMI, diastolic BP, and systolic BP; ",
      "x-axis labels follow the cumulative order of mortality EPOCH clocks"
    ),
    caption = paste0(
      sample_caption,
      " Dashed horizontal line marks the clinical baseline. ",
      "Solid line = mortality EPOCH; dotted line = matched BAG. ",
      "Step-wise significance symbols are based on sequential likelihood-ratio tests: ",
      "* P<0.05, ** P<0.01, *** P<0.001."
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      size = 13.5,
      color = "black",
      face = "plain",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 15,
      color = "black",
      face = "plain"
    ),
    axis.title.x = element_text(
      size = 17,
      face = "bold",
      color = "black",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 17,
      face = "bold",
      color = "black",
      margin = margin(r = 10)
    ),
    axis.line = element_line(
      color = "black",
      linewidth = 0.65
    ),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.6
    ),
    plot.title = element_text(
      size = 17,
      face = "bold",
      color = "black"
    ),
    plot.subtitle = element_text(
      size = 12,
      color = "black"
    ),
    plot.caption = element_text(
      size = 9.2,
      color = "black",
      hjust = 0,
      margin = margin(t = 10)
    ),
    legend.position = "bottom",
    legend.text = element_text(
      size = 13,
      color = "black"
    ),
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    strip.background = element_blank(),
    plot.margin = margin(12, 24, 12, 12)
  )

if (length(unique(plot_dat$metric_label)) > 1L) {
  p <- p +
    facet_wrap(~ metric_label, nrow = 1, scales = "free_y")
}

if (SHOW_POINT_CINDEX_LABELS) {
  p <- p +
    geom_text(
      aes(label = cindex_label),
      size = 3.6,
      color = "black",
      vjust = -1.2,
      show.legend = FALSE
    )
}

if (SHOW_FINAL_LABELS) {
  p <- p +
    geom_text(
      data = final_label_df,
      aes(
        x = cumulative_step + 0.25,
        y = c_index_plot,
        label = final_label,
        color = score_label
      ),
      hjust = 0,
      vjust = 0.5,
      size = 3.4,
      fontface = "bold",
      show.legend = FALSE
    )
}

if (SHOW_STEPWISE_SIGNIFICANCE && nrow(sig_label_df) > 0L) {
  p <- p +
    geom_text(
      data = sig_label_df,
      aes(
        x = cumulative_step,
        y = sig_y,
        label = stepwise_sig,
        color = score_label
      ),
      size = 4.0,
      fontface = "bold",
      show.legend = FALSE
    )
}

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------

safe_disease <- safe_filename(paste0(disease_code_clean, "_", DISEASE_NAME))
metric_tag <- safe_filename(PLOT_METRIC)

pdf_path <- file.path(
  OUTPUT_DIR,
  paste0(
    "HeadToHead_Cumulative_Cindex_EPOCH_vs_BAG_PM_",
    safe_disease,
    "_",
    metric_tag,
    ".pdf"
  )
)

png_path <- file.path(
  OUTPUT_DIR,
  paste0(
    "HeadToHead_Cumulative_Cindex_EPOCH_vs_BAG_PM_",
    safe_disease,
    "_",
    metric_tag,
    ".png"
  )
)

plot_data_path <- file.path(
  OUTPUT_DIR,
  paste0(
    "HeadToHead_Cumulative_Cindex_EPOCH_vs_BAG_PM_",
    safe_disease,
    "_",
    metric_tag,
    "_plot_data.tsv"
  )
)

summary_path <- file.path(
  OUTPUT_DIR,
  paste0(
    "HeadToHead_Cumulative_Cindex_EPOCH_vs_BAG_PM_",
    safe_disease,
    "_",
    metric_tag,
    "_summary.tsv"
  )
)

comparison_path <- file.path(
  OUTPUT_DIR,
  paste0(
    "HeadToHead_Cumulative_Cindex_EPOCH_vs_BAG_PM_",
    safe_disease,
    "_",
    metric_tag,
    "_comparison.tsv"
  )
)

ggsave(
  filename = pdf_path,
  plot = p,
  width = PDF_WIDTH,
  height = PDF_HEIGHT,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)

ggsave(
  filename = png_path,
  plot = p,
  width = PDF_WIDTH,
  height = PDF_HEIGHT,
  units = "in",
  dpi = PNG_DPI,
  bg = "white",
  limitsize = FALSE
)

plot_export <- plot_dat %>%
  transmute(
    disease_id = disease_code_clean,
    disease_name = DISEASE_NAME,
    score_type,
    score_label,
    metric,
    metric_label,
    cumulative_step,
    epoch_x_axis_label = epoch_x_labels[as.character(cumulative_step)],
    added_clock,
    n_clocks,
    N,
    N_case,
    N_noncase,
    event_rate,
    median_followup_years,
    c_index = c_index_plot,
    base_c_index = base_c_index_plot,
    delta_c_index_vs_base = delta_vs_base_plot,
    delta_c_index_vs_previous = delta_vs_previous_plot,
    sequential_lr_p_vs_previous = sequential_p,
    stepwise_significance = stepwise_sig,
    status
  )

readr::write_tsv(plot_export, plot_data_path)
readr::write_tsv(summary_by_curve, summary_path)
readr::write_tsv(comparison_summary, comparison_path)

message("Disease: ", DISEASE_NAME, " (", disease_code_clean, ")")
message("EPOCH input: ", EPOCH_INPUT_FILE)
message("BAG input: ", BAG_INPUT_FILE)
message("Plot metric: ", PLOT_METRIC)
message("Wrote PDF: ", pdf_path)
message("Wrote PNG: ", png_path)
message("Wrote plot data: ", plot_data_path)
message("Wrote summary: ", summary_path)
message("Wrote comparison: ", comparison_path)

print(summary_by_curve)
print(comparison_summary)