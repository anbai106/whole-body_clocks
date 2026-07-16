#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
  library(patchwork)
  library(scales)
})

# ============================================================
# NHANES Model 2 non-disease-input mortality EPOCH performance
#
# Input files from Step II:
#   nhanes_model2_epoch_scores.tsv
#   nhanes_model2_performance.tsv
#   nhanes_model2_alpha_path_validation.tsv
#   nhanes_model2_coefficients.tsv
#
# Main panels:
#   A. Discrimination across temporal splits
#   B. Temporal C-index gap relative to held-out test
#   C. EPOCH acceleration distribution across splits
#   D. Held-out test cumulative mortality by EPOCH quartile
#   E. Observed 3-year mortality by EPOCH quartile
#   F. 3-year calibration using full-model risk score
# ============================================================

# -----------------------------
# Arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) return(default)
  args[idx + 1]
}

INDIR <- get_arg(
  "--indir",
  "/Users/hao/Dropbox/NHANES/output_dir/model2_nondisease_mortality_epoch"
)

OUTDIR <- get_arg(
  "--outdir",
  file.path(INDIR, "figures")
)

HORIZON_YEARS <- as.numeric(get_arg("--horizon_years", "3"))
BASE_FAMILY <- get_arg("--base_family", "Times New Roman")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("Input directory:  ", INDIR)
message("Output directory: ", OUTDIR)
message("Horizon years:    ", HORIZON_YEARS)

# -----------------------------
# EPOCH-style color palette
# -----------------------------
split_cols <- c(
  train = "#3E95B5",
  validation = "#F39C12",
  test = "#78A65A",
  all = "#4D4D4D"
)

quartile_cols <- c(
  "Q1 lowest" = "#66C2A5",
  "Q2" = "#8DA0CB",
  "Q3" = "#E5C494",
  "Q4 highest" = "#FC8D62"
)

score_cols <- c(
  "Full model LP" = "#3E95B5",
  "Feature-only EPOCH" = "#F39C12",
  "EPOCH acceleration" = "#78A65A"
)

coef_cols <- c(
  Positive = "#F39C12",
  Negative = "#3E95B5"
)

theme_epoch <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.35),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.35),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, color = "#2C3E50"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey80"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = base_size + 3)
    )
}

save_plot_pdf_png <- function(plot, basename, width, height) {
  pdf_file <- file.path(OUTDIR, paste0(basename, ".pdf"))
  png_file <- file.path(OUTDIR, paste0(basename, ".png"))
  
  pdf_device <- if (capabilities("cairo")) cairo_pdf else pdf
  
  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    device = pdf_device
  )
  
  ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  
  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
}

# -----------------------------
# Input files
# -----------------------------
score_file <- file.path(INDIR, "nhanes_model2_epoch_scores.tsv")
perf_file <- file.path(INDIR, "nhanes_model2_performance.tsv")
coef_file <- file.path(INDIR, "nhanes_model2_coefficients.tsv")

if (!file.exists(score_file)) stop("Missing score file: ", score_file)
if (!file.exists(perf_file)) stop("Missing performance file: ", perf_file)

scores <- fread(score_file)
perf <- fread(perf_file)

split_levels <- c("train", "validation", "test", "all")
split_labels <- c(
  train = "Train\n1999-2010",
  validation = "Validation\n2011-2014",
  test = "Test\n2015-2018",
  all = "All"
)

score_labels <- c(
  cindex_lp_total = "Full model LP",
  cindex_lp_feature_only = "Feature-only EPOCH",
  cindex_acceleration_z = "EPOCH acceleration"
)

scores[, split := factor(as.character(split), levels = split_levels)]
perf[, split := factor(as.character(split), levels = split_levels)]

num_cols <- c(
  "death",
  "followup_years_exm",
  "RIDAGEYR",
  "mortality_epoch_lp_total",
  "mortality_epoch_lp_feature_only",
  "mortality_epoch_year_equivalent",
  "mortality_epoch_acceleration_years",
  "mortality_epoch_acceleration_z"
)

for (v in intersect(num_cols, names(scores))) {
  scores[, (v) := suppressWarnings(as.numeric(get(v)))]
}

scores <- scores[
  !is.na(death) &
    !is.na(followup_years_exm) &
    followup_years_exm > 0 &
    !is.na(mortality_epoch_acceleration_z) &
    !is.na(mortality_epoch_lp_total)
]

# ============================================================
# Helper functions
# ============================================================
tidy_survfit <- function(fit, strata_prefix = "") {
  s <- summary(fit)
  
  dt <- data.table(
    time = s$time,
    surv = s$surv,
    n_risk = s$n.risk,
    n_event = s$n.event,
    strata = if (is.null(s$strata)) "All" else as.character(s$strata)
  )
  
  if (strata_prefix != "") {
    dt[, strata := sub(strata_prefix, "", strata)]
  }
  
  dt[, cum_mortality := 1 - surv]
  dt
}

km_risk_at_time <- function(dt, horizon_years = 3) {
  dt <- dt[
    !is.na(followup_years_exm) &
      followup_years_exm > 0 &
      !is.na(death)
  ]
  
  if (nrow(dt) < 20) {
    return(data.table(
      risk = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      n = nrow(dt),
      deaths = sum(dt$death == 1, na.rm = TRUE)
    ))
  }
  
  fit <- tryCatch({
    survfit(Surv(followup_years_exm, death) ~ 1, data = dt)
  }, error = function(e) NULL)
  
  if (is.null(fit)) {
    return(data.table(
      risk = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      n = nrow(dt),
      deaths = sum(dt$death == 1, na.rm = TRUE)
    ))
  }
  
  s <- tryCatch({
    summary(fit, times = horizon_years, extend = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(s) || length(s$surv) == 0) {
    return(data.table(
      risk = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      n = nrow(dt),
      deaths = sum(dt$death == 1, na.rm = TRUE)
    ))
  }
  
  surv <- s$surv[1]
  
  if (!is.null(s$lower) && length(s$lower) > 0 && !is.na(s$lower[1])) {
    lower_surv <- s$lower[1]
    upper_surv <- s$upper[1]
    risk_lower <- 1 - upper_surv
    risk_upper <- 1 - lower_surv
  } else {
    se <- ifelse(length(s$std.err) > 0, s$std.err[1], NA_real_)
    if (is.na(se)) {
      risk_lower <- NA_real_
      risk_upper <- NA_real_
    } else {
      lower_surv <- max(0, surv - 1.96 * se)
      upper_surv <- min(1, surv + 1.96 * se)
      risk_lower <- 1 - upper_surv
      risk_upper <- 1 - lower_surv
    }
  }
  
  data.table(
    risk = 1 - surv,
    lower = risk_lower,
    upper = risk_upper,
    n = nrow(dt),
    deaths = sum(dt$death == 1, na.rm = TRUE)
  )
}

cox_hr_one <- function(dt, split_name, score_var) {
  dt2 <- dt[
    !is.na(get(score_var)) &
      !is.na(followup_years_exm) &
      followup_years_exm > 0 &
      !is.na(death)
  ]
  
  if (nrow(dt2) < 50 || sum(dt2$death == 1, na.rm = TRUE) < 10) {
    return(data.table(
      split = split_name,
      n = nrow(dt2),
      deaths = sum(dt2$death == 1, na.rm = TRUE),
      hr = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      p = NA_real_
    ))
  }
  
  fit <- coxph(
    as.formula(paste0("Surv(followup_years_exm, death) ~ ", score_var)),
    data = dt2
  )
  
  sm <- summary(fit)
  
  data.table(
    split = split_name,
    n = nrow(dt2),
    deaths = sum(dt2$death == 1, na.rm = TRUE),
    hr = sm$conf.int[1, "exp(coef)"],
    lower = sm$conf.int[1, "lower .95"],
    upper = sm$conf.int[1, "upper .95"],
    p = sm$coefficients[1, "Pr(>|z|)"]
  )
}

get_cumhaz_at_horizon <- function(cox_fit, horizon_years) {
  bh <- survival::basehaz(cox_fit, centered = FALSE)
  bh <- as.data.table(bh)
  
  if (!"time" %in% names(bh)) {
    stop("basehaz() output does not contain a 'time' column. Columns found: ",
         paste(names(bh), collapse = ", "))
  }
  
  hazard_candidates <- c("hazard", "cumhaz", "cum_hazard", "basehaz")
  hazard_col <- hazard_candidates[hazard_candidates %in% names(bh)][1]
  
  if (is.na(hazard_col)) {
    stop("basehaz() output does not contain a cumulative hazard column. Columns found: ",
         paste(names(bh), collapse = ", "))
  }
  
  bh <- bh[
    is.finite(time) &
      is.finite(get(hazard_col))
  ]
  
  if (nrow(bh) == 0) {
    stop("basehaz() output has no finite time/hazard rows.")
  }
  
  setorder(bh, time)
  
  idx <- which(bh$time <= horizon_years)
  if (length(idx) == 0) {
    use_idx <- 1
  } else {
    use_idx <- max(idx)
  }
  
  as.numeric(bh[[hazard_col]][use_idx])
}

safe_exp <- function(x) {
  exp(pmin(pmax(x, -50), 50))
}

# ============================================================
# Panel A: Discrimination across temporal splits
# FIXED: use geom_rect() with explicit baseline so bars render correctly
# ============================================================
perf_long <- melt(
  perf,
  id.vars = c("split", "n", "deaths", "median_followup_years"),
  measure.vars = intersect(names(score_labels), names(perf)),
  variable.name = "score_type",
  value.name = "cindex"
)

perf_long[, score_type := factor(
  as.character(score_type),
  levels = names(score_labels),
  labels = score_labels
)]

perf_long[, split_label := factor(
  split_labels[as.character(split)],
  levels = split_labels[split_levels]
)]

pA_dt <- copy(perf_long[split %in% c("train", "validation", "test")])
pA_dt[, x_id := as.numeric(split_label)]

cindex_base <- 0.55
cindex_top <- 0.90
bar_width <- 0.62

pA <- ggplot(pA_dt) +
  geom_rect(
    aes(
      xmin = x_id - bar_width / 2,
      xmax = x_id + bar_width / 2,
      ymin = cindex_base,
      ymax = cindex,
      fill = split
    ),
    color = NA
  ) +
  geom_text(
    aes(
      x = x_id,
      y = pmin(cindex + 0.012, cindex_top - 0.006),
      label = sprintf("%.3f", cindex)
    ),
    size = 3.0,
    fontface = "bold"
  ) +
  facet_wrap(~ score_type, nrow = 1) +
  scale_fill_manual(values = split_cols, breaks = c("train", "validation", "test")) +
  scale_x_continuous(
    breaks = seq_along(split_labels[c("train", "validation", "test")]),
    labels = split_labels[c("train", "validation", "test")],
    limits = c(0.45, 3.55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(cindex_base, cindex_top),
    breaks = seq(cindex_base, cindex_top, 0.05),
    expand = expansion(mult = c(0, 0.04))
  ) +
  labs(
    title = "Discrimination across temporal splits",
    subtitle = "C-index in train, validation, and held-out test cycles",
    x = NULL,
    y = "Harrell C-index",
    fill = NULL
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Panel B: Temporal C-index gap relative to held-out test
# ============================================================
perf_gap <- dcast(
  perf_long[split %in% c("train", "validation", "test")],
  score_type ~ split,
  value.var = "cindex"
)

gap_dt <- rbindlist(list(
  data.table(
    score_type = perf_gap$score_type,
    comparison = "Train - test",
    gap = perf_gap$train - perf_gap$test,
    fill_group = "train"
  ),
  data.table(
    score_type = perf_gap$score_type,
    comparison = "Validation - test",
    gap = perf_gap$validation - perf_gap$test,
    fill_group = "validation"
  )
))

gap_dt[, score_type := factor(score_type, levels = score_labels)]
gap_dt[, comparison := factor(comparison, levels = c("Train - test", "Validation - test"))]

pB <- ggplot(gap_dt, aes(x = comparison, y = gap, fill = fill_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_col(width = 0.62, color = NA) +
  geom_text(
    aes(label = sprintf("%+.3f", gap)),
    vjust = ifelse(gap_dt$gap >= 0, -0.45, 1.25),
    size = 3.0,
    fontface = "bold"
  ) +
  facet_wrap(~ score_type, nrow = 1) +
  scale_fill_manual(values = split_cols, breaks = c("train", "validation")) +
  labs(
    title = "Temporal generalization gap",
    subtitle = "Positive values suggest higher apparent performance than held-out test",
    x = NULL,
    y = "C-index difference",
    fill = NULL
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Panel C: EPOCH acceleration score distribution
# ============================================================
dist_dt <- scores[split %in% c("train", "validation", "test")]

xlim_dist <- quantile(
  dist_dt$mortality_epoch_acceleration_z,
  probs = c(0.005, 0.995),
  na.rm = TRUE
)

pC <- ggplot(
  dist_dt,
  aes(x = mortality_epoch_acceleration_z, color = split, fill = split)
) +
  geom_density(alpha = 0.22, linewidth = 0.75, adjust = 1.1) +
  scale_color_manual(values = split_cols, breaks = c("train", "validation", "test")) +
  scale_fill_manual(values = split_cols, breaks = c("train", "validation", "test")) +
  coord_cartesian(xlim = xlim_dist) +
  labs(
    title = "EPOCH acceleration distribution",
    subtitle = "Train, validation, and test distributions should broadly overlap",
    x = "NHANES mortality EPOCH acceleration, z score",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Define train-based quartiles and deciles
# ============================================================
train_z <- scores[split == "train", mortality_epoch_acceleration_z]

qcuts <- quantile(train_z, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
scores[, epoch_q_train := cut(
  mortality_epoch_acceleration_z,
  breaks = c(-Inf, qcuts, Inf),
  labels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"),
  include.lowest = TRUE
)]

dcuts <- quantile(train_z, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE)
scores[, epoch_decile_train := cut(
  mortality_epoch_acceleration_z,
  breaks = c(-Inf, dcuts, Inf),
  labels = paste0("D", 1:10),
  include.lowest = TRUE
)]

# ============================================================
# Panel D: Held-out test cumulative mortality by quartile
# ============================================================
km_test <- scores[split == "test" & !is.na(epoch_q_train)]

fit_test <- survfit(
  Surv(followup_years_exm, death) ~ epoch_q_train,
  data = km_test
)

km_dt <- tidy_survfit(fit_test, strata_prefix = "^epoch_q_train=")
km_dt[, strata := factor(strata, levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))]

pD <- ggplot(km_dt, aes(x = time, y = cum_mortality, color = strata)) +
  geom_step(linewidth = 0.95) +
  scale_color_manual(values = quartile_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0, min(6, max(km_dt$time, na.rm = TRUE)))) +
  labs(
    title = "Held-out test mortality separation",
    subtitle = "Kaplan-Meier cumulative mortality by training-defined EPOCH quartile",
    x = "Years after baseline MEC exam",
    y = "Cumulative mortality",
    color = "EPOCH quartile"
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Panel E: Observed mortality risk by quartile and split
# ============================================================
risk_quartile <- rbindlist(lapply(c("train", "validation", "test"), function(sp) {
  rbindlist(lapply(c("Q1 lowest", "Q2", "Q3", "Q4 highest"), function(q) {
    dt_sub <- scores[split == sp & epoch_q_train == q]
    out <- km_risk_at_time(dt_sub, horizon_years = HORIZON_YEARS)
    out[, split := sp]
    out[, quartile := q]
    out
  }))
}))

risk_quartile[, split := factor(split, levels = c("train", "validation", "test"))]
risk_quartile[, quartile := factor(quartile, levels = c("Q1 lowest", "Q2", "Q3", "Q4 highest"))]

pE <- ggplot(
  risk_quartile,
  aes(x = quartile, y = risk, color = split, group = split)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.16, linewidth = 0.5) +
  scale_color_manual(values = split_cols, breaks = c("train", "validation", "test")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = paste0("Observed ", HORIZON_YEARS, "-year mortality by EPOCH quartile"),
    subtitle = "Quartiles are defined in training and applied to validation/test cycles",
    x = NULL,
    y = paste0("Observed ", HORIZON_YEARS, "-year mortality risk"),
    color = NULL
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Panel F: Calibration at horizon using full model LP
# ============================================================
train_scores <- scores[
  split == "train" &
    !is.na(followup_years_exm) &
    followup_years_exm > 0 &
    !is.na(death) &
    is.finite(mortality_epoch_lp_total)
]

cal_fit <- coxph(
  Surv(followup_years_exm, death) ~ mortality_epoch_lp_total,
  data = train_scores
)

base_haz_horizon <- get_cumhaz_at_horizon(cal_fit, HORIZON_YEARS)
coef_lp <- as.numeric(coef(cal_fit)[["mortality_epoch_lp_total"]])

message("Calibration model coefficient for LP: ", signif(coef_lp, 5))
message("Baseline cumulative hazard at ", HORIZON_YEARS, " years: ", signif(base_haz_horizon, 5))

scores[, pred_horizon_risk := {
  eta <- coef_lp * mortality_epoch_lp_total
  1 - exp(-base_haz_horizon * safe_exp(eta))
}]

scores[!is.finite(pred_horizon_risk), pred_horizon_risk := NA_real_]
scores[pred_horizon_risk < 0, pred_horizon_risk := 0]
scores[pred_horizon_risk > 1, pred_horizon_risk := 1]

train_pred <- scores[split == "train" & is.finite(pred_horizon_risk), pred_horizon_risk]
pred_cuts <- unique(quantile(train_pred, probs = seq(0.1, 0.9, 0.1), na.rm = TRUE))

scores[, pred_risk_decile := cut(
  pred_horizon_risk,
  breaks = c(-Inf, pred_cuts, Inf),
  labels = paste0("D", seq_len(length(pred_cuts) + 1)),
  include.lowest = TRUE
)]

decile_levels <- levels(scores$pred_risk_decile)

calib_dt <- rbindlist(lapply(c("train", "validation", "test"), function(sp) {
  rbindlist(lapply(decile_levels, function(d) {
    dt_sub <- scores[split == sp & pred_risk_decile == d]
    obs <- km_risk_at_time(dt_sub, horizon_years = HORIZON_YEARS)
    
    data.table(
      split = sp,
      decile = d,
      mean_predicted = mean(dt_sub$pred_horizon_risk, na.rm = TRUE),
      observed = obs$risk,
      lower = obs$lower,
      upper = obs$upper,
      n = obs$n,
      deaths = obs$deaths
    )
  }))
}))

calib_dt[, split := factor(split, levels = c("train", "validation", "test"))]
calib_dt <- calib_dt[is.finite(mean_predicted) & is.finite(observed)]

max_cal <- max(
  calib_dt$mean_predicted,
  calib_dt$upper,
  calib_dt$observed,
  na.rm = TRUE
)

if (!is.finite(max_cal) || max_cal <= 0) {
  max_cal <- 0.05
}

pF <- ggplot(
  calib_dt,
  aes(x = mean_predicted, y = observed, color = split, group = split)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 2.4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, linewidth = 0.45, alpha = 0.80) +
  scale_color_manual(values = split_cols, breaks = c("train", "validation", "test")) +
  scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0, max_cal * 1.10), ylim = c(0, max_cal * 1.10)) +
  labs(
    title = paste0("Calibration at ", HORIZON_YEARS, " years"),
    subtitle = "Observed Kaplan-Meier mortality risk versus predicted Cox risk by decile",
    x = paste0("Mean predicted ", HORIZON_YEARS, "-year mortality risk"),
    y = paste0("Observed ", HORIZON_YEARS, "-year mortality risk"),
    color = NULL
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Supplementary panel: HR per 1-SD acceleration
# ============================================================
hr_dt <- rbindlist(lapply(split_levels, function(sp) {
  dt_sp <- if (sp == "all") scores else scores[split == sp]
  cox_hr_one(dt_sp, sp, "mortality_epoch_acceleration_z")
}))

hr_dt[, split_label := factor(
  split_labels[as.character(split)],
  levels = rev(split_labels[split_levels])
)]

hr_dt[, label := sprintf("HR %.2f [%.2f, %.2f]", hr, lower, upper)]

pHR <- ggplot(hr_dt, aes(x = hr, y = split_label, color = split)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 0.70) +
  geom_point(size = 3.2) +
  geom_text(
    aes(label = label),
    hjust = -0.04,
    size = 3.0,
    color = "black"
  ) +
  scale_color_manual(values = split_cols, breaks = split_levels) +
  scale_x_log10() +
  coord_cartesian(xlim = c(
    min(hr_dt$lower, na.rm = TRUE) * 0.90,
    max(hr_dt$upper, na.rm = TRUE) * 1.45
  )) +
  labs(
    title = "Mortality risk per 1-SD EPOCH acceleration",
    subtitle = "Cox model using residualized EPOCH acceleration",
    x = "Hazard ratio, log scale",
    y = NULL,
    color = NULL
  ) +
  theme_epoch(base_size = 11)

# ============================================================
# Supplementary panel: Top nonzero coefficients
# ============================================================
if (file.exists(coef_file)) {
  coef_dt <- fread(coef_file)
  
  if (!"source_variable" %in% names(coef_dt)) {
    coef_dt[, source_variable := transformed_feature]
  }
  
  coef_top <- coef_dt[
    is_epoch_feature == TRUE &
      nonzero == TRUE &
      is.finite(beta) &
      beta != 0
  ]
  
  coef_top[, abs_beta := abs(beta)]
  coef_top[, direction := ifelse(beta >= 0, "Positive", "Negative")]
  coef_top <- coef_top[order(-abs_beta)][1:min(.N, 25)]
  
  coef_top[, feature_label := source_variable]
  coef_top[, feature_label := make.unique(feature_label)]
  coef_top[, feature_label := factor(feature_label, levels = rev(feature_label))]
  
  pCoef <- ggplot(coef_top, aes(x = beta, y = feature_label, fill = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = coef_cols) +
    labs(
      title = "Top nonzero penalized EPOCH features",
      subtitle = "Largest absolute Cox elastic-net coefficients",
      x = "Cox elastic-net coefficient",
      y = NULL,
      fill = "Direction"
    ) +
    theme_epoch(base_size = 10)
} else {
  coef_top <- data.table()
  pCoef <- ggplot() +
    annotate("text", x = 0, y = 0, label = "Coefficient file not found", size = 5) +
    theme_void() +
    labs(title = "Top nonzero penalized EPOCH features")
}

# ============================================================
# Assemble main and supplementary figures
# ============================================================
main_fig <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    title = "NHANES non-disease-input mortality EPOCH model performance",
    subtitle = paste0(
      "Temporal train/validation/test evaluation using NHANES 1999-2018 linked mortality data; ",
      "risk horizon = ", HORIZON_YEARS, " years"
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, family = BASE_FAMILY),
      plot.subtitle = element_text(size = 11, family = BASE_FAMILY)
    )
  )

supp_fig <- (pHR | pCoef) +
  plot_annotation(
    title = "Supplementary NHANES mortality EPOCH summaries",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, family = BASE_FAMILY)
    )
  )

save_plot_pdf_png(
  main_fig,
  "NHANES_Model2_mortality_EPOCH_performance_main",
  width = 14,
  height = 15
)

save_plot_pdf_png(
  supp_fig,
  "NHANES_Model2_mortality_EPOCH_performance_supplement",
  width = 14,
  height = 6
)

# ============================================================
# Save plotting tables
# ============================================================
fwrite(perf_long, file.path(OUTDIR, "plot_table_panel_A_cindex.tsv"), sep = "\t")
fwrite(gap_dt, file.path(OUTDIR, "plot_table_panel_B_temporal_gap.tsv"), sep = "\t")
fwrite(km_dt, file.path(OUTDIR, "plot_table_panel_D_test_km_quartiles.tsv"), sep = "\t")
fwrite(risk_quartile, file.path(OUTDIR, "plot_table_panel_E_observed_risk_quartiles.tsv"), sep = "\t")
fwrite(calib_dt, file.path(OUTDIR, "plot_table_panel_F_calibration.tsv"), sep = "\t")
fwrite(hr_dt, file.path(OUTDIR, "plot_table_supp_HR_epoch_acceleration.tsv"), sep = "\t")

if (nrow(coef_top) > 0) {
  fwrite(coef_top, file.path(OUTDIR, "plot_table_supp_top_coefficients.tsv"), sep = "\t")
}

message("\nDone.")
message("Main figure:")
message("  ", file.path(OUTDIR, "NHANES_Model2_mortality_EPOCH_performance_main.pdf"))
message("  ", file.path(OUTDIR, "NHANES_Model2_mortality_EPOCH_performance_main.png"))
message("Supplementary figure:")
message("  ", file.path(OUTDIR, "NHANES_Model2_mortality_EPOCH_performance_supplement.pdf"))
message("  ", file.path(OUTDIR, "NHANES_Model2_mortality_EPOCH_performance_supplement.png"))