# ============================================================
# Brain MRI mortality clock performance visualization
# RStudio-ready direct-run script
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(jsonlite)
  library(patchwork)
  library(scales)
  library(glue)
})

# ============================================================
# 1. Set paths and plotting options
# ============================================================

clock_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/heart_mri_mortality_clock"

prediction_file <- file.path(clock_dir, "heart_mri_mortality_clock_predictions.tsv")
performance_file <- file.path(clock_dir, "heart_mri_mortality_clock_performance.json")

# Choose the calibration horizon.
# Use 5 if risk_5y exists and is not empty.
# If your follow-up is short, you can change this to 3.
horizon_years <- 5

# Bootstrap iterations for C-index confidence intervals.
# Use 200 for quick testing, 500 or 1000 for final figures.
n_boot <- 500

# Output files
out_pdf <- file.path(clock_dir, "heart_mri_mortality_clock_model_performance.pdf")
out_png <- file.path(clock_dir, "heart_mri_mortality_clock_model_performance.png")
out_overfit <- file.path(clock_dir, "heart_mri_mortality_clock_overfitting_summary.tsv")
out_quartile <- file.path(clock_dir, "heart_mri_mortality_clock_risk_quartile_summary.tsv")

stopifnot(file.exists(prediction_file))
stopifnot(file.exists(performance_file))

message("Reading: ", prediction_file)
message("Reading: ", performance_file)

# ============================================================
# 2. Load data
# ============================================================

pred <- readr::read_tsv(
  prediction_file,
  show_col_types = FALSE,
  progress = FALSE
)

perf <- jsonlite::fromJSON(performance_file)

risk_score_col <- "heart_mri_mortality_risk_score"
accel_z_col <- "heart_mri_mortality_clock_acceleration_z"
accel_year_col <- "heart_mri_mortality_clock_acceleration_years"
clock_age_col <- "heart_mri_mortality_clock_age_years"

if (!accel_year_col %in% colnames(pred)) {
  pred[[accel_year_col]] <- NA_real_
}
if (!clock_age_col %in% colnames(pred)) {
  pred[[clock_age_col]] <- NA_real_
}

required_cols <- c(
  "participant_id",
  "time_years",
  "event",
  "split",
  "age_at_imaging",
  "sex",
  risk_score_col,
  accel_z_col
)

missing_cols <- setdiff(required_cols, colnames(pred))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

pred <- pred %>%
  mutate(
    event = case_when(
      event %in% c(TRUE, "True", "TRUE", "true", "1", 1) ~ TRUE,
      event %in% c(FALSE, "False", "FALSE", "false", "0", 0) ~ FALSE,
      TRUE ~ as.logical(event)
    ),
    split = factor(split, levels = c("train", "validation", "test")),
    risk_score = .data[[risk_score_col]],
    clock_accel_z = .data[[accel_z_col]],
    clock_accel_years = .data[[accel_year_col]],
    clock_age_years = .data[[clock_age_col]]
  ) %>%
  filter(
    !is.na(time_years),
    !is.na(event),
    !is.na(split),
    !is.na(risk_score),
    !is.na(clock_accel_z)
  )

message("N total: ", nrow(pred))
message("Deaths total: ", sum(pred$event, na.rm = TRUE))
message("Median follow-up: ", round(median(pred$time_years, na.rm = TRUE), 2), " years")

# ============================================================
# 3. Detect usable absolute-risk column
# ============================================================

risk_col_requested <- paste0("risk_", horizon_years, "y")
risk_cols_available <- grep("^risk_[0-9.]+y$", colnames(pred), value = TRUE)

risk_col <- NULL

if (
  risk_col_requested %in% colnames(pred) &&
  any(is.finite(pred[[risk_col_requested]]), na.rm = TRUE)
) {
  risk_col <- risk_col_requested
} else {
  usable_risk_cols <- risk_cols_available[
    purrr::map_lgl(
      risk_cols_available,
      ~ any(is.finite(pred[[.x]]), na.rm = TRUE)
    )
  ]

  if (length(usable_risk_cols) > 0) {
    risk_col <- usable_risk_cols[[1]]
    horizon_years <- as.numeric(stringr::str_match(risk_col, "^risk_([0-9.]+)y$")[, 2])
    message("Requested absolute-risk column unavailable. Using: ", risk_col)
  } else {
    warning("No usable absolute-risk column found. Calibration panel will be skipped.")
  }
}

# ============================================================
# 4. Plot theme
# ============================================================

theme_clock <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        color = "#566573",
        margin = margin(b = 8)
      ),
      plot.caption = element_text(size = base_size - 3, color = "#7B7D7D"),
      axis.title = element_text(face = "bold", size = base_size - 1),
      axis.text = element_text(color = "#2C3E50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#EAECEE", linewidth = 0.35),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

split_cols <- c(
  "train" = "#2E86AB",
  "validation" = "#F18F01",
  "test" = "#6A994E"
)

quartile_cols <- c(
  "Q1 lowest risk" = "#1B9E77",
  "Q2" = "#7570B3",
  "Q3" = "#E6AB02",
  "Q4 highest risk" = "#D95F02"
)

# ============================================================
# 5. Helper functions
# ============================================================

calc_cindex <- function(df) {
  if (nrow(df) < 5 || sum(df$event, na.rm = TRUE) < 2) {
    return(NA_real_)
  }

  # survival::concordance assumes larger predictor means longer survival.
  # Since larger risk_score means higher mortality hazard, use -risk_score.
  fit <- survival::concordance(
    survival::Surv(time_years, event) ~ I(-risk_score),
    data = df
  )

  as.numeric(fit$concordance)
}

boot_cindex <- function(df, B = 500, seed = 2026) {
  set.seed(seed)

  observed <- calc_cindex(df)

  if (nrow(df) < 20 || sum(df$event, na.rm = TRUE) < 5) {
    return(tibble(
      cindex = observed,
      cindex_lower = NA_real_,
      cindex_upper = NA_real_,
      n_boot_ok = 0
    ))
  }

  vals <- replicate(B, {
    idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
    d <- df[idx, , drop = FALSE]

    if (sum(d$event, na.rm = TRUE) < 2) {
      return(NA_real_)
    }

    calc_cindex(d)
  })

  vals <- vals[is.finite(vals)]

  tibble(
    cindex = observed,
    cindex_lower = as.numeric(quantile(vals, 0.025, na.rm = TRUE)),
    cindex_upper = as.numeric(quantile(vals, 0.975, na.rm = TRUE)),
    n_boot_ok = length(vals)
  )
}

km_risk_at_time <- function(df, tau) {
  if (nrow(df) < 5 || sum(df$event, na.rm = TRUE) < 1) {
    return(tibble(
      observed_risk = NA_real_,
      observed_lower = NA_real_,
      observed_upper = NA_real_
    ))
  }

  fit <- survival::survfit(
    survival::Surv(time_years, event) ~ 1,
    data = df,
    conf.type = "log-log"
  )

  ss <- summary(fit, times = tau, extend = TRUE)

  if (length(ss$surv) == 0) {
    return(tibble(
      observed_risk = NA_real_,
      observed_lower = NA_real_,
      observed_upper = NA_real_
    ))
  }

  tibble(
    observed_risk = 1 - ss$surv[[1]],
    observed_lower = 1 - ss$upper[[1]],
    observed_upper = 1 - ss$lower[[1]]
  )
}

tidy_survfit <- function(fit) {
  ss <- summary(fit)

  tibble(
    time = ss$time,
    survival = ss$surv,
    lower = ss$lower,
    upper = ss$upper,
    strata = if (is.null(ss$strata)) "All" else as.character(ss$strata)
  ) %>%
    mutate(
      strata = stringr::str_replace(strata, "^risk_quartile=", "")
    )
}

# ============================================================
# 6. C-index and overfitting
# ============================================================

message("Computing bootstrap C-index confidence intervals...")

cindex_tbl <- pred %>%
  group_by(split) %>%
  group_modify(~ boot_cindex(.x, B = n_boot, seed = 2026)) %>%
  ungroup() %>%
  mutate(split = factor(split, levels = c("train", "validation", "test")))

json_cindex_tbl <- tibble(
  split = factor(c("train", "validation", "test"), levels = c("train", "validation", "test")),
  cindex_json = c(
    ifelse(is.null(perf$cindex_train), NA_real_, perf$cindex_train),
    ifelse(is.null(perf$cindex_validation), NA_real_, perf$cindex_validation),
    ifelse(is.null(perf$cindex_test), NA_real_, perf$cindex_test)
  )
)

cindex_tbl <- cindex_tbl %>%
  left_join(json_cindex_tbl, by = "split")

test_cindex <- cindex_tbl %>%
  filter(split == "test") %>%
  pull(cindex)

overfit_tbl <- cindex_tbl %>%
  mutate(
    optimism_vs_test = cindex - test_cindex,
    n = purrr::map_int(as.character(split), ~ sum(pred$split == .x)),
    events = purrr::map_int(as.character(split), ~ sum(pred$split == .x & pred$event)),
    overfitting_flag = case_when(
      split == "train" & optimism_vs_test > 0.05 ~ "Possible overfitting",
      split == "train" & optimism_vs_test > 0.02 ~ "Mild optimism",
      split == "train" ~ "Low optimism",
      TRUE ~ ""
    )
  )

readr::write_tsv(overfit_tbl, out_overfit)

# ============================================================
# 7. Risk quartiles using training-set cutoffs
# ============================================================

train_risk <- pred %>%
  filter(split == "train") %>%
  pull(risk_score)

risk_breaks <- quantile(
  train_risk,
  probs = seq(0, 1, by = 0.25),
  na.rm = TRUE,
  type = 8
)

risk_breaks[1] <- -Inf
risk_breaks[length(risk_breaks)] <- Inf

if (length(unique(risk_breaks)) < length(risk_breaks)) {
  warning("Non-unique training risk quartile breaks. Using within-split quartiles.")
  pred <- pred %>%
    group_by(split) %>%
    mutate(
      risk_quartile = ntile(risk_score, 4),
      risk_quartile = factor(
        risk_quartile,
        levels = 1:4,
        labels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
      )
    ) %>%
    ungroup()
} else {
  pred <- pred %>%
    mutate(
      risk_quartile = cut(
        risk_score,
        breaks = risk_breaks,
        include.lowest = TRUE,
        labels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
      )
    )
}

# ============================================================
# 8. Observed mortality by risk quartile
# ============================================================

quartile_risk_tbl <- pred %>%
  filter(!is.na(risk_quartile)) %>%
  group_by(split, risk_quartile) %>%
  group_modify(~ {
    d <- .x
    km <- km_risk_at_time(d, horizon_years)

    tibble(
      n = nrow(d),
      events = sum(d$event, na.rm = TRUE),
      mean_risk_score = mean(d$risk_score, na.rm = TRUE),
      mean_clock_accel_z = mean(d$clock_accel_z, na.rm = TRUE)
    ) %>%
      bind_cols(km)
  }) %>%
  ungroup()

readr::write_tsv(quartile_risk_tbl, out_quartile)

# ============================================================
# 9. Calibration table
# ============================================================

if (!is.null(risk_col)) {
  cal_input <- pred %>%
    mutate(predicted_risk = .data[[risk_col]]) %>%
    filter(is.finite(predicted_risk))

  if (nrow(cal_input) > 0) {
    cal_tbl <- cal_input %>%
      group_by(split) %>%
      mutate(cal_bin = ntile(predicted_risk, 10)) %>%
      ungroup() %>%
      group_by(split, cal_bin) %>%
      group_modify(~ {
        d <- .x
        km <- km_risk_at_time(d, horizon_years)

        tibble(
          n = nrow(d),
          events = sum(d$event, na.rm = TRUE),
          predicted_risk = mean(d$predicted_risk, na.rm = TRUE)
        ) %>%
          bind_cols(km)
      }) %>%
      ungroup()

    out_cal <- file.path(
      clock_dir,
      paste0("heart_mri_mortality_clock_calibration_", horizon_years, "y.tsv")
    )
    readr::write_tsv(cal_tbl, out_cal)
  } else {
    cal_tbl <- tibble()
  }
} else {
  cal_tbl <- tibble()
}

# ============================================================
# 10. Test-set Kaplan-Meier by risk quartile
# ============================================================

test_for_km <- pred %>%
  filter(split == "test", !is.na(risk_quartile))

if (sum(test_for_km$event, na.rm = TRUE) >= 2) {
  km_fit <- survival::survfit(
    survival::Surv(time_years, event) ~ risk_quartile,
    data = test_for_km
  )

  km_df <- tidy_survfit(km_fit)
} else {
  km_df <- tibble()
  warning("Too few test-set events for Kaplan-Meier by quartile.")
}

# ============================================================
# 11. Plot panels
# ============================================================

p_cindex <- cindex_tbl %>%
  ggplot(aes(x = split, y = cindex, fill = split)) +
  geom_col(width = 0.62, alpha = 0.94, color = "white", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = cindex_lower, ymax = cindex_upper),
    width = 0.14,
    linewidth = 0.65,
    color = "#1C2833"
  ) +
  geom_text(
    aes(label = sprintf("%.3f", cindex)),
    vjust = -0.55,
    size = 3.8,
    fontface = "bold",
    color = "#1C2833"
  ) +
  scale_fill_manual(values = split_cols) +
  coord_cartesian(
    ylim = c(
      0.45,
      max(cindex_tbl$cindex_upper, cindex_tbl$cindex, na.rm = TRUE) + 0.05
    )
  ) +
  labs(
    title = "A. Discrimination across splits",
    subtitle = glue("C-index with bootstrap 95% CI, B = {n_boot}"),
    x = NULL,
    y = "C-index"
  ) +
  theme_clock()

gap_tbl <- overfit_tbl %>%
  filter(split %in% c("train", "validation")) %>%
  mutate(
    comparison = case_when(
      split == "train" ~ "Train - test",
      split == "validation" ~ "Validation - test",
      TRUE ~ as.character(split)
    )
  )

p_gap <- gap_tbl %>%
  ggplot(aes(x = comparison, y = optimism_vs_test, fill = split)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#7B7D7D", linewidth = 0.5) +
  geom_col(width = 0.55, alpha = 0.94, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%+.3f", optimism_vs_test)),
    vjust = ifelse(gap_tbl$optimism_vs_test >= 0, -0.55, 1.25),
    size = 3.8,
    fontface = "bold"
  ) +
  scale_fill_manual(values = split_cols) +
  labs(
    title = "B. Optimism gap",
    subtitle = "Large positive train-test gap suggests overfitting",
    x = NULL,
    y = "C-index difference"
  ) +
  theme_clock()

p_density <- pred %>%
  ggplot(aes(x = risk_score, fill = split, color = split)) +
  geom_density(alpha = 0.22, linewidth = 0.75) +
  scale_fill_manual(values = split_cols) +
  scale_color_manual(values = split_cols) +
  labs(
    title = "C. Risk-score distribution",
    subtitle = "Train, validation, and test distributions should broadly overlap",
    x = "Brain MRI mortality risk score",
    y = "Density"
  ) +
  theme_clock()

if (nrow(km_df) > 0) {
  p_km <- km_df %>%
    ggplot(aes(x = time, y = 1 - survival, color = strata, fill = strata)) +
    geom_step(linewidth = 0.9) +
    geom_ribbon(
      aes(ymin = 1 - upper, ymax = 1 - lower),
      alpha = 0.11,
      color = NA
    ) +
    scale_color_manual(values = quartile_cols) +
    scale_fill_manual(values = quartile_cols) +
    scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    labs(
      title = "D. Test-set mortality separation",
      subtitle = "Kaplan-Meier cumulative mortality by training-defined risk quartile",
      x = "Years after imaging",
      y = "Cumulative mortality"
    ) +
    theme_clock()
} else {
  p_km <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Kaplan-Meier panel unavailable: too few test-set events.",
      size = 4.5,
      fontface = "bold"
    ) +
    theme_void() +
    labs(title = "D. Test-set mortality separation")
}

p_quartile <- quartile_risk_tbl %>%
  ggplot(aes(x = risk_quartile, y = observed_risk, color = split, group = split)) +
  geom_pointrange(
    aes(ymin = observed_lower, ymax = observed_upper),
    position = position_dodge(width = 0.45),
    linewidth = 0.55,
    size = 0.45
  ) +
  geom_line(
    position = position_dodge(width = 0.45),
    linewidth = 0.65,
    alpha = 0.85
  ) +
  scale_color_manual(values = split_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title = glue("E. Observed {horizon_years}-year mortality by risk quartile"),
    subtitle = "Monotonic validation/test patterns support generalization",
    x = NULL,
    y = glue("Observed {horizon_years}-year mortality")
  ) +
  theme_clock() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

if (nrow(cal_tbl) > 0) {
  max_axis <- max(
    cal_tbl$predicted_risk,
    cal_tbl$observed_upper,
    cal_tbl$observed_risk,
    na.rm = TRUE
  )

  p_cal <- cal_tbl %>%
    ggplot(aes(x = predicted_risk, y = observed_risk, color = split)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "#7B7D7D",
      linewidth = 0.55
    ) +
    geom_errorbar(
      aes(ymin = observed_lower, ymax = observed_upper),
      width = 0,
      alpha = 0.65,
      linewidth = 0.45
    ) +
    geom_point(size = 2.2, alpha = 0.9) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.75, alpha = 0.8) +
    scale_color_manual(values = split_cols) +
    scale_x_continuous(
      labels = percent_format(accuracy = 0.1),
      limits = c(0, max_axis * 1.08)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.1),
      limits = c(0, max_axis * 1.08)
    ) +
    labs(
      title = glue("F. Calibration at {horizon_years} years"),
      subtitle = "Observed Kaplan-Meier risk versus predicted Cox risk by decile",
      x = glue("Mean predicted {horizon_years}-year mortality"),
      y = glue("Observed {horizon_years}-year mortality")
    ) +
    theme_clock()
} else {
  p_cal <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Calibration skipped: no usable absolute-risk column.",
      size = 4.5,
      fontface = "bold"
    ) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(
      title = "F. Calibration unavailable",
      x = NULL,
      y = NULL
    ) +
    theme_void()
}

# ============================================================
# 12. Final multi-panel figure
# ============================================================

admin_censor_date <- ifelse(
  is.null(perf$admin_censor_date),
  "not recorded",
  perf$admin_censor_date
)

subtitle_text <- glue(
  "N = {nrow(pred)}, deaths = {sum(pred$event)}, ",
  "median follow-up = {round(median(pred$time_years, na.rm = TRUE), 2)} years; ",
  "administrative censor date = {admin_censor_date}"
)

final_fig <-
  (p_cindex | p_gap) /
  (p_density | p_km) /
  (p_quartile | p_cal) +
  plot_annotation(
    title = "Brain MRI mortality clock performance",
    subtitle = subtitle_text,
    caption = "Risk quartiles are defined from the training set and applied to validation/test. Higher risk score indicates greater mortality proximity.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 19, color = "#17202A"),
      plot.subtitle = element_text(size = 11, color = "#566573"),
      plot.caption = element_text(size = 9, color = "#7B7D7D")
    )
  )

print(final_fig)

# ============================================================
# 13. Save outputs
# ============================================================

if (capabilities("cairo")) {
  ggsave(
    filename = out_pdf,
    plot = final_fig,
    width = 15,
    height = 17,
    device = cairo_pdf
  )
} else {
  ggsave(
    filename = out_pdf,
    plot = final_fig,
    width = 15,
    height = 17
  )
}

ggsave(
  filename = out_png,
  plot = final_fig,
  width = 15,
  height = 17,
  dpi = 320
)

message("Saved figure:")
message("  ", out_pdf)
message("  ", out_png)

message("Saved summary tables:")
message("  ", out_overfit)
message("  ", out_quartile)

if (exists("out_cal")) {
  message("  ", out_cal)
}

message("\n===== Overfitting summary =====")
print(
  overfit_tbl %>%
    select(
      split,
      n,
      events,
      cindex,
      cindex_lower,
      cindex_upper,
      optimism_vs_test,
      overfitting_flag,
      cindex_json
    )
)

message("\nInterpretation guide:")
message("  - Train C-index much higher than validation/test indicates overfitting.")
message("  - Small train-test optimism suggests acceptable generalization.")
message("  - KM curves should separate monotonically from Q1 to Q4 in the test set.")
message("  - Calibration points should lie near the diagonal; systematic deviation suggests miscalibration.")