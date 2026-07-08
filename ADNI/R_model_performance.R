# ============================================================
# Plot ADNI brain MRI AD L'EPOCH performance and group effects
#
# Main outputs:
#   1. Model C-index performance plot
#   2. Test-set delta C-index plot
#   3. Test-set KM cumulative incidence by training-defined risk quartile
#   4. AD L'EPOCH acceleration-years by conversion group
#      with raw pairwise Welch t-test P-values and sample sizes
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(survival)
  library(scales)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

pred_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/adni_brain_mri_ad_lepoch_predictions.tsv"

out_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Read predictions
# ------------------------------------------------------------

pred <- readr::read_tsv(pred_file, show_col_types = FALSE)

required_cols <- c(
  "PTID",
  "split",
  "time_years",
  "event",
  "event_or_censor_dx",
  "adni_brain_mri_ad_lepoch_risk_score",
  "adni_brain_mri_ad_lepoch_acceleration_years"
)

missing_cols <- setdiff(required_cols, colnames(pred))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 3. Clean event and outcome group variables
# ------------------------------------------------------------

pred2 <- pred %>%
  mutate(
    event = case_when(
      event %in% c(TRUE, "TRUE", "True", "true", 1, "1") ~ TRUE,
      event %in% c(FALSE, "FALSE", "False", "false", 0, "0") ~ FALSE,
      TRUE ~ NA
    ),
    event_or_censor_dx = as.character(event_or_censor_dx),
    event_or_censor_dx = toupper(event_or_censor_dx),
    
    conversion_group = case_when(
      event == TRUE & event_or_censor_dx == "MCI" ~ "CN to MCI",
      event == TRUE & event_or_censor_dx == "AD"  ~ "CN to AD",
      event == FALSE ~ "Censored / non-event",
      TRUE ~ "Other / unclear"
    ),
    conversion_group = factor(
      conversion_group,
      levels = c("Censored / non-event", "CN to MCI", "CN to AD", "Other / unclear")
    ),
    
    split = factor(
      split,
      levels = c("train", "validation", "test")
    )
  ) %>%
  filter(!is.na(time_years), !is.na(event))

plot_group_df <- pred2 %>%
  filter(conversion_group %in% c("Censored / non-event", "CN to MCI", "CN to AD")) %>%
  mutate(
    conversion_group = factor(
      conversion_group,
      levels = c("Censored / non-event", "CN to MCI", "CN to AD")
    )
  ) %>%
  droplevels()

# ------------------------------------------------------------
# 4. Define model risk-score columns
# ------------------------------------------------------------

model_cols <- c(
  "risk_score_M0_age_sex",
  "risk_score_M1_covariate_baseline",
  "risk_score_M2_brain_mri_only",
  "risk_score_M3_full_covariates_plus_brain_mri"
)

model_cols <- model_cols[model_cols %in% colnames(pred2)]

if (length(model_cols) == 0) {
  stop("No risk_score_M* columns found in predictions file.")
}

model_labels <- c(
  "risk_score_M0_age_sex" = "M0: age + sex",
  "risk_score_M1_covariate_baseline" = "M1: covariates",
  "risk_score_M2_brain_mri_only" = "M2: brain MRI only",
  "risk_score_M3_full_covariates_plus_brain_mri" = "M3: covariates + brain MRI"
)

primary_m3_risk_col <- if ("risk_score_M3_full_covariates_plus_brain_mri" %in% colnames(pred2)) {
  "risk_score_M3_full_covariates_plus_brain_mri"
} else {
  "adni_brain_mri_ad_lepoch_risk_score"
}

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

format_p <- function(p) {
  case_when(
    is.na(p) ~ "P = NA",
    p < 2.2e-16 ~ "P < 2.2e-16",
    p < 0.001 ~ paste0("P = ", formatC(p, format = "e", digits = 2)),
    TRUE ~ paste0("P = ", signif(p, 3))
  )
}

calc_cindex <- function(df, risk_col) {
  d <- df %>%
    select(time_years, event, all_of(risk_col)) %>%
    rename(risk = all_of(risk_col)) %>%
    mutate(
      time_years = as.numeric(time_years),
      risk = as.numeric(risk)
    ) %>%
    filter(!is.na(time_years), !is.na(event), !is.na(risk))
  
  if (nrow(d) < 10 || sum(d$event, na.rm = TRUE) < 2) {
    return(tibble(
      n = nrow(d),
      n_events = sum(d$event, na.rm = TRUE),
      cindex = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_
    ))
  }
  
  cc <- survival::concordance(
    survival::Surv(time_years, event) ~ risk,
    data = d,
    reverse = TRUE
  )
  
  cval <- as.numeric(cc$concordance)
  se <- sqrt(as.numeric(cc$var))
  
  tibble(
    n = nrow(d),
    n_events = sum(d$event, na.rm = TRUE),
    cindex = cval,
    se = se,
    ci_low = pmax(0, cval - 1.96 * se),
    ci_high = pmin(1, cval + 1.96 * se)
  )
}

# ------------------------------------------------------------
# 6. Compute and plot C-index by model and split
# ------------------------------------------------------------

cindex_tbl <- expand_grid(
  split = levels(droplevels(pred2$split)),
  risk_col = model_cols
) %>%
  filter(!is.na(split)) %>%
  mutate(
    perf = map2(split, risk_col, ~ {
      calc_cindex(pred2 %>% filter(split == .x), .y)
    })
  ) %>%
  unnest(perf) %>%
  mutate(
    model = recode(risk_col, !!!model_labels),
    model = factor(model, levels = unname(model_labels[model_cols])),
    split = factor(split, levels = c("train", "validation", "test"))
  )

readr::write_tsv(
  cindex_tbl,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cindex_by_split.tsv")
)

vg_model_palette <- c(
  "M0: age + sex" = "#8C6D31",
  "M1: covariates" = "#2A6F9E",
  "M2: brain MRI only" = "#DDAA33",
  "M3: covariates + brain MRI" = "#1F3B73"
)

p_cindex <- ggplot(
  cindex_tbl,
  aes(x = split, y = cindex, color = model, group = model)
) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 3.0) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.10,
    linewidth = 0.55,
    alpha = 0.90
  ) +
  scale_color_manual(values = vg_model_palette) +
  scale_y_continuous(
    name = "C-index",
    limits = c(0.45, 1.00),
    breaks = seq(0.5, 1.0, by = 0.1),
    labels = number_format(accuracy = 0.01)
  ) +
  scale_x_discrete(name = NULL) +
  labs(
    title = "ADNI brain MRI AD L’EPOCH model performance",
    subtitle = "C-index by split. Higher Cox risk scores indicate greater proximity to MCI/AD conversion.",
    color = NULL,
    caption = "Error bars show approximate 95% confidence intervals from survival::concordance."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cindex_by_split.pdf"),
  p_cindex,
  width = 8.5,
  height = 5.5
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cindex_by_split.png"),
  p_cindex,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Bootstrap delta C-index on test split
# ------------------------------------------------------------

calc_cindex_single <- function(df, risk_col) {
  out <- calc_cindex(df, risk_col)
  out$cindex[[1]]
}

bootstrap_delta_cindex <- function(df, model_a, model_b, n_boot = 1000, seed = 20260707) {
  set.seed(seed)
  
  d <- df %>%
    select(time_years, event, all_of(model_a), all_of(model_b)) %>%
    rename(
      risk_a = all_of(model_a),
      risk_b = all_of(model_b)
    ) %>%
    mutate(
      time_years = as.numeric(time_years),
      risk_a = as.numeric(risk_a),
      risk_b = as.numeric(risk_b)
    ) %>%
    filter(!is.na(time_years), !is.na(event), !is.na(risk_a), !is.na(risk_b))
  
  if (nrow(d) < 10 || sum(d$event) < 2) {
    return(tibble(
      model_a = model_a,
      model_b = model_b,
      n = nrow(d),
      n_events = sum(d$event),
      cindex_a = NA_real_,
      cindex_b = NA_real_,
      delta_cindex = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_one_sided_delta_le_0 = NA_real_,
      p_two_sided = NA_real_
    ))
  }
  
  c_a <- calc_cindex_single(d, "risk_a")
  c_b <- calc_cindex_single(d, "risk_b")
  delta_obs <- c_a - c_b
  
  boot_delta <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(d)), size = nrow(d), replace = TRUE)
    db <- d[idx, , drop = FALSE]
    
    if (sum(db$event) < 2) {
      return(NA_real_)
    }
    
    ca <- tryCatch(calc_cindex_single(db, "risk_a"), error = function(e) NA_real_)
    cb <- tryCatch(calc_cindex_single(db, "risk_b"), error = function(e) NA_real_)
    
    ca - cb
  })
  
  boot_delta <- boot_delta[is.finite(boot_delta)]
  
  if (length(boot_delta) < 20) {
    ci_low <- NA_real_
    ci_high <- NA_real_
    p_le_0 <- NA_real_
    p_two <- NA_real_
  } else {
    ci_low <- as.numeric(quantile(boot_delta, 0.025, na.rm = TRUE))
    ci_high <- as.numeric(quantile(boot_delta, 0.975, na.rm = TRUE))
    p_le_0 <- mean(boot_delta <= 0, na.rm = TRUE)
    p_ge_0 <- mean(boot_delta >= 0, na.rm = TRUE)
    p_two <- min(1, 2 * min(p_le_0, p_ge_0))
  }
  
  tibble(
    model_a = model_a,
    model_b = model_b,
    n = nrow(d),
    n_events = sum(d$event),
    cindex_a = c_a,
    cindex_b = c_b,
    delta_cindex = delta_obs,
    ci_low = ci_low,
    ci_high = ci_high,
    p_one_sided_delta_le_0 = p_le_0,
    p_two_sided = p_two
  )
}

test_df <- pred2 %>% filter(split == "test")

delta_list <- list()

if (
  all(c(
    "risk_score_M3_full_covariates_plus_brain_mri",
    "risk_score_M1_covariate_baseline"
  ) %in% colnames(test_df))
) {
  delta_list[["M3_vs_M1"]] <- bootstrap_delta_cindex(
    test_df,
    model_a = "risk_score_M3_full_covariates_plus_brain_mri",
    model_b = "risk_score_M1_covariate_baseline",
    n_boot = 1000
  )
}

if (
  all(c(
    "risk_score_M3_full_covariates_plus_brain_mri",
    "risk_score_M2_brain_mri_only"
  ) %in% colnames(test_df))
) {
  delta_list[["M3_vs_M2"]] <- bootstrap_delta_cindex(
    test_df,
    model_a = "risk_score_M3_full_covariates_plus_brain_mri",
    model_b = "risk_score_M2_brain_mri_only",
    n_boot = 1000
  )
}

if (length(delta_list) > 0) {
  delta_tbl <- bind_rows(delta_list, .id = "comparison") %>%
    mutate(
      comparison_label = recode(
        comparison,
        "M3_vs_M1" = "M3 full − M1 covariates",
        "M3_vs_M2" = "M3 full − M2 brain MRI only"
      ),
      comparison_label = factor(
        comparison_label,
        levels = c("M3 full − M1 covariates", "M3 full − M2 brain MRI only")
      )
    )
  
  readr::write_tsv(
    delta_tbl,
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_delta_cindex_bootstrap.tsv")
  )
  
  p_delta <- ggplot(
    delta_tbl,
    aes(x = delta_cindex, y = comparison_label)
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5,
      color = "grey50"
    ) +
    geom_errorbarh(
      aes(xmin = ci_low, xmax = ci_high),
      height = 0.18,
      linewidth = 0.65,
      color = "#1F3B73"
    ) +
    geom_point(
      size = 3.4,
      color = "#1F3B73"
    ) +
    scale_x_continuous(
      name = "ΔC-index on test set",
      labels = number_format(accuracy = 0.001)
    ) +
    scale_y_discrete(name = NULL) +
    labs(
      title = "Incremental predictive value of ADNI brain MRI AD L’EPOCH",
      subtitle = "Positive ΔC-index indicates better performance for the full M3 L’EPOCH model.",
      caption = "Intervals are paired bootstrap 95% confidence intervals on the test split."
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.5, hjust = 0),
      axis.text = element_text(color = "black"),
      axis.title.x = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_delta_cindex.pdf"),
    p_delta,
    width = 7.5,
    height = 3.8
  )
  
  ggsave(
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_delta_cindex.png"),
    p_delta,
    width = 7.5,
    height = 3.8,
    dpi = 300
  )
} else {
  delta_tbl <- tibble()
}

# ------------------------------------------------------------
# 8. Test-set KM cumulative incidence plot
# ------------------------------------------------------------

km_train <- pred2 %>%
  filter(split == "train") %>%
  mutate(risk_score_for_quartile = as.numeric(.data[[primary_m3_risk_col]])) %>%
  filter(!is.na(risk_score_for_quartile))

km_test <- pred2 %>%
  filter(split == "test") %>%
  mutate(
    risk_score_for_quartile = as.numeric(.data[[primary_m3_risk_col]]),
    time_years = as.numeric(time_years)
  ) %>%
  filter(!is.na(risk_score_for_quartile), !is.na(time_years), !is.na(event))

if (nrow(km_train) >= 10 && nrow(km_test) >= 10 && sum(km_test$event, na.rm = TRUE) >= 2) {
  
  q_breaks <- quantile(
    km_train$risk_score_for_quartile,
    probs = c(0, 0.25, 0.50, 0.75, 1),
    na.rm = TRUE,
    type = 8
  )
  
  q_breaks <- unique(as.numeric(q_breaks))
  
  if (length(q_breaks) < 5) {
    warning("Training risk-score quartile breaks are not unique. Falling back to ntile within test set.")
    
    km_test <- km_test %>%
      mutate(
        risk_quartile = ntile(risk_score_for_quartile, 4),
        risk_quartile = factor(
          risk_quartile,
          levels = 1:4,
          labels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
        )
      )
    
    quartile_cut_tbl <- tibble(
      method = "fallback_test_ntile_due_to_nonunique_training_breaks",
      q0 = NA_real_,
      q25 = NA_real_,
      q50 = NA_real_,
      q75 = NA_real_,
      q100 = NA_real_
    )
  } else {
    q_breaks[1] <- -Inf
    q_breaks[length(q_breaks)] <- Inf
    
    km_test <- km_test %>%
      mutate(
        risk_quartile = cut(
          risk_score_for_quartile,
          breaks = q_breaks,
          include.lowest = TRUE,
          labels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
        )
      )
    
    quartile_cut_tbl <- tibble(
      method = "training_defined_quartiles",
      q0 = quantile(km_train$risk_score_for_quartile, 0, na.rm = TRUE, type = 8),
      q25 = quantile(km_train$risk_score_for_quartile, 0.25, na.rm = TRUE, type = 8),
      q50 = quantile(km_train$risk_score_for_quartile, 0.50, na.rm = TRUE, type = 8),
      q75 = quantile(km_train$risk_score_for_quartile, 0.75, na.rm = TRUE, type = 8),
      q100 = quantile(km_train$risk_score_for_quartile, 1, na.rm = TRUE, type = 8)
    )
  }
  
  km_test <- km_test %>%
    filter(!is.na(risk_quartile)) %>%
    mutate(
      risk_quartile = factor(
        risk_quartile,
        levels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
      )
    )
  
  readr::write_tsv(
    quartile_cut_tbl,
    file.path(out_dir, "adni_brain_mri_ad_lepoch_training_defined_risk_quartile_cutpoints.tsv")
  )
  
  readr::write_tsv(
    km_test %>%
      select(
        PTID,
        split,
        time_years,
        event,
        event_or_censor_dx,
        all_of(primary_m3_risk_col),
        risk_quartile
      ),
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_risk_quartiles.tsv")
  )
  
  km_fit <- survival::survfit(
    survival::Surv(time_years, event) ~ risk_quartile,
    data = km_test
  )
  
  km_sum <- summary(km_fit)
  
  km_df <- tibble(
    time = km_sum$time,
    n_risk = km_sum$n.risk,
    n_event = km_sum$n.event,
    survival = km_sum$surv,
    lower_survival = km_sum$lower,
    upper_survival = km_sum$upper,
    strata = km_sum$strata
  ) %>%
    mutate(
      risk_quartile = str_replace(as.character(strata), "^risk_quartile=", ""),
      risk_quartile = factor(
        risk_quartile,
        levels = c("Q1 lowest risk", "Q2", "Q3", "Q4 highest risk")
      ),
      cumulative_incidence = 1 - survival,
      cumulative_incidence_low = 1 - upper_survival,
      cumulative_incidence_high = 1 - lower_survival
    )
  
  km_zero <- km_test %>%
    count(risk_quartile, name = "n_risk") %>%
    mutate(
      time = 0,
      n_event = 0,
      survival = 1,
      lower_survival = 1,
      upper_survival = 1,
      strata = paste0("risk_quartile=", risk_quartile),
      cumulative_incidence = 0,
      cumulative_incidence_low = 0,
      cumulative_incidence_high = 0
    )
  
  km_df <- bind_rows(km_zero, km_df) %>%
    arrange(risk_quartile, time)
  
  logrank <- survival::survdiff(
    survival::Surv(time_years, event) ~ risk_quartile,
    data = km_test
  )
  
  logrank_p <- pchisq(logrank$chisq, df = length(logrank$n) - 1, lower.tail = FALSE)
  
  km_group_summary <- km_test %>%
    group_by(risk_quartile) %>%
    summarise(
      n = n(),
      n_events = sum(event, na.rm = TRUE),
      event_rate = n_events / n,
      median_followup_years = median(time_years, na.rm = TRUE),
      mean_risk_score = mean(risk_score_for_quartile, na.rm = TRUE),
      .groups = "drop"
    )
  
  readr::write_tsv(
    km_group_summary,
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_km_risk_quartile_summary.tsv")
  )
  
  readr::write_tsv(
    tibble(logrank_p = logrank_p),
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_km_logrank_p.tsv")
  )
  
  km_palette <- c(
    "Q1 lowest risk" = "#0F766E",
    "Q2" = "#6D597A",
    "Q3" = "#DDAA33",
    "Q4 highest risk" = "#C65D35"
  )
  
  p_km_test <- ggplot(
    km_df,
    aes(
      x = time,
      y = cumulative_incidence,
      color = risk_quartile,
      fill = risk_quartile
    )
  ) +
    geom_ribbon(
      aes(
        ymin = pmax(0, cumulative_incidence_low),
        ymax = pmin(1, cumulative_incidence_high)
      ),
      alpha = 0.13,
      color = NA
    ) +
    geom_step(linewidth = 1.05) +
    scale_color_manual(values = km_palette) +
    scale_fill_manual(values = km_palette) +
    scale_x_continuous(
      name = "Years after baseline",
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      name = "Cumulative incidence of MCI/AD conversion",
      labels = percent_format(accuracy = 0.1),
      expand = expansion(mult = c(0.00, 0.08))
    ) +
    labs(
      title = "Test-set AD L’EPOCH risk separation",
      subtitle = paste0(
        "Kaplan–Meier cumulative incidence by training-defined M3 risk quartile; log-rank ",
        format_p(logrank_p)
      ),
      color = "Risk quartile",
      fill = "Risk quartile",
      caption = "Quartile cutpoints were estimated in the training set and applied to the test set."
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10.5),
      plot.caption = element_text(size = 8.5, hjust = 0),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_km_cumulative_incidence_by_training_risk_quartile.pdf"),
    p_km_test,
    width = 8.5,
    height = 5.8
  )
  
  ggsave(
    file.path(out_dir, "adni_brain_mri_ad_lepoch_test_km_cumulative_incidence_by_training_risk_quartile.png"),
    p_km_test,
    width = 8.5,
    height = 5.8,
    dpi = 300
  )
} else {
  warning("Insufficient training/test data or test-set events for KM quartile plot.")
}

# ------------------------------------------------------------
# 9. Acceleration-years group summary
# ------------------------------------------------------------

accel_col <- "adni_brain_mri_ad_lepoch_acceleration_years"

plot_group_df <- plot_group_df %>%
  mutate(
    acceleration_years = as.numeric(.data[[accel_col]])
  )

group_summary <- plot_group_df %>%
  group_by(conversion_group) %>%
  summarise(
    n_total = n(),
    n_nonmissing = sum(!is.na(acceleration_years)),
    mean = mean(acceleration_years, na.rm = TRUE),
    sd = sd(acceleration_years, na.rm = TRUE),
    median = median(acceleration_years, na.rm = TRUE),
    q1 = quantile(acceleration_years, 0.25, na.rm = TRUE),
    q3 = quantile(acceleration_years, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  group_summary,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_group_summary.tsv")
)

x_label_map <- group_summary %>%
  mutate(
    x_label = paste0(
      as.character(conversion_group),
      "\n(n = ",
      n_nonmissing,
      ")"
    )
  ) %>%
  select(conversion_group, x_label) %>%
  deframe()

# ------------------------------------------------------------
# 10. Pairwise Welch two-sample t-tests using raw P-values
# ------------------------------------------------------------

ttest_df <- plot_group_df %>%
  filter(!is.na(acceleration_years)) %>%
  filter(conversion_group %in% c("Censored / non-event", "CN to MCI", "CN to AD")) %>%
  mutate(
    conversion_group = factor(
      conversion_group,
      levels = c("Censored / non-event", "CN to MCI", "CN to AD")
    )
  ) %>%
  droplevels()

pairwise_comparisons <- tibble(
  group_1 = c("Censored / non-event", "Censored / non-event", "CN to MCI"),
  group_2 = c("CN to MCI", "CN to AD", "CN to AD")
)

pairwise_ttest_tbl <- pairwise_comparisons %>%
  mutate(
    test = map2(group_1, group_2, ~ {
      x <- ttest_df %>%
        filter(conversion_group == .x) %>%
        pull(acceleration_years)
      
      y <- ttest_df %>%
        filter(conversion_group == .y) %>%
        pull(acceleration_years)
      
      x <- x[is.finite(x)]
      y <- y[is.finite(y)]
      
      if (length(x) < 2 || length(y) < 2) {
        return(tibble(
          n_1 = length(x),
          n_2 = length(y),
          mean_1 = mean(x, na.rm = TRUE),
          mean_2 = mean(y, na.rm = TRUE),
          mean_diff_group2_minus_group1 = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
          t_statistic = NA_real_,
          df = NA_real_,
          p_raw = NA_real_
        ))
      }
      
      tt <- t.test(y, x, var.equal = FALSE)
      
      tibble(
        n_1 = length(x),
        n_2 = length(y),
        mean_1 = mean(x, na.rm = TRUE),
        mean_2 = mean(y, na.rm = TRUE),
        mean_diff_group2_minus_group1 = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
        t_statistic = unname(tt$statistic),
        df = unname(tt$parameter),
        p_raw = tt$p.value
      )
    })
  ) %>%
  unnest(test) %>%
  mutate(
    p_adj_bh = p.adjust(p_raw, method = "BH"),
    p_label_raw = case_when(
      is.na(p_raw) ~ "P = NA",
      p_raw < 2.2e-16 ~ "P < 2.2e-16",
      p_raw < 0.001 ~ paste0("P = ", formatC(p_raw, format = "e", digits = 2)),
      TRUE ~ paste0("P = ", signif(p_raw, 3))
    ),
    p_label_bh = case_when(
      is.na(p_adj_bh) ~ "BH P = NA",
      p_adj_bh < 2.2e-16 ~ "BH P < 2.2e-16",
      p_adj_bh < 0.001 ~ paste0("BH P = ", formatC(p_adj_bh, format = "e", digits = 2)),
      TRUE ~ paste0("BH P = ", signif(p_adj_bh, 3))
    )
  )

readr::write_tsv(
  pairwise_ttest_tbl,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_pairwise_ttest_raw_p.tsv")
)

# Optional Wilcoxon sensitivity.
pairwise_wilcox_raw <- pairwise.wilcox.test(
  x = plot_group_df$acceleration_years,
  g = plot_group_df$conversion_group,
  p.adjust.method = "none",
  exact = FALSE
)

pairwise_wilcox_tbl <- as.data.frame(as.table(pairwise_wilcox_raw$p.value)) %>%
  rename(
    group_1 = Var1,
    group_2 = Var2,
    p_raw_wilcox = Freq
  ) %>%
  filter(!is.na(p_raw_wilcox)) %>%
  mutate(
    p_adj_bh_wilcox = p.adjust(p_raw_wilcox, method = "BH")
  )

readr::write_tsv(
  pairwise_wilcox_tbl,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_pairwise_wilcox_raw_p.tsv")
)

# ------------------------------------------------------------
# 11. Plot acceleration years with raw t-test P-values and sample sizes
# ------------------------------------------------------------

group_palette <- c(
  "Censored / non-event" = "#8C6D31",
  "CN to MCI" = "#2A6F9E",
  "CN to AD" = "#A23E48"
)

y_min <- min(ttest_df$acceleration_years, na.rm = TRUE)
y_max <- max(ttest_df$acceleration_years, na.rm = TRUE)
y_range <- y_max - y_min

if (!is.finite(y_range) || y_range == 0) {
  y_range <- 1
}

bracket_tbl <- pairwise_ttest_tbl %>%
  mutate(
    group_1 = factor(group_1, levels = levels(ttest_df$conversion_group)),
    group_2 = factor(group_2, levels = levels(ttest_df$conversion_group)),
    x1 = as.numeric(group_1),
    x2 = as.numeric(group_2),
    y = y_max + row_number() * 0.11 * y_range,
    y_text = y + 0.025 * y_range
  )

p_accel <- ggplot(
  ttest_df,
  aes(
    x = conversion_group,
    y = acceleration_years,
    fill = conversion_group,
    color = conversion_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey50"
  ) +
  geom_violin(
    width = 0.80,
    alpha = 0.28,
    linewidth = 0.35,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.45,
    size = 1.5,
    stroke = 0
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.2,
    fill = "white",
    color = "black"
  ) +
  geom_segment(
    data = bracket_tbl,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.45,
    color = "black"
  ) +
  geom_segment(
    data = bracket_tbl,
    aes(x = x1, xend = x1, y = y - 0.025 * y_range, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.45,
    color = "black"
  ) +
  geom_segment(
    data = bracket_tbl,
    aes(x = x2, xend = x2, y = y - 0.025 * y_range, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.45,
    color = "black"
  ) +
  geom_text(
    data = bracket_tbl,
    aes(x = (x1 + x2) / 2, y = y_text, label = p_label_raw),
    inherit.aes = FALSE,
    size = 3.0,
    color = "black"
  ) +
  scale_fill_manual(values = group_palette, guide = "none") +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_x_discrete(
    name = NULL,
    labels = x_label_map
  ) +
  scale_y_continuous(
    name = "AD L’EPOCH acceleration years",
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.05, 0.30))
  ) +
  labs(
    title = "ADNI brain MRI AD L’EPOCH acceleration by conversion outcome",
    subtitle = "Pairwise Welch two-sample t-tests are annotated using raw P-values.",
    caption = "Sample sizes in x-axis labels are participants with non-missing acceleration years; diamond indicates group median."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_by_conversion_group_with_raw_ttests_and_n.pdf"),
  p_accel,
  width = 8.2,
  height = 6.3
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_by_conversion_group_with_raw_ttests_and_n.png"),
  p_accel,
  width = 8.2,
  height = 6.3,
  dpi = 300
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_by_conversion_group.pdf"),
  p_accel,
  width = 8.2,
  height = 6.3
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_acceleration_years_by_conversion_group.png"),
  p_accel,
  width = 8.2,
  height = 6.3,
  dpi = 300
)

# ------------------------------------------------------------
# 12. Optional: risk score by conversion group with sample sizes
# ------------------------------------------------------------

risk_col <- "adni_brain_mri_ad_lepoch_risk_score"

p_risk <- ggplot(
  plot_group_df,
  aes(
    x = conversion_group,
    y = as.numeric(.data[[risk_col]]),
    fill = conversion_group,
    color = conversion_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey50"
  ) +
  geom_violin(
    width = 0.80,
    alpha = 0.28,
    linewidth = 0.35,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.45,
    size = 1.5,
    stroke = 0
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 3.2,
    fill = "white",
    color = "black"
  ) +
  scale_fill_manual(values = group_palette, guide = "none") +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_x_discrete(
    name = NULL,
    labels = x_label_map
  ) +
  scale_y_continuous(
    name = "AD L’EPOCH Cox risk score",
    labels = number_format(accuracy = 0.01)
  ) +
  labs(
    title = "ADNI brain MRI AD L’EPOCH risk score by conversion outcome",
    subtitle = "Higher values indicate greater model-estimated proximity to MCI/AD conversion.",
    caption = "Sample sizes in x-axis labels are participants with non-missing acceleration years; diamond indicates group median."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_risk_score_by_conversion_group.pdf"),
  p_risk,
  width = 7.8,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_risk_score_by_conversion_group.png"),
  p_risk,
  width = 7.8,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 13. Print key summaries
# ------------------------------------------------------------

message("Outputs saved to: ", out_dir)

message("Group counts:")
print(
  plot_group_df %>%
    count(conversion_group, event_or_censor_dx, event)
)

message("Acceleration-years group summary:")
print(group_summary)

message("Pairwise Welch t-test results, raw P-values annotated in figure:")
print(pairwise_ttest_tbl)

message("C-index summary:")
print(cindex_tbl)

if (exists("delta_tbl")) {
  message("Test-set delta C-index:")
  print(delta_tbl)
}

if (exists("km_group_summary")) {
  message("Test-set KM risk quartile summary:")
  print(km_group_summary)
}