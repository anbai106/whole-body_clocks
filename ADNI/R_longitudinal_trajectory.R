# ============================================================
# Plot longitudinal ADNI brain MRI AD L'EPOCH trajectories
# CN scans only, before MCI/AD conversion
#
# Revised:
#   1) p_spaghetti includes participant sample size in facet labels.
#   2) p_spaghetti includes population-level linear trend line.
#   3) p_spaghetti annotates slope beta, standardized beta, and raw P-value.
#
# Input:
#   adni_brain_mri_ad_lepoch_longitudinal_cn_only_predictions.tsv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

long_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/results_brain_mri_ad_lepoch_longitudinal_cn_only"

pred_file <- file.path(
  long_dir,
  "adni_brain_mri_ad_lepoch_longitudinal_cn_only_predictions.tsv"
)

out_dir <- file.path(long_dir, "trajectory_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 2. Read longitudinal predictions
# ------------------------------------------------------------

df <- readr::read_tsv(pred_file, show_col_types = FALSE)

required_cols <- c(
  "PTID",
  "Visit_Code",
  "scan_dx",
  "years_since_selected_baseline",
  "conversion_group",
  "scan_relation_to_event",
  "event_from_selected_baseline",
  "adni_brain_mri_ad_lepoch_acceleration_years",
  "adni_brain_mri_ad_lepoch_acceleration_z",
  "adni_brain_mri_ad_lepoch_risk_score"
)

missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ------------------------------------------------------------
# 3. Clean variables
# ------------------------------------------------------------

plot_df <- df %>%
  mutate(
    PTID = as.character(PTID),
    years_since_selected_baseline = as.numeric(years_since_selected_baseline),
    acceleration_years = as.numeric(adni_brain_mri_ad_lepoch_acceleration_years),
    acceleration_z = as.numeric(adni_brain_mri_ad_lepoch_acceleration_z),
    risk_score = as.numeric(adni_brain_mri_ad_lepoch_risk_score),
    
    scan_dx = toupper(as.character(scan_dx)),
    scan_relation_to_event = as.character(scan_relation_to_event),
    conversion_group = as.character(conversion_group),
    
    conversion_group = case_when(
      conversion_group %in% c("CN to MCI", "CN to AD", "Censored / non-event") ~ conversion_group,
      TRUE ~ "Other"
    ),
    conversion_group = factor(
      conversion_group,
      levels = c("Censored / non-event", "CN to MCI", "CN to AD", "Other")
    ),
    
    scan_relation_to_event = factor(
      scan_relation_to_event,
      levels = c(
        "selected_baseline",
        "pre_event_CN_followup",
        "censored_CN_followup"
      )
    ),
    
    scan_dx = factor(scan_dx, levels = c("CN", "MCI", "AD"))
  ) %>%
  filter(
    !is.na(years_since_selected_baseline),
    !is.na(acceleration_years),
    conversion_group %in% c("Censored / non-event", "CN to MCI", "CN to AD"),
    scan_dx == "CN",
    scan_relation_to_event %in% c(
      "selected_baseline",
      "pre_event_CN_followup",
      "censored_CN_followup"
    )
  ) %>%
  arrange(PTID, years_since_selected_baseline) %>%
  droplevels()

# Safety check: no MCI/AD scans should remain.
bad_scans <- plot_df %>%
  filter(scan_dx %in% c("MCI", "AD"))

if (nrow(bad_scans) > 0) {
  stop("MCI/AD scans are present in the CN-only longitudinal prediction file.")
}

# Keep people with at least two qualified CN scans for within-subject trajectory plots.
longitudinal_subjects <- plot_df %>%
  group_by(PTID) %>%
  summarise(
    n_scans = n(),
    conversion_group = first(conversion_group),
    .groups = "drop"
  ) %>%
  filter(n_scans >= 2)

traj_df <- plot_df %>%
  semi_join(longitudinal_subjects, by = "PTID") %>%
  droplevels()

# ------------------------------------------------------------
# 4. Palettes
# ------------------------------------------------------------

group_palette <- c(
  "Censored / non-event" = "#8C6D31",
  "CN to MCI" = "#2A6F9E",
  "CN to AD" = "#A23E48"
)

relation_shape <- c(
  "selected_baseline" = 21,
  "pre_event_CN_followup" = 16,
  "censored_CN_followup" = 16
)

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

format_beta <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    abs(x) < 0.001 ~ formatC(x, format = "e", digits = 2),
    TRUE ~ sprintf("%.3f", x)
  )
}

# Fit population-level trend within each conversion group.
# Primary model:
#   random-intercept mixed model if nlme is available:
#      acceleration_years ~ years_since_selected_baseline + random intercept for PTID
#
# Fallback:
#   ordinary linear model.
#
# Effect size:
#   beta = slope in acceleration-years per year
#   std_beta = beta * SD(time) / SD(acceleration_years)
fit_group_trend <- function(d, group_name) {
  d <- d %>%
    filter(
      is.finite(years_since_selected_baseline),
      is.finite(acceleration_years)
    )
  
  n_subjects <- d %>% distinct(PTID) %>% nrow()
  n_scans <- nrow(d)
  
  if (
    n_subjects < 2 ||
    n_scans < 4 ||
    n_distinct(d$years_since_selected_baseline) < 2
  ) {
    return(tibble(
      conversion_group = as.character(group_name),
      n_subjects = n_subjects,
      n_scans = n_scans,
      model = "insufficient data",
      beta_slope_per_year = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      std_beta = NA_real_
    ))
  }
  
  sd_x <- sd(d$years_since_selected_baseline, na.rm = TRUE)
  sd_y <- sd(d$acceleration_years, na.rm = TRUE)
  
  # Try mixed model first.
  if (requireNamespace("nlme", quietly = TRUE)) {
    fit_lme <- tryCatch(
      {
        nlme::lme(
          acceleration_years ~ years_since_selected_baseline,
          random = ~ 1 | PTID,
          data = d,
          method = "REML",
          na.action = na.omit,
          control = nlme::lmeControl(
            opt = "optim",
            msMaxIter = 200,
            returnObject = TRUE
          )
        )
      },
      error = function(e) NULL
    )
    
    if (!is.null(fit_lme)) {
      tt <- summary(fit_lme)$tTable
      
      if ("years_since_selected_baseline" %in% rownames(tt)) {
        beta <- unname(tt["years_since_selected_baseline", "Value"])
        se <- unname(tt["years_since_selected_baseline", "Std.Error"])
        p <- unname(tt["years_since_selected_baseline", "p-value"])
        df_lme <- unname(tt["years_since_selected_baseline", "DF"])
        
        ci_low <- beta - qt(0.975, df = df_lme) * se
        ci_high <- beta + qt(0.975, df = df_lme) * se
        
        std_beta <- ifelse(
          is.finite(sd_x) && is.finite(sd_y) && sd_y > 0,
          beta * sd_x / sd_y,
          NA_real_
        )
        
        return(tibble(
          conversion_group = as.character(group_name),
          n_subjects = n_subjects,
          n_scans = n_scans,
          model = "linear mixed model: random intercept for PTID",
          beta_slope_per_year = beta,
          se = se,
          ci_low = ci_low,
          ci_high = ci_high,
          p_value = p,
          std_beta = std_beta
        ))
      }
    }
  }
  
  # Fallback ordinary linear model.
  fit_lm <- lm(acceleration_years ~ years_since_selected_baseline, data = d)
  sm <- summary(fit_lm)
  tt <- sm$coefficients
  
  beta <- unname(tt["years_since_selected_baseline", "Estimate"])
  se <- unname(tt["years_since_selected_baseline", "Std. Error"])
  p <- unname(tt["years_since_selected_baseline", "Pr(>|t|)"])
  df_lm <- fit_lm$df.residual
  
  ci_low <- beta - qt(0.975, df = df_lm) * se
  ci_high <- beta + qt(0.975, df = df_lm) * se
  
  std_beta <- ifelse(
    is.finite(sd_x) && is.finite(sd_y) && sd_y > 0,
    beta * sd_x / sd_y,
    NA_real_
  )
  
  tibble(
    conversion_group = as.character(group_name),
    n_subjects = n_subjects,
    n_scans = n_scans,
    model = "ordinary linear model fallback",
    beta_slope_per_year = beta,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high,
    p_value = p,
    std_beta = std_beta
  )
}

# ------------------------------------------------------------
# 6. Basic summaries
# ------------------------------------------------------------

scan_summary <- plot_df %>%
  group_by(conversion_group, scan_relation_to_event, scan_dx) %>%
  summarise(
    n_scans = n(),
    n_subjects = n_distinct(PTID),
    mean_acceleration_years = mean(acceleration_years, na.rm = TRUE),
    sd_acceleration_years = sd(acceleration_years, na.rm = TRUE),
    median_acceleration_years = median(acceleration_years, na.rm = TRUE),
    .groups = "drop"
  )

subject_summary <- plot_df %>%
  group_by(PTID) %>%
  summarise(
    conversion_group = first(conversion_group),
    event_from_selected_baseline = first(event_from_selected_baseline),
    n_scans = n(),
    min_year = min(years_since_selected_baseline, na.rm = TRUE),
    max_year = max(years_since_selected_baseline, na.rm = TRUE),
    followup_span_years = max_year - min_year,
    baseline_acceleration_years = acceleration_years[which.min(abs(years_since_selected_baseline))][1],
    last_acceleration_years = acceleration_years[which.max(years_since_selected_baseline)][1],
    change_last_minus_baseline = last_acceleration_years - baseline_acceleration_years,
    .groups = "drop"
  )

readr::write_tsv(
  scan_summary,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_longitudinal_scan_summary.tsv")
)

readr::write_tsv(
  subject_summary,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_longitudinal_subject_summary.tsv")
)

# Sample sizes for facet labels, using only subjects with >=2 scans in p_spaghetti.
facet_n_tbl <- traj_df %>%
  group_by(conversion_group) %>%
  summarise(
    n_subjects = n_distinct(PTID),
    n_scans = n(),
    .groups = "drop"
  ) %>%
  mutate(
    facet_label = paste0(
      as.character(conversion_group),
      "\nN = ",
      n_subjects,
      " participants; ",
      n_scans,
      " CN scans"
    )
  )

facet_label_map <- facet_n_tbl %>%
  select(conversion_group, facet_label) %>%
  mutate(conversion_group = as.character(conversion_group)) %>%
  deframe()

readr::write_tsv(
  facet_n_tbl,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_spaghetti_facet_sample_sizes.tsv")
)

# ------------------------------------------------------------
# 7. Population-level trend models for p_spaghetti
# ------------------------------------------------------------

trend_tbl <- traj_df %>%
  group_split(conversion_group) %>%
  map_dfr(~ fit_group_trend(.x, unique(.x$conversion_group)[1])) %>%
  mutate(
    conversion_group = factor(
      conversion_group,
      levels = c("Censored / non-event", "CN to MCI", "CN to AD")
    ),
    trend_label = paste0(
      "N = ", n_subjects, "; scans = ", n_scans,
      "\nβ = ", format_beta(beta_slope_per_year), " y/y",
      "\nstd β = ", format_beta(std_beta),
      "\n", format_p(p_value)
    )
  )

readr::write_tsv(
  trend_tbl,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_trend_models.tsv")
)

# Annotation positions for each facet.
trend_annot_tbl <- traj_df %>%
  group_by(conversion_group) %>%
  summarise(
    x_min = min(years_since_selected_baseline, na.rm = TRUE),
    x_max = max(years_since_selected_baseline, na.rm = TRUE),
    y_min = min(acceleration_years, na.rm = TRUE),
    y_max = max(acceleration_years, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_range = if_else(is.finite(x_max - x_min) & (x_max - x_min) > 0, x_max - x_min, 1),
    y_range = if_else(is.finite(y_max - y_min) & (y_max - y_min) > 0, y_max - y_min, 1),
    x = x_min + 0.04 * x_range,
    y = y_max - 0.04 * y_range
  ) %>%
  left_join(trend_tbl, by = "conversion_group")

# ------------------------------------------------------------
# 8. Plot 1: individual CN-only trajectories by conversion group
# ------------------------------------------------------------

p_spaghetti <- ggplot(
  traj_df,
  aes(
    x = years_since_selected_baseline,
    y = acceleration_years,
    group = PTID,
    color = conversion_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_line(
    alpha = 0.25,
    linewidth = 0.55
  ) +
  geom_point(
    aes(shape = scan_relation_to_event),
    alpha = 0.78,
    size = 1.9,
    stroke = 0.35
  ) +
  
  # Linear population-level fitting line.
  # The annotated beta/P-value comes from the mixed model above if nlme is available.
  geom_smooth(
    aes(group = conversion_group),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.15,
    alpha = 0.16
  ) +
  
  # Trend annotation.
  geom_label(
    data = trend_annot_tbl,
    aes(
      x = x,
      y = y,
      label = trend_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    label.size = 0.25,
    label.padding = unit(0.18, "lines"),
    fill = "white",
    alpha = 0.88,
    color = "black"
  ) +
  
  facet_wrap(
    ~ conversion_group,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(facet_label_map)
  ) +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_shape_manual(values = relation_shape, drop = FALSE) +
  scale_x_continuous(
    name = "Years since selected CN imaging baseline",
    breaks = pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    name = "AD L’EPOCH acceleration years",
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "Within-subject longitudinal AD L’EPOCH trajectories",
    subtitle = "CN-labeled MRI scans only; scans at or after MCI/AD conversion are excluded.",
    shape = "Scan relation",
    caption = paste0(
      "Each thin line is one participant with at least two CN MRI scans. ",
      "Thick line shows linear population trend. ",
      "β is acceleration-years change per calendar year; P tests whether β differs from zero."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_individual_trajectories_acceleration_years_with_trend_stats.pdf"),
  p_spaghetti,
  width = 12.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_individual_trajectories_acceleration_years_with_trend_stats.png"),
  p_spaghetti,
  width = 12.5,
  height = 5.6,
  dpi = 300
)

# Also save to original convenient filename.
ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_individual_trajectories_acceleration_years.pdf"),
  p_spaghetti,
  width = 12.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_individual_trajectories_acceleration_years.png"),
  p_spaghetti,
  width = 12.5,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Plot 2: population-level mean trajectory by time bin
# ------------------------------------------------------------

bin_width <- 0.5

pop_df <- plot_df %>%
  mutate(
    time_bin = round(years_since_selected_baseline / bin_width) * bin_width
  ) %>%
  group_by(conversion_group, time_bin) %>%
  summarise(
    n_scans = n(),
    n_subjects = n_distinct(PTID),
    mean_acceleration_years = mean(acceleration_years, na.rm = TRUE),
    sd_acceleration_years = sd(acceleration_years, na.rm = TRUE),
    se_acceleration_years = sd_acceleration_years / sqrt(n_scans),
    ci_low = mean_acceleration_years - 1.96 * se_acceleration_years,
    ci_high = mean_acceleration_years + 1.96 * se_acceleration_years,
    .groups = "drop"
  ) %>%
  filter(n_scans >= 2)

readr::write_tsv(
  pop_df,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_timebin_summary.tsv")
)

p_population <- ggplot(
  pop_df,
  aes(
    x = time_bin,
    y = mean_acceleration_years,
    color = conversion_group,
    fill = conversion_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(aes(size = n_subjects), alpha = 0.90) +
  scale_color_manual(values = group_palette) +
  scale_fill_manual(values = group_palette) +
  scale_size_continuous(
    name = "N subjects",
    range = c(1.8, 4.2)
  ) +
  scale_x_continuous(
    name = "Years since selected CN imaging baseline",
    breaks = pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    name = "Mean AD L’EPOCH acceleration years",
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Population-level AD L’EPOCH trajectory",
    subtitle = paste0("CN scans only; mean acceleration years in ", bin_width, "-year bins."),
    color = NULL,
    fill = NULL,
    caption = "Shaded region is approximate 95% CI. Bins with fewer than 2 scans are not shown."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_trajectory_acceleration_years.pdf"),
  p_population,
  width = 8.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_trajectory_acceleration_years.png"),
  p_population,
  width = 8.5,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Plot 3: within-subject delta from selected baseline
# ------------------------------------------------------------

baseline_tbl <- plot_df %>%
  group_by(PTID) %>%
  arrange(abs(years_since_selected_baseline), .by_group = TRUE) %>%
  summarise(
    baseline_acceleration_years = first(acceleration_years),
    baseline_year = first(years_since_selected_baseline),
    .groups = "drop"
  )

delta_df <- plot_df %>%
  left_join(baseline_tbl, by = "PTID") %>%
  mutate(
    delta_acceleration_years = acceleration_years - baseline_acceleration_years
  )

delta_pop_df <- delta_df %>%
  mutate(
    time_bin = round(years_since_selected_baseline / bin_width) * bin_width
  ) %>%
  group_by(conversion_group, time_bin) %>%
  summarise(
    n_scans = n(),
    n_subjects = n_distinct(PTID),
    mean_delta_acceleration_years = mean(delta_acceleration_years, na.rm = TRUE),
    sd_delta_acceleration_years = sd(delta_acceleration_years, na.rm = TRUE),
    se_delta_acceleration_years = sd_delta_acceleration_years / sqrt(n_scans),
    ci_low = mean_delta_acceleration_years - 1.96 * se_delta_acceleration_years,
    ci_high = mean_delta_acceleration_years + 1.96 * se_delta_acceleration_years,
    .groups = "drop"
  ) %>%
  filter(n_scans >= 2)

readr::write_tsv(
  delta_df,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_delta_from_baseline.tsv")
)

readr::write_tsv(
  delta_pop_df,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_delta_timebin_summary.tsv")
)

p_delta_population <- ggplot(
  delta_pop_df,
  aes(
    x = time_bin,
    y = mean_delta_acceleration_years,
    color = conversion_group,
    fill = conversion_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(aes(size = n_subjects), alpha = 0.90) +
  scale_color_manual(values = group_palette) +
  scale_fill_manual(values = group_palette) +
  scale_size_continuous(
    name = "N subjects",
    range = c(1.8, 4.2)
  ) +
  scale_x_continuous(
    name = "Years since selected CN imaging baseline",
    breaks = pretty_breaks(n = 6)
  ) +
  scale_y_continuous(
    name = "Mean change in AD L’EPOCH acceleration years",
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Within-subject change in AD L’EPOCH acceleration",
    subtitle = "CN scans only; each subject is centered at selected CN imaging baseline.",
    color = NULL,
    fill = NULL,
    caption = "Positive values indicate increasing AD L’EPOCH acceleration relative to the participant's baseline scan."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    legend.position = "top",
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_delta_from_baseline.pdf"),
  p_delta_population,
  width = 8.5,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_population_delta_from_baseline.png"),
  p_delta_population,
  width = 8.5,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Plot 4: subject-level linear slopes
# ------------------------------------------------------------

slope_tbl <- plot_df %>%
  group_by(PTID) %>%
  filter(n() >= 2) %>%
  summarise(
    conversion_group = first(conversion_group),
    n_scans = n(),
    followup_span_years = max(years_since_selected_baseline, na.rm = TRUE) -
      min(years_since_selected_baseline, na.rm = TRUE),
    slope_acceleration_years_per_year = {
      fit <- lm(acceleration_years ~ years_since_selected_baseline)
      unname(coef(fit)[["years_since_selected_baseline"]])
    },
    intercept = {
      fit <- lm(acceleration_years ~ years_since_selected_baseline)
      unname(coef(fit)[["(Intercept)"]])
    },
    .groups = "drop"
  ) %>%
  filter(is.finite(slope_acceleration_years_per_year))

readr::write_tsv(
  slope_tbl,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_subject_level_slopes.tsv")
)

slope_summary <- slope_tbl %>%
  group_by(conversion_group) %>%
  summarise(
    n_subjects = n(),
    mean_slope = mean(slope_acceleration_years_per_year, na.rm = TRUE),
    sd_slope = sd(slope_acceleration_years_per_year, na.rm = TRUE),
    median_slope = median(slope_acceleration_years_per_year, na.rm = TRUE),
    q1 = quantile(slope_acceleration_years_per_year, 0.25, na.rm = TRUE),
    q3 = quantile(slope_acceleration_years_per_year, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  slope_summary,
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_subject_level_slope_summary.tsv")
)

slope_x_labels <- slope_summary %>%
  mutate(
    label = paste0(as.character(conversion_group), "\n(n = ", n_subjects, ")")
  ) %>%
  select(conversion_group, label) %>%
  deframe()

p_slope <- ggplot(
  slope_tbl,
  aes(
    x = conversion_group,
    y = slope_acceleration_years_per_year,
    fill = conversion_group,
    color = conversion_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey55"
  ) +
  geom_violin(
    width = 0.80,
    alpha = 0.25,
    trim = FALSE,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.85,
    linewidth = 0.45
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.60,
    size = 1.8
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
    labels = slope_x_labels
  ) +
  scale_y_continuous(
    name = "Subject-level slope\n(acceleration years per year)",
    labels = number_format(accuracy = 0.1)
  ) +
  labs(
    title = "Subject-level AD L’EPOCH trajectory slopes",
    subtitle = "CN scans only; each point is a participant with at least two qualified CN MRI scans.",
    caption = "Positive slopes indicate increasing AD L’EPOCH acceleration before conversion or censoring."
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10.5),
    plot.caption = element_text(size = 8.5, hjust = 0),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_subject_level_slope_distribution.pdf"),
  p_slope,
  width = 7.8,
  height = 5.6
)

ggsave(
  file.path(out_dir, "adni_brain_mri_ad_lepoch_cn_only_subject_level_slope_distribution.png"),
  p_slope,
  width = 7.8,
  height = 5.6,
  dpi = 300
)

# ------------------------------------------------------------
# 12. Print summaries
# ------------------------------------------------------------

message("Outputs saved to: ", out_dir)

message("CN-only longitudinal scan summary:")
print(scan_summary)

message("Subject summary:")
print(subject_summary %>% count(conversion_group, n_scans))

message("Population-level trend models for p_spaghetti:")
print(trend_tbl)

message("Slope summary:")
print(slope_summary)