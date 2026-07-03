# ============================================================
# Longitudinal pulmonary L'EPOCH (de)acceleration analysis
# across UKB proteomics instances
#
# Instance 0: baseline, 2006-2010
# Instance 2: imaging visit, 2014+
# Instance 3: repeat imaging, 2019+
#
# Main y-axis variable:
# pulmonary_proteomics_mortality_clock_acceleration_years
#
# Groups:
# 1) Whole sample
# 2) Cases: participants with death event during administrative follow-up
# 3) Non-cases
#
# Note:
# Cross-instance calibration may exaggerate absolute gaps.
# Therefore, interpret the plots as longitudinal shifts in L'EPOCH
# acceleration rather than exact biological age acceleration.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(scales)
})

# -----------------------------
# 1. Input files
# -----------------------------

baseline_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Pulmonary_proteomics_mortality_clock/pulmonary_proteomics_mortality_clock_predictions.tsv"

instance2_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/Pulmonary/pulmonary_proteomics_mortality_clock_apply_instance_2_0_predictions.tsv"

instance3_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/Pulmonary/pulmonary_proteomics_mortality_clock_apply_instance_3_0_predictions.tsv"

outdir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/Pulmonary/longitudinal_acceleration_analysis_by_case_status"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

accel_col <- "pulmonary_proteomics_mortality_clock_acceleration_years"
z_col     <- "pulmonary_proteomics_mortality_clock_acceleration_z"
clock_col <- "pulmonary_proteomics_mortality_clock_age_years"

# -----------------------------
# 2. Van Gogh-style palette
# -----------------------------

instance_palette <- c(
  "Instance 0\nBaseline"       = "#2F4B7C",
  "Instance 2\nImaging"        = "#E0B43B",
  "Instance 3\nRepeat imaging" = "#3A7D6B"
)

instance_palette_fill <- alpha(instance_palette, 0.75)

group_palette <- c(
  "Whole sample" = "#2F4B7C",
  "Cases"        = "#B85C38",
  "Non-cases"    = "#3A7D6B"
)

# -----------------------------
# 3. Helper functions
# -----------------------------

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-300) return("<1e-300")
  if (p < 0.001) return(format(p, scientific = TRUE, digits = 2))
  formatC(p, format = "f", digits = 3)
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

safe_wilcox_p <- function(x, y) {
  tryCatch(
    wilcox.test(x, y, paired = TRUE)$p.value,
    error = function(e) NA_real_
  )
}

parse_event_column <- function(x) {
  x_chr <- as.character(x)
  
  case_when(
    x_chr %in% c("TRUE", "True", "true", "T", "1", "1.0") ~ TRUE,
    x_chr %in% c("FALSE", "False", "false", "F", "0", "0.0") ~ FALSE,
    TRUE ~ NA
  )
}

# -----------------------------
# 4. Read and standardize each instance
# -----------------------------

read_clock_instance <- function(file, instance_label, visit_order, expected_instance_year) {
  x <- fread(file)
  
  if (!"participant_id" %in% names(x)) {
    stop("participant_id column not found in: ", file)
  }
  
  if (!accel_col %in% names(x)) {
    stop(accel_col, " not found in: ", file)
  }
  
  optional_cols <- c(
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "event",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    accel_col,
    z_col,
    clock_col
  )
  
  for (cc in optional_cols) {
    if (!cc %in% names(x)) {
      x[[cc]] <- NA
    }
  }
  
  keep_cols <- c(
    "participant_id",
    "sample_date",
    "death_date",
    "admin_censor_date",
    "end_date",
    "event",
    "age_at_baseline",
    "age_at_imaging",
    "sex",
    accel_col,
    z_col,
    clock_col
  )
  
  x <- x[, ..keep_cols]
  
  x <- x %>%
    mutate(
      application_instance = instance_label,
      visit_order = visit_order,
      expected_instance_year = expected_instance_year,
      
      sample_date = as.Date(sample_date),
      death_date = as.Date(death_date),
      admin_censor_date = as.Date(admin_censor_date),
      end_date = as.Date(end_date),
      
      sample_year = as.numeric(format(sample_date, "%Y")),
      
      chronological_age = as.numeric(age_at_baseline),
      clock_acceleration_years = as.numeric(.data[[accel_col]]),
      clock_acceleration_z = as.numeric(.data[[z_col]]),
      clock_age_years = as.numeric(.data[[clock_col]]),
      
      event_from_column = parse_event_column(event),
      event_from_dates = case_when(
        !is.na(death_date) & !is.na(sample_date) & !is.na(admin_censor_date) ~
          death_date > sample_date & death_date <= admin_censor_date,
        !is.na(death_date) & !is.na(sample_date) & is.na(admin_censor_date) ~
          death_date > sample_date,
        TRUE ~ NA
      ),
      event = case_when(
        !is.na(event_from_column) ~ event_from_column,
        !is.na(event_from_dates) ~ event_from_dates,
        TRUE ~ FALSE
      )
    ) %>%
    select(-event_from_column, -event_from_dates)
  
  x
}

# -----------------------------
# 5. Read instance 0, 2, and 3
# -----------------------------

dat0 <- read_clock_instance(
  file = baseline_file,
  instance_label = "0_0",
  visit_order = 0,
  expected_instance_year = 2008
)

dat2 <- read_clock_instance(
  file = instance2_file,
  instance_label = "2_0",
  visit_order = 1,
  expected_instance_year = 2014
)

dat3 <- read_clock_instance(
  file = instance3_file,
  instance_label = "3_0",
  visit_order = 2,
  expected_instance_year = 2019
)

dat_all <- bind_rows(dat0, dat2, dat3) %>%
  mutate(
    application_instance = factor(application_instance, levels = c("0_0", "2_0", "3_0")),
    instance_label = recode(
      as.character(application_instance),
      "0_0" = "Instance 0\nBaseline",
      "2_0" = "Instance 2\nImaging",
      "3_0" = "Instance 3\nRepeat imaging"
    ),
    instance_label = factor(
      instance_label,
      levels = c(
        "Instance 0\nBaseline",
        "Instance 2\nImaging",
        "Instance 3\nRepeat imaging"
      )
    )
  )

# -----------------------------
# 6. Define common population across all 3 instances
# -----------------------------

common_ids <- dat_all %>%
  filter(!is.na(clock_acceleration_years)) %>%
  distinct(participant_id, application_instance) %>%
  count(participant_id, name = "n_instances") %>%
  filter(n_instances == 3) %>%
  pull(participant_id)

dat_common <- dat_all %>%
  filter(participant_id %in% common_ids) %>%
  filter(!is.na(clock_acceleration_years)) %>%
  arrange(participant_id, visit_order)

cat("N in instance 0:", n_distinct(dat0$participant_id), "\n")
cat("N in instance 2:", n_distinct(dat2$participant_id), "\n")
cat("N in instance 3:", n_distinct(dat3$participant_id), "\n")
cat("N common across 0, 2, and 3:", length(common_ids), "\n")

# -----------------------------
# 7. Define case/non-case status
# -----------------------------
# Case = death event during administrative follow-up.
# The case definition is based preferentially on instance 0 because the
# baseline prediction file was generated with the mortality survival outcome.

case_status_tbl <- dat_common %>%
  filter(application_instance == "0_0") %>%
  group_by(participant_id) %>%
  summarise(
    event_baseline = any(event %in% TRUE, na.rm = TRUE),
    death_date_baseline = suppressWarnings(min(death_date, na.rm = TRUE)),
    sample_date_baseline = suppressWarnings(min(sample_date, na.rm = TRUE)),
    admin_censor_date_baseline = suppressWarnings(min(admin_censor_date, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    death_date_baseline = as.Date(ifelse(is.infinite(death_date_baseline), NA, death_date_baseline), origin = "1970-01-01"),
    sample_date_baseline = as.Date(ifelse(is.infinite(sample_date_baseline), NA, sample_date_baseline), origin = "1970-01-01"),
    admin_censor_date_baseline = as.Date(ifelse(is.infinite(admin_censor_date_baseline), NA, admin_censor_date_baseline), origin = "1970-01-01"),
    
    event_from_dates_baseline = case_when(
      !is.na(death_date_baseline) & !is.na(sample_date_baseline) & !is.na(admin_censor_date_baseline) ~
        death_date_baseline > sample_date_baseline & death_date_baseline <= admin_censor_date_baseline,
      !is.na(death_date_baseline) & !is.na(sample_date_baseline) & is.na(admin_censor_date_baseline) ~
        death_date_baseline > sample_date_baseline,
      TRUE ~ FALSE
    ),
    
    event_any = event_baseline | event_from_dates_baseline,
    case_status = if_else(event_any, "Cases", "Non-cases")
  ) %>%
  select(participant_id, event_any, case_status)

dat_common <- dat_common %>%
  left_join(case_status_tbl, by = "participant_id") %>%
  mutate(
    case_status = factor(case_status, levels = c("Cases", "Non-cases"))
  )

cat("N cases among common population:", n_distinct(dat_common$participant_id[dat_common$case_status == "Cases"]), "\n")
cat("N non-cases among common population:", n_distinct(dat_common$participant_id[dat_common$case_status == "Non-cases"]), "\n")

# Create plotting/analysis dataset with whole sample + subgroup rows
dat_plot <- bind_rows(
  dat_common %>% mutate(analysis_group = "Whole sample"),
  dat_common %>% mutate(analysis_group = as.character(case_status))
) %>%
  mutate(
    analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
  )

fwrite(
  dat_common,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_common_population_long.tsv"),
  sep = "\t"
)

fwrite(
  dat_plot,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_plot_dataset_with_groups.tsv"),
  sep = "\t"
)

# -----------------------------
# 8. Summary statistics by group and instance
# -----------------------------

summary_tbl <- dat_plot %>%
  group_by(analysis_group, application_instance, instance_label) %>%
  summarise(
    n_rows = n(),
    n_participants = n_distinct(participant_id),
    
    mean_acceleration_years = mean(clock_acceleration_years, na.rm = TRUE),
    sd_acceleration_years = sd(clock_acceleration_years, na.rm = TRUE),
    median_acceleration_years = median(clock_acceleration_years, na.rm = TRUE),
    q25_acceleration_years = quantile(clock_acceleration_years, 0.25, na.rm = TRUE),
    q75_acceleration_years = quantile(clock_acceleration_years, 0.75, na.rm = TRUE),
    min_acceleration_years = min(clock_acceleration_years, na.rm = TRUE),
    max_acceleration_years = max(clock_acceleration_years, na.rm = TRUE),
    
    mean_chronological_age = mean(chronological_age, na.rm = TRUE),
    sd_chronological_age = sd(chronological_age, na.rm = TRUE),
    mean_clock_age = mean(clock_age_years, na.rm = TRUE),
    sd_clock_age = sd(clock_age_years, na.rm = TRUE),
    
    .groups = "drop"
  )

fwrite(
  summary_tbl,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_instance_summary_by_group.tsv"),
  sep = "\t"
)

print(summary_tbl)

# -----------------------------
# 9. Wide table and paired deltas by group
# -----------------------------

dat_wide <- dat_plot %>%
  select(
    participant_id,
    analysis_group,
    application_instance,
    clock_acceleration_years,
    chronological_age,
    clock_age_years
  ) %>%
  pivot_wider(
    names_from = application_instance,
    values_from = c(clock_acceleration_years, chronological_age, clock_age_years)
  ) %>%
  mutate(
    delta_accel_2_minus_0 = clock_acceleration_years_2_0 - clock_acceleration_years_0_0,
    delta_accel_3_minus_2 = clock_acceleration_years_3_0 - clock_acceleration_years_2_0,
    delta_accel_3_minus_0 = clock_acceleration_years_3_0 - clock_acceleration_years_0_0,
    
    delta_chrono_age_2_minus_0 = chronological_age_2_0 - chronological_age_0_0,
    delta_chrono_age_3_minus_2 = chronological_age_3_0 - chronological_age_2_0,
    delta_chrono_age_3_minus_0 = chronological_age_3_0 - chronological_age_0_0,
    
    delta_clock_age_2_minus_0 = clock_age_years_2_0 - clock_age_years_0_0,
    delta_clock_age_3_minus_2 = clock_age_years_3_0 - clock_age_years_2_0,
    delta_clock_age_3_minus_0 = clock_age_years_3_0 - clock_age_years_0_0
  )

fwrite(
  dat_wide,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_common_population_wide_deltas_by_group.tsv"),
  sep = "\t"
)

delta_summary_tbl <- dat_wide %>%
  group_by(analysis_group) %>%
  summarise(
    n = n(),
    
    mean_accel_0 = mean(clock_acceleration_years_0_0, na.rm = TRUE),
    mean_accel_2 = mean(clock_acceleration_years_2_0, na.rm = TRUE),
    mean_accel_3 = mean(clock_acceleration_years_3_0, na.rm = TRUE),
    
    mean_delta_accel_2_minus_0 = mean(delta_accel_2_minus_0, na.rm = TRUE),
    sd_delta_accel_2_minus_0 = sd(delta_accel_2_minus_0, na.rm = TRUE),
    median_delta_accel_2_minus_0 = median(delta_accel_2_minus_0, na.rm = TRUE),
    p_wilcox_accel_2_vs_0 = safe_wilcox_p(clock_acceleration_years_2_0, clock_acceleration_years_0_0),
    
    mean_delta_accel_3_minus_2 = mean(delta_accel_3_minus_2, na.rm = TRUE),
    sd_delta_accel_3_minus_2 = sd(delta_accel_3_minus_2, na.rm = TRUE),
    median_delta_accel_3_minus_2 = median(delta_accel_3_minus_2, na.rm = TRUE),
    p_wilcox_accel_3_vs_2 = safe_wilcox_p(clock_acceleration_years_3_0, clock_acceleration_years_2_0),
    
    mean_delta_accel_3_minus_0 = mean(delta_accel_3_minus_0, na.rm = TRUE),
    sd_delta_accel_3_minus_0 = sd(delta_accel_3_minus_0, na.rm = TRUE),
    median_delta_accel_3_minus_0 = median(delta_accel_3_minus_0, na.rm = TRUE),
    p_wilcox_accel_3_vs_0 = safe_wilcox_p(clock_acceleration_years_3_0, clock_acceleration_years_0_0),
    
    mean_delta_chrono_age_2_minus_0 = mean(delta_chrono_age_2_minus_0, na.rm = TRUE),
    mean_delta_chrono_age_3_minus_2 = mean(delta_chrono_age_3_minus_2, na.rm = TRUE),
    mean_delta_chrono_age_3_minus_0 = mean(delta_chrono_age_3_minus_0, na.rm = TRUE),
    
    mean_delta_clock_age_2_minus_0 = mean(delta_clock_age_2_minus_0, na.rm = TRUE),
    mean_delta_clock_age_3_minus_2 = mean(delta_clock_age_3_minus_2, na.rm = TRUE),
    mean_delta_clock_age_3_minus_0 = mean(delta_clock_age_3_minus_0, na.rm = TRUE),
    
    .groups = "drop"
  )

fwrite(
  delta_summary_tbl,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_delta_summary_by_group.tsv"),
  sep = "\t"
)

print(delta_summary_tbl)

# -----------------------------
# 10. LMM trend tests by group
# -----------------------------

run_lmm_by_group <- function(df, group_name) {
  n_pid <- n_distinct(df$participant_id)
  n_obs <- nrow(df)
  n_visits <- n_distinct(df$visit_order)
  
  if (n_pid < 5 || n_obs < 10 || n_visits < 2) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "visit_order",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      label = paste0(group_name, ": LMM unavailable")
    ))
  }
  
  fit <- lmer(
    clock_acceleration_years ~ visit_order + (1 | participant_id),
    data = df
  )
  
  tt <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "visit_order")
  
  data.frame(
    analysis_group = group_name,
    n_participants = n_pid,
    n_observations = n_obs,
    term = tt$term[1],
    estimate = tt$estimate[1],
    std.error = tt$std.error[1],
    statistic = tt$statistic[1],
    df = if ("df" %in% names(tt)) tt$df[1] else NA_real_,
    p.value = tt$p.value[1],
    conf.low = tt$conf.low[1],
    conf.high = tt$conf.high[1],
    label = paste0(
      group_name,
      ": LMM \u03B2 = ", fmt_num(tt$estimate[1], 2),
      " acceleration-years/visit (95% CI ",
      fmt_num(tt$conf.low[1], 2), ", ",
      fmt_num(tt$conf.high[1], 2),
      "), P = ", fmt_p(tt$p.value[1])
    )
  )
}

lmm_tbl <- bind_rows(
  run_lmm_by_group(dat_plot %>% filter(analysis_group == "Whole sample"), "Whole sample"),
  run_lmm_by_group(dat_plot %>% filter(analysis_group == "Cases"), "Cases"),
  run_lmm_by_group(dat_plot %>% filter(analysis_group == "Non-cases"), "Non-cases")
) %>%
  mutate(
    analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
  )

fwrite(
  lmm_tbl,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_lmm_visit_order_trend_by_group.tsv"),
  sep = "\t"
)

print(lmm_tbl)

# Optional sample-year LMM by group
run_lmm_year_by_group <- function(df, group_name) {
  df <- df %>% filter(!is.na(sample_year))
  n_pid <- n_distinct(df$participant_id)
  n_obs <- nrow(df)
  
  if (n_pid < 5 || n_obs < 10 || n_distinct(df$sample_year) < 2) {
    return(data.frame(
      analysis_group = group_name,
      n_participants = n_pid,
      n_observations = n_obs,
      term = "sample_year",
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      df = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    ))
  }
  
  fit <- lmer(
    clock_acceleration_years ~ sample_year + (1 | participant_id),
    data = df
  )
  
  tt <- broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "sample_year")
  
  data.frame(
    analysis_group = group_name,
    n_participants = n_pid,
    n_observations = n_obs,
    term = tt$term[1],
    estimate = tt$estimate[1],
    std.error = tt$std.error[1],
    statistic = tt$statistic[1],
    df = if ("df" %in% names(tt)) tt$df[1] else NA_real_,
    p.value = tt$p.value[1],
    conf.low = tt$conf.low[1],
    conf.high = tt$conf.high[1]
  )
}

lmm_year_tbl <- bind_rows(
  run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Whole sample"), "Whole sample"),
  run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Cases"), "Cases"),
  run_lmm_year_by_group(dat_plot %>% filter(analysis_group == "Non-cases"), "Non-cases")
) %>%
  mutate(
    analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases"))
  )

fwrite(
  lmm_year_tbl,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_lmm_sample_year_trend_by_group.tsv"),
  sep = "\t"
)

# -----------------------------
# 11. Annotation positions and labels
# -----------------------------

global_y_max <- max(dat_plot$clock_acceleration_years, na.rm = TRUE)
global_y_min <- min(dat_plot$clock_acceleration_years, na.rm = TRUE)
global_y_range <- global_y_max - global_y_min

ann_tbl <- lmm_tbl %>%
  mutate(
    x = 2,
    y = global_y_max + 0.16 * global_y_range,
    label_short = paste0(
      "LMM \u03B2 = ", fmt_num(estimate, 2),
      " years/visit\nP = ", vapply(p.value, fmt_p, character(1))
    )
  )

delta_ann_tbl <- delta_summary_tbl %>%
  mutate(
    analysis_group = factor(analysis_group, levels = c("Whole sample", "Cases", "Non-cases")),
    x = 2,
    y = global_y_max + 0.08 * global_y_range,
    delta_label = paste0(
      "\u0394 acceleration: 2-0 = ", fmt_num(mean_delta_accel_2_minus_0, 1),
      "; 3-2 = ", fmt_num(mean_delta_accel_3_minus_2, 1),
      "; 3-0 = ", fmt_num(mean_delta_accel_3_minus_0, 1), " years"
    )
  )

caution_ann_tbl <- data.frame(
  analysis_group = factor(c("Whole sample", "Cases", "Non-cases"), levels = c("Whole sample", "Cases", "Non-cases")),
  x = 2,
  y = global_y_max + 0.01 * global_y_range,
  caution_label = "Caution: cross-instance calibration may exaggerate absolute gaps"
)

# -----------------------------
# 12. Distribution plot of acceleration years
# -----------------------------

p_dist_accel <- ggplot(
  dat_plot,
  aes(x = instance_label, y = clock_acceleration_years, fill = instance_label, color = instance_label)
) +
  geom_violin(trim = FALSE, alpha = 0.28, linewidth = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.70, linewidth = 0.5) +
  geom_jitter(width = 0.08, alpha = 0.12, size = 0.65, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2.8, shape = 21, fill = "white", color = "black") +
  geom_text(
    data = ann_tbl,
    aes(x = x, y = y, label = label_short),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  geom_text(
    data = delta_ann_tbl,
    aes(x = x, y = y, label = delta_label),
    inherit.aes = FALSE,
    size = 3.0
  ) +
  geom_text(
    data = caution_ann_tbl,
    aes(x = x, y = y, label = caution_label),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "italic",
    color = "grey30"
  ) +
  facet_wrap(~ analysis_group, nrow = 1) +
  scale_fill_manual(values = instance_palette_fill, guide = "none") +
  scale_color_manual(values = instance_palette, guide = "none") +
  coord_cartesian(
    ylim = c(global_y_min, global_y_max + 0.25 * global_y_range),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Pulmonary L'EPOCH (de)acceleration (years)",
    title = "Pulmonary proteomics L'EPOCH (de)acceleration across UKB instances",
    subtitle = paste0("Common participants across instance 0, 2, and 3: N = ", length(common_ids))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_distribution_by_case_status.pdf"),
  p_dist_accel,
  width = 14.5,
  height = 6.2
)

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_distribution_by_case_status.png"),
  p_dist_accel,
  width = 14.5,
  height = 6.2,
  dpi = 300
)

# -----------------------------
# 13. Mean +/- 95% CI trend plot of acceleration years
# -----------------------------

mean_se_tbl <- dat_plot %>%
  group_by(analysis_group, visit_order, instance_label) %>%
  summarise(
    n = n(),
    n_participants = n_distinct(participant_id),
    mean = mean(clock_acceleration_years, na.rm = TRUE),
    se = sd(clock_acceleration_years, na.rm = TRUE) / sqrt(n()),
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se,
    .groups = "drop"
  )

fwrite(
  mean_se_tbl,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_mean_se_by_instance_and_group.tsv"),
  sep = "\t"
)

p_mean_accel <- ggplot(
  mean_se_tbl,
  aes(x = visit_order, y = mean, group = 1)
) +
  geom_line(linewidth = 1.3, color = "#1C1C1C") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = alpha("#6C8EBF", 0.20), color = NA) +
  geom_point(aes(fill = instance_label), shape = 21, size = 4.2, color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.06, linewidth = 0.6) +
  geom_text(
    data = ann_tbl,
    aes(x = x, y = y, label = label_short),
    inherit.aes = FALSE,
    size = 3.2,
    fontface = "bold"
  ) +
  geom_text(
    data = delta_ann_tbl,
    aes(x = x, y = y, label = delta_label),
    inherit.aes = FALSE,
    size = 3.0
  ) +
  geom_text(
    data = caution_ann_tbl,
    aes(x = x, y = y, label = caution_label),
    inherit.aes = FALSE,
    size = 2.8,
    fontface = "italic",
    color = "grey30"
  ) +
  facet_wrap(~ analysis_group, nrow = 1) +
  scale_fill_manual(values = instance_palette_fill, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("Instance 0\nBaseline", "Instance 2\nImaging", "Instance 3\nRepeat imaging")
  ) +
  coord_cartesian(
    ylim = c(global_y_min, global_y_max + 0.25 * global_y_range),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Mean pulmonary L'EPOCH (de)acceleration (years)",
    title = "Mean pulmonary L'EPOCH (de)acceleration trajectory",
    subtitle = paste0("Common participants across all three instances: N = ", length(common_ids))
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_mean_trend_by_case_status.pdf"),
  p_mean_accel,
  width = 14.5,
  height = 6.2
)

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_mean_trend_by_case_status.png"),
  p_mean_accel,
  width = 14.5,
  height = 6.2,
  dpi = 300
)

# -----------------------------
# 14. Paired spaghetti plot of acceleration years
# -----------------------------

set.seed(2026)

plot_id_tbl <- dat_plot %>%
  distinct(analysis_group, participant_id) %>%
  group_by(analysis_group) %>%
  group_modify(~ {
    n_take <- min(250, nrow(.x))
    .x[sample(seq_len(nrow(.x)), n_take), , drop = FALSE]
  }) %>%
  ungroup()

dat_spaghetti <- dat_plot %>%
  inner_join(plot_id_tbl, by = c("analysis_group", "participant_id"))

mean_traj_tbl <- dat_plot %>%
  group_by(analysis_group, visit_order, instance_label) %>%
  summarise(
    mean = mean(clock_acceleration_years, na.rm = TRUE),
    .groups = "drop"
  )

p_spaghetti_accel <- ggplot(
  dat_spaghetti,
  aes(x = visit_order, y = clock_acceleration_years, group = participant_id)
) +
  geom_line(alpha = 0.10, color = "grey50") +
  geom_point(alpha = 0.15, size = 0.7, color = "grey50") +
  geom_line(
    data = mean_traj_tbl,
    aes(x = visit_order, y = mean),
    inherit.aes = FALSE,
    linewidth = 1.4,
    color = "#1C1C1C"
  ) +
  geom_point(
    data = mean_traj_tbl,
    aes(x = visit_order, y = mean, fill = instance_label),
    inherit.aes = FALSE,
    shape = 21,
    color = "black",
    size = 4
  ) +
  facet_wrap(~ analysis_group, nrow = 1) +
  scale_fill_manual(values = instance_palette_fill, guide = "none") +
  scale_x_continuous(
    breaks = c(0, 1, 2),
    labels = c("Instance 0\nBaseline", "Instance 2\nImaging", "Instance 3\nRepeat imaging")
  ) +
  labs(
    x = NULL,
    y = "Pulmonary L'EPOCH (de)acceleration (years)",
    title = "Within-person pulmonary L'EPOCH (de)acceleration trajectories",
    subtitle = "Thin lines show sampled participants; black line shows mean trajectory"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 10, face = "bold")
  )

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_spaghetti_by_case_status.pdf"),
  p_spaghetti_accel,
  width = 14.5,
  height = 6.2
)

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_spaghetti_by_case_status.png"),
  p_spaghetti_accel,
  width = 14.5,
  height = 6.2,
  dpi = 300
)

# -----------------------------
# 15. QC plot: acceleration vs chronological age
# -----------------------------

p_accel_vs_chrono <- ggplot(
  dat_plot,
  aes(x = chronological_age, y = clock_acceleration_years, color = instance_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_wrap(~ analysis_group, nrow = 1) +
  scale_color_manual(values = instance_palette, name = "Instance") +
  labs(
    x = "Chronological age at proteomics assessment (years)",
    y = "Pulmonary L'EPOCH (de)acceleration (years)",
    title = "Pulmonary L'EPOCH (de)acceleration versus chronological age",
    subtitle = "Dashed line indicates zero acceleration; calibration across instances may affect absolute offsets"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_vs_chronological_age_by_case_status.pdf"),
  p_accel_vs_chrono,
  width = 14.5,
  height = 6.2
)

ggsave(
  file.path(outdir, "pulmonary_LEPOCH_acceleration_years_vs_chronological_age_by_case_status.png"),
  p_accel_vs_chrono,
  width = 14.5,
  height = 6.2,
  dpi = 300
)

# -----------------------------
# 16. Save figure annotation text
# -----------------------------

annotation_tbl <- lmm_tbl %>%
  select(
    analysis_group,
    n_participants,
    n_observations,
    estimate,
    conf.low,
    conf.high,
    p.value,
    label
  )

fwrite(
  annotation_tbl,
  file.path(outdir, "pulmonary_LEPOCH_acceleration_figure_annotation_text_by_group.tsv"),
  sep = "\t"
)

# -----------------------------
# 17. Print key files
# -----------------------------

cat("\nFinished longitudinal pulmonary L'EPOCH (de)acceleration analysis.\n")
cat("Output directory:\n", outdir, "\n\n")
cat("Main outputs:\n")
cat("  pulmonary_LEPOCH_acceleration_common_population_long.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_plot_dataset_with_groups.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_instance_summary_by_group.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_common_population_wide_deltas_by_group.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_delta_summary_by_group.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_lmm_visit_order_trend_by_group.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_lmm_sample_year_trend_by_group.tsv\n")
cat("  pulmonary_LEPOCH_acceleration_years_distribution_by_case_status.pdf/png\n")
cat("  pulmonary_LEPOCH_acceleration_years_mean_trend_by_case_status.pdf/png\n")
cat("  pulmonary_LEPOCH_acceleration_years_spaghetti_by_case_status.pdf/png\n")
cat("  pulmonary_LEPOCH_acceleration_years_vs_chronological_age_by_case_status.pdf/png\n")