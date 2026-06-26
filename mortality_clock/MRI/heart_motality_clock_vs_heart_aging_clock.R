# ============================================================
# Scatter plot:
# Heart_MRIBAG vs heart MRI mortality clock acceleration
# Revised elegant version for RStudio
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(scales)
  library(grid)
})

# ============================================================
# 1. Set file paths
# ============================================================

mortality_clock_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/heart_mri_mortality_clock/heart_mri_mortality_clock_predictions.tsv"

# Update if needed
heart_age_file <- "/Users/hao/cubic-home/Reproducibile_paper/SleepAging/data/MomoBAG.tsv"

outdir <- dirname(mortality_clock_file)

out_pdf <- file.path(outdir, "scatter_Heart_MRIBAG_vs_heart_mri_mortality_clock_acceleration_years.pdf")
out_png <- file.path(outdir, "scatter_Heart_MRIBAG_vs_heart_mri_mortality_clock_acceleration_years.png")
out_tsv <- file.path(outdir, "scatter_Heart_MRIBAG_vs_heart_mri_mortality_clock_acceleration_years_data.tsv")
out_stat <- file.path(outdir, "scatter_Heart_MRIBAG_vs_mortality_clock_correlation_stats.tsv")

# ============================================================
# 2. Helper function to read csv/tsv automatically
# ============================================================

read_table_auto <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
  
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    readr::read_tsv(path, show_col_types = FALSE)
  }
}

# ============================================================
# 3. Read data
# ============================================================

mort <- read_table_auto(mortality_clock_file)
bag  <- read_table_auto(heart_age_file)

# Harmonize participant ID column if needed
if (!"participant_id" %in% colnames(bag)) {
  if ("eid" %in% colnames(bag)) {
    bag <- bag %>% rename(participant_id = eid)
  } else if ("id" %in% colnames(bag)) {
    bag <- bag %>% rename(participant_id = id)
  } else {
    stop("Could not find participant_id/eid/id column in the Heart_MRIBAG file.")
  }
}

# Check required columns
if (!"participant_id" %in% colnames(mort)) {
  stop("mortality clock file does not contain participant_id.")
}

if (!"heart_mri_mortality_clock_acceleration_years" %in% colnames(mort)) {
  stop("mortality clock file does not contain heart_mri_mortality_clock_acceleration_years.")
}

# If your BAG file uses another column name, update here
if (!"Heart_MRIBAG" %in% colnames(bag)) {
  message("Available columns in heart age file:")
  print(colnames(bag))
  stop("Heart age file does not contain Heart_MRIBAG. Please rename the column or update the code.")
}

# ============================================================
# 4. Merge and clean
# ============================================================

df <- mort %>%
  select(
    participant_id,
    split,
    age_at_imaging,
    sex,
    heart_mri_mortality_clock_acceleration_years,
    heart_mri_mortality_clock_acceleration_z,
    heart_mri_mortality_risk_score
  ) %>%
  inner_join(
    bag %>%
      select(participant_id, Heart_MRIBAG),
    by = "participant_id"
  ) %>%
  mutate(
    Heart_MRIBAG = as.numeric(Heart_MRIBAG),
    heart_mri_mortality_clock_acceleration_years =
      as.numeric(heart_mri_mortality_clock_acceleration_years),
    split = factor(split, levels = c("train", "validation", "test"))
  ) %>%
  filter(
    is.finite(Heart_MRIBAG),
    is.finite(heart_mri_mortality_clock_acceleration_years)
  )

message("Merged N = ", nrow(df))
message("Participants by split:")
print(table(df$split, useNA = "ifany"))

readr::write_tsv(df, out_tsv)

# ============================================================
# 5. Correlation statistics
# ============================================================

pearson_test <- cor.test(
  df$Heart_MRIBAG,
  df$heart_mri_mortality_clock_acceleration_years,
  method = "pearson"
)

spearman_test <- cor.test(
  df$Heart_MRIBAG,
  df$heart_mri_mortality_clock_acceleration_years,
  method = "spearman",
  exact = FALSE
)

lm_fit <- lm(
  heart_mri_mortality_clock_acceleration_years ~ Heart_MRIBAG,
  data = df
)

lm_summary <- summary(lm_fit)

pearson_r <- unname(pearson_test$estimate)
pearson_p <- pearson_test$p.value
spearman_rho <- unname(spearman_test$estimate)
spearman_p <- spearman_test$p.value
lm_beta <- coef(lm_fit)[["Heart_MRIBAG"]]
lm_r2 <- lm_summary$r.squared

cor_tbl <- tibble(
  n = nrow(df),
  pearson_r = pearson_r,
  pearson_p = pearson_p,
  spearman_rho = spearman_rho,
  spearman_p = spearman_p,
  lm_beta_years_per_Heart_MRIBAG_year = lm_beta,
  lm_r2 = lm_r2
)

readr::write_tsv(cor_tbl, out_stat)
print(cor_tbl)

# ============================================================
# 6. Plot settings
# ============================================================

# Colors similar in spirit to the screenshot
point_color <- "#4B006E"   # deep purple
line_color  <- "#000000"   # black
title_color <- "#3E82A8"   # blue-ish for top annotation/title

theme_elegant <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", color = "#222222"),
      axis.text = element_text(color = "#333333"),
      axis.line = element_line(linewidth = 0.8, color = "#222222"),
      axis.ticks = element_line(linewidth = 0.7, color = "#222222"),
      plot.title = element_text(
        face = "plain",
        color = title_color,
        size = base_size + 2,
        hjust = 0.02,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        color = title_color,
        size = base_size + 1,
        hjust = 0.02,
        margin = margin(b = 8)
      ),
      plot.margin = margin(14, 14, 10, 10),
      legend.position = "none"
    )
}

# ============================================================
# 7. Annotation text
# ============================================================

stat_text <- paste0(
  "R = ", round(pearson_r, 2),
  "; P = ", format.pval(pearson_p, digits = 3, eps = 1e-300),
  "; R\u00b2 = ", round(lm_r2, 3)
)

# ============================================================
# 8. Create extra top space so annotation does not overlap points
# ============================================================

x_min <- min(df$Heart_MRIBAG, na.rm = TRUE)
x_max <- max(df$Heart_MRIBAG, na.rm = TRUE)
y_min <- min(df$heart_mri_mortality_clock_acceleration_years, na.rm = TRUE)
y_max <- max(df$heart_mri_mortality_clock_acceleration_years, na.rm = TRUE)

x_range <- x_max - x_min
y_range <- y_max - y_min

# Add upper y-axis space for top annotation
y_upper_plot <- y_max + 0.20 * y_range
y_lower_plot <- y_min - 0.04 * y_range

annot_x <- x_min + 0.04 * x_range
annot_y <- y_max + 0.14 * y_range

# ============================================================
# 9. Elegant plot
# ============================================================

point_color <- "#4B006E"   # deep purple
line_color  <- "#000000"   # black
stat_color  <- "#3E82A8"   # blue annotation

theme_elegant <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "#17202A"),
      axis.title = element_text(face = "bold", color = "#222222"),
      axis.text = element_text(color = "#333333"),
      axis.line = element_line(linewidth = 0.8, color = "#222222"),
      axis.ticks = element_line(linewidth = 0.7, color = "#222222"),
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(12, 16, 10, 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

p <- ggplot(
  df,
  aes(
    x = Heart_MRIBAG,
    y = heart_mri_mortality_clock_acceleration_years
  )
) +
  geom_point(
    color = point_color,
    alpha = 0.55,
    size = 2.3,
    shape = 16,
    stroke = 0
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 1.8,
    color = line_color
  ) +
  annotate(
    "text",
    x = annot_x,
    y = annot_y,
    label = stat_text,
    hjust = 0,
    vjust = 1,
    size = 5.2,
    color = stat_color,
    fontface = "italic"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(y_lower_plot, y_upper_plot),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Heart MRIBAG, years",
    y = "Heart MRI mortality clock, years"
  ) +
  theme_elegant(base_size = 15)

print(p)

# ============================================================
# 10. Save
# ============================================================

ggsave(
  filename = out_pdf,
  plot = p,
  width = 5.4,
  height = 4.6,
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = p,
  width = 5.4,
  height = 4.6,
  dpi = 400
)

message("Saved:")
message("  ", out_pdf)
message("  ", out_png)
message("  ", out_tsv)
message("  ", out_stat)