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

x_min <- min(df$Brain_PhenoBAG, na.rm = TRUE)
x_max <- max(df$Brain_PhenoBAG, na.rm = TRUE)
y_min <- min(df$brain_mri_mortality_clock_acceleration_years, na.rm = TRUE)
y_max <- max(df$brain_mri_mortality_clock_acceleration_years, na.rm = TRUE)

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
    x = Brain_PhenoBAG,
    y = brain_mri_mortality_clock_acceleration_years
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
    x = "Brain MRIBAG, years",
    y = "Brain MRI mortality clock, years"
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