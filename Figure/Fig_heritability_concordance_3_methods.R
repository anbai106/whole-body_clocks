# ============================================================
# 5. Pairwise consistency statistics
# ============================================================

format_pvalue <- function(p) {
  if (!is.finite(p)) {
    return("P = NA")
  } else if (p < 2.2e-16) {
    return("P < 2.2e-16")
  } else if (p < 0.001) {
    return(paste0("P = ", formatC(p, format = "e", digits = 2)))
  } else {
    return(paste0("P = ", signif(p, 3)))
  }
}

method_pairs <- combn(method_order, 2, simplify = FALSE)

pairwise_stats <- map_dfr(method_pairs, function(pair) {
  
  m1 <- pair[1]
  m2 <- pair[2]
  
  x_col <- paste0("h2_", m1)
  y_col <- paste0("h2_", m2)
  
  x <- df_wide[[x_col]]
  y <- df_wide[[y_col]]
  
  ok <- is.finite(x) & is.finite(y)
  x_ok <- x[ok]
  y_ok <- y[ok]
  diff <- x_ok - y_ok
  
  pearson <- safe_cor_test(x, y, method = "pearson")
  spearman <- safe_cor_test(x, y, method = "spearman")
  
  paired_t <- if (length(diff) >= 3) {
    suppressWarnings(t.test(x_ok, y_ok, paired = TRUE))
  } else {
    NULL
  }
  
  paired_w <- if (length(diff) >= 3) {
    suppressWarnings(wilcox.test(x_ok, y_ok, paired = TRUE, exact = FALSE))
  } else {
    NULL
  }
  
  tibble(
    method_1 = m1,
    method_2 = m2,
    method_1_label = unname(method_labels[m1]),
    method_2_label = unname(method_labels[m2]),
    comparison = paste(unname(method_labels[m1]), "vs", unname(method_labels[m2])),
    n_clocks = sum(ok),
    
    pearson_r = pearson$estimate,
    pearson_p = pearson$p_value,
    
    spearman_rho = spearman$estimate,
    spearman_p = spearman$p_value,
    
    concordance_correlation = lin_ccc(x, y),
    
    mean_h2_method_1 = mean(x_ok, na.rm = TRUE),
    mean_h2_method_2 = mean(y_ok, na.rm = TRUE),
    
    mean_difference_method1_minus_method2 = mean(diff, na.rm = TRUE),
    sd_difference = sd(diff, na.rm = TRUE),
    median_difference = median(diff, na.rm = TRUE),
    
    mean_absolute_difference = mean(abs(diff), na.rm = TRUE),
    root_mean_square_difference = sqrt(mean(diff^2, na.rm = TRUE)),
    
    paired_t_p = ifelse(is.null(paired_t), NA_real_, paired_t$p.value),
    paired_wilcoxon_p = ifelse(is.null(paired_w), NA_real_, paired_w$p.value),
    
    annotation_label = paste0(
      "Pearson r = ", sprintf("%.2f", pearson$estimate), "\n",
      format_pvalue(pearson$p_value)
    )
  )
})

write_tsv(
  pairwise_stats,
  file.path(output_dir, "h2_3method_pairwise_consistency_statistics.tsv")
)

print(pairwise_stats)


# ============================================================
# 8. Plot 2: pairwise method scatter plots with Pearson annotation
# ============================================================

df_scatter <- map_dfr(method_pairs, function(pair) {
  
  m1 <- pair[1]
  m2 <- pair[2]
  
  x_col <- paste0("h2_", m1)
  y_col <- paste0("h2_", m2)
  
  tibble(
    mortality_clock = df_wide$mortality_clock_chr,
    modality = df_wide$modality,
    organ = df_wide$organ,
    method_1 = m1,
    method_2 = m2,
    method_1_label = unname(method_labels[m1]),
    method_2_label = unname(method_labels[m2]),
    comparison = paste(unname(method_labels[m1]), "vs", unname(method_labels[m2])),
    h2_method_1 = df_wide[[x_col]],
    h2_method_2 = df_wide[[y_col]]
  )
}) %>%
  drop_na(h2_method_1, h2_method_2)

# Annotation position for each facet
scatter_annotation <- df_scatter %>%
  group_by(comparison) %>%
  summarise(
    x_min = min(h2_method_1, na.rm = TRUE),
    x_max = max(h2_method_1, na.rm = TRUE),
    y_min = min(h2_method_2, na.rm = TRUE),
    y_max = max(h2_method_2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_range = x_max - x_min,
    y_range = y_max - y_min,
    x_pos = x_min + 0.05 * ifelse(x_range == 0, 1, x_range),
    y_pos = y_max - 0.08 * ifelse(y_range == 0, 1, y_range)
  ) %>%
  left_join(
    pairwise_stats %>%
      select(comparison, annotation_label),
    by = "comparison"
  )

p_scatter <- ggplot(
  df_scatter,
  aes(x = h2_method_1, y = h2_method_2, shape = modality)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "#444444",
    linewidth = 0.55
  ) +
  geom_point(
    size = 3,
    alpha = 0.9,
    color = "#1F4E79"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.65,
    color = "#DDAA33"
  ) +
  geom_label(
    data = scatter_annotation,
    aes(
      x = x_pos,
      y = y_pos,
      label = annotation_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.4,
    fontface = "bold",
    color = "#1F2A44",
    fill = "#FBF7EF",
    label.size = 0.25,
    label.padding = unit(0.18, "lines")
  ) +
  facet_wrap(
    ~ comparison,
    scales = "free",
    nrow = 1
  ) +
  labs(
    title = "Pairwise agreement of h² estimates across methods",
    x = "Method 1 h²",
    y = "Method 2 h²",
    shape = "Modality"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      color = "#1F2A44"
    ),
    strip.text = element_text(
      face = "bold",
      size = 9
    ),
    legend.position = "top",
    plot.background = element_rect(
      fill = "#FBF7EF",
      color = NA
    )
  )

ggsave(
  file.path(output_dir, "h2_3method_pairwise_scatter.pdf"),
  p_scatter,
  width = 15,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_pairwise_scatter.png"),
  p_scatter,
  width = 15,
  height = 5,
  dpi = 450
)


# ============================================================
# 11. Combined summary figure
# ============================================================

combined_plot <- p_line / p_scatter / p_bland +
  plot_layout(
    heights = c(1.25, 1.10, 1.05)
  ) +
  plot_annotation(
    title = "Cross-method consistency of SNP heritability estimates",
    subtitle = "Comparison of h² across raw-genotype GCTA/HEreg, summary-level SBayesS, and summary-level LDSC",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        color = "#1F2A44"
      ),
      plot.subtitle = element_text(
        size = 12,
        color = "#4A4A4A"
      ),
      plot.background = element_rect(
        fill = "#FBF7EF",
        color = NA
      )
    )
  )

ggsave(
  file.path(output_dir, "h2_3method_consistency_combined_figure.pdf"),
  combined_plot,
  width = 16,
  height = 15,
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "h2_3method_consistency_combined_figure.png"),
  combined_plot,
  width = 16,
  height = 15,
  dpi = 450
)