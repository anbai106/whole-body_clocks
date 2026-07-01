library(tidyverse)
library(patchwork)

# ============================================================
# 1. Input and output
# ============================================================

input_tsv <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result/GCTB_SBayesS_parameters_mortality_clocks.tsv"

output_dir <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/Result"

out_pdf <- file.path(output_dir, "GCTB_SBayesS_parameters_mortality_clocks_3panel.pdf")
out_png <- file.path(output_dir, "GCTB_SBayesS_parameters_mortality_clocks_3panel.png")

df <- read_tsv(input_tsv, show_col_types = FALSE)


# ============================================================
# 2. Define desired clock order
# ============================================================

clock_order <- c(
  # 7 MRI mortality clocks
  "adipose_mri",
  "brain_mri",
  "heart_mri",
  "kidney_mri",
  "liver_mri",
  "pancreas_mri",
  "spleen_mri",

  # 11 proteomics mortality clocks
  "Brain_proteomics",
  "Endocrine_proteomics",
  "Eye_proteomics",
  "Heart_proteomics",
  "Hepatic_proteomics",
  "Immune_proteomics",
  "Pulmonary_proteomics",
  "Renal_proteomics",
  "Reproductive_female_proteomics",
  "Reproductive_male_proteomics",
  "Skin_proteomics",

  # 4 metabolomics mortality clocks
  "Digestive_metabolomics",
  "Endocrine_metabolomics",
  "Hepatic_metabolomics",
  "Immune_metabolomics"
)


# ============================================================
# 3. Organ parsing and display labels
# ============================================================

df_plot_base <- df %>%
  filter(status == "ok") %>%
  mutate(
    mortality_clock = factor(mortality_clock, levels = clock_order),
    modality = factor(modality, levels = c("MRI", "Proteomics", "Metabolomics")),

    organ_raw = mortality_clock %>%
      as.character() %>%
      str_remove("_mortality_clock$") %>%
      str_remove("_mri$") %>%
      str_remove("_proteomics$") %>%
      str_remove("_metabolomics$"),

    organ_clean = case_when(
      str_to_lower(organ_raw) %in% c("liver", "hepatic") ~ "Liver/Hepatic",
      str_to_lower(organ_raw) %in% c("kidney", "renal") ~ "Kidney/Renal",
      str_to_lower(organ_raw) %in% c("heart", "cardiac") ~ "Heart",
      str_to_lower(organ_raw) %in% c("lung", "pulmonary") ~ "Pulmonary",
      str_to_lower(organ_raw) == "reproductive_female" ~ "Reproductive female",
      str_to_lower(organ_raw) == "reproductive_male" ~ "Reproductive male",
      TRUE ~ str_to_title(str_replace_all(organ_raw, "_", " "))
    ),

    clock_label = case_when(
      modality == "MRI" ~ str_replace_all(as.character(mortality_clock), "_mri", ""),
      modality == "Proteomics" ~ str_replace_all(as.character(mortality_clock), "_proteomics", ""),
      modality == "Metabolomics" ~ str_replace_all(as.character(mortality_clock), "_metabolomics", ""),
      TRUE ~ as.character(mortality_clock)
    ),

    clock_label = str_replace_all(clock_label, "_", " "),
    clock_label = str_to_title(clock_label),
    clock_label = factor(clock_label, levels = unique(clock_label[order(mortality_clock)]))
  ) %>%
  arrange(mortality_clock)


# ============================================================
# 4. Van Gogh-inspired organ palette
#    Liver and hepatic share the same color through organ_clean.
# ============================================================

vangogh_organ_colors <- c(
  "Adipose" = "#DDAA33",              # sunflower ochre
  "Brain" = "#1F4E79",                # starry-night blue
  "Heart" = "#B55239",                # burnt sienna
  "Kidney/Renal" = "#3F7F6B",         # cypress green
  "Liver/Hepatic" = "#C69214",        # golden wheat
  "Pancreas" = "#6D5BA6",             # iris violet
  "Spleen" = "#2C5C8A",               # deep cobalt
  "Endocrine" = "#E1B94F",            # sunflower yellow
  "Eye" = "#4A6FA5",                  # blue iris
  "Immune" = "#5E8C61",               # olive green
  "Pulmonary" = "#4C90A8",            # sky teal
  "Reproductive female" = "#B47BA5",  # muted rose
  "Reproductive male" = "#7B6FB0",    # violet blue
  "Skin" = "#C77C3A",                 # warm orange
  "Digestive" = "#8A6F3D"             # earthy ochre brown
)

# Add fallback colors if any organ names are not covered.
missing_organs <- setdiff(unique(df_plot_base$organ_clean), names(vangogh_organ_colors))
if (length(missing_organs) > 0) {
  fallback_cols <- scales::hue_pal()(length(missing_organs))
  names(fallback_cols) <- missing_organs
  vangogh_organ_colors <- c(vangogh_organ_colors, fallback_cols)
}


# ============================================================
# 5. Long-format table for three parameters
# ============================================================

df_long <- bind_rows(
  df_plot_base %>%
    transmute(
      mortality_clock,
      clock_label,
      modality,
      organ_clean,
      parameter = "SNP heritability (h²)",
      estimate = h2_mean,
      se = h2_se
    ),

  df_plot_base %>%
    transmute(
      mortality_clock,
      clock_label,
      modality,
      organ_clean,
      parameter = "Selection signature (S)",
      estimate = S,
      se = S_se
    ),

  df_plot_base %>%
    transmute(
      mortality_clock,
      clock_label,
      modality,
      organ_clean,
      parameter = "Polygenicity (Pi)",
      estimate = Pi,
      se = Pi_se
    )
) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c(
        "SNP heritability (h²)",
        "Selection signature (S)",
        "Polygenicity (Pi)"
      )
    ),
    x_id = as.numeric(mortality_clock)
  )


# ============================================================
# 6. Modality block annotation
# ============================================================

modality_blocks <- df_plot_base %>%
  distinct(mortality_clock, modality) %>%
  mutate(x_id = as.numeric(mortality_clock)) %>%
  group_by(modality) %>%
  summarise(
    xmin = min(x_id) - 0.5,
    xmax = max(x_id) + 0.5,
    xmid = mean(range(x_id)),
    .groups = "drop"
  )

separator_df <- tibble(
  xintercept = c(7.5, 18.5)
)


# ============================================================
# 7. Plot function
# ============================================================

make_parameter_panel <- function(dat, parameter_name, y_lab, add_x_text = FALSE) {

  dat_sub <- dat %>%
    filter(parameter == parameter_name)

  y_range <- range(
    c(dat_sub$estimate - dat_sub$se, dat_sub$estimate + dat_sub$se),
    na.rm = TRUE
  )

  y_pad <- diff(y_range) * 0.18
  if (!is.finite(y_pad) || y_pad == 0) {
    y_pad <- 0.05
  }

  y_top <- y_range[2] + y_pad
  y_bottom <- y_range[1] - y_pad

  p <- ggplot(dat_sub, aes(x = mortality_clock, y = estimate, fill = organ_clean)) +

    # Modality background bands
    geom_rect(
      data = modality_blocks,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "#F6F1E6",
      alpha = 0.55
    ) +

    # Zero line
    geom_hline(
      yintercept = 0,
      linewidth = 0.45,
      color = "#36454F",
      linetype = "dashed"
    ) +

    # Modality separators
    geom_vline(
      data = separator_df,
      aes(xintercept = xintercept),
      inherit.aes = FALSE,
      linewidth = 0.45,
      color = "#3B3B3B",
      linetype = "dotted"
    ) +

    # Bars
    geom_col(
      width = 0.72,
      color = "#1D1D1D",
      linewidth = 0.22,
      alpha = 0.94
    ) +

    # Standard errors
    geom_errorbar(
      aes(ymin = estimate - se, ymax = estimate + se),
      width = 0.22,
      linewidth = 0.48,
      color = "#1D1D1D"
    ) +

    # Modality labels
    geom_text(
      data = modality_blocks,
      aes(x = xmid, y = y_top, label = modality),
      inherit.aes = FALSE,
      size = 4.0,
      fontface = "bold",
      color = "#283747"
    ) +

    scale_fill_manual(values = vangogh_organ_colors, name = "Organ/system") +

    scale_x_discrete(
      limits = clock_order,
      labels = df_plot_base$clock_label[match(clock_order, as.character(df_plot_base$mortality_clock))]
    ) +

    coord_cartesian(ylim = c(y_bottom, y_top), clip = "off") +

    labs(
      x = NULL,
      y = y_lab,
      title = parameter_name
    ) +

    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        color = "#1F2A44",
        margin = margin(b = 8)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 12,
        color = "#1F2A44",
        margin = margin(r = 8)
      ),
      axis.text.y = element_text(
        size = 10,
        color = "#2F2F2F"
      ),
      axis.line = element_line(
        linewidth = 0.35,
        color = "#2F2F2F"
      ),
      axis.ticks = element_line(
        linewidth = 0.35,
        color = "#2F2F2F"
      ),
      panel.grid.major.y = element_line(
        linewidth = 0.25,
        color = "#E3D7C4"
      ),
      panel.grid.minor.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 9),
      plot.margin = margin(t = 12, r = 15, b = 8, l = 8)
    )

  if (add_x_text) {
    p <- p +
      theme(
        axis.text.x = element_text(
          angle = 55,
          hjust = 1,
          vjust = 1,
          size = 9,
          color = "#2F2F2F"
        )
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  return(p)
}


# ============================================================
# 8. Build three horizontal panels
# ============================================================

p_h2 <- make_parameter_panel(
  dat = df_long,
  parameter_name = "SNP heritability (h²)",
  y_lab = "h² estimate ± SE",
  add_x_text = FALSE
)

p_s <- make_parameter_panel(
  dat = df_long,
  parameter_name = "Selection signature (S)",
  y_lab = "S estimate ± SE",
  add_x_text = FALSE
)

p_pi <- make_parameter_panel(
  dat = df_long,
  parameter_name = "Polygenicity (Pi)",
  y_lab = "Pi estimate ± SE",
  add_x_text = TRUE
)

final_plot <- (p_h2 / p_s / p_pi) +
  plot_layout(heights = c(1, 1, 1), guides = "collect") &
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "#FBF7EF", color = NA)
  )

final_plot <- final_plot +
  plot_annotation(
    title = "SBayesS genetic architecture of 22 mortality L’EPOCH clocks",
    subtitle = "SNP heritability, selection signature, and polygenicity across MRI, proteomics, and metabolomics mortality clocks",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        color = "#1F2A44",
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 12.5,
        color = "#4A4A4A",
        margin = margin(b = 12)
      )
    )
  )


# ============================================================
# 9. Save
# ============================================================

ggsave(
  filename = out_pdf,
  plot = final_plot,
  width = 16,
  height = 10,
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = final_plot,
  width = 16,
  height = 10,
  dpi = 450
)

print(final_plot)