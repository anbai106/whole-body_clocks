#!/usr/bin/env Rscript

# Integrated within- and cross-system temporal coupling network
#
# Cross-system arrows:
#   Baseline EPOCH in system A -> annualized delta EPOCH in system B
#   All off-diagonal associations are displayed.
#   Bonferroni-significant edges are solid; all other edges are dashed.
#
# Within-system self-loops:
#   Baseline EPOCH in system A -> follow-up EPOCH in the same system
#   Each system is represented by a closed-arrow self-loop.
#   Self-loop width is proportional to the within-system persistence coefficient.
#
# The negative diagonal coefficients from the delta-change model are intentionally
# excluded because they are susceptible to mathematical coupling and regression
# to the mean.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
  library(grid)
  library(tibble)
})

# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

change_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/metabolomics_temporal_coupling_network/metabolomics_temporal_coupling_coefficients_change.tsv"

followup_file <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/metabolomics_temporal_coupling_network/metabolomics_temporal_coupling_coefficients_followup.tsv"

out_prefix <- "metabolomics_temporal_coupling_integrated_network"

# All off-diagonal cross-system edges are shown.
# Set this above zero only if you later want to suppress extremely small effects.
minimum_abs_beta <- 0

# Fixed node order and coordinates for a reproducible layout.
system_order <- c("Endocrine", "Digestive", "Hepatic", "Immune")

node_positions <- tibble(
  name = system_order,
  x = c(0, -1.40, 1.40, 0),
  y = c(1.25, 0, 0, -1.25)
)

# Plotting ranges.
cross_edge_width_range <- c(0.7, 5.5)
self_loop_width_range  <- c(1.8, 5.5)

node_size <- 15
node_radius_mm <- 8.7

# Self-loop curvature. Larger values make the loop more pronounced.
self_loop_strength <- 1.25

# -----------------------------------------------------------------------------
# 2. Read and validate files
# -----------------------------------------------------------------------------

required_columns <- c(
  "outcome_system",
  "predictor_system",
  "edge_type",
  "beta",
  "p",
  "q_fdr_bh",
  "fdr_significant",
  "bonferroni_significant"
)

read_results <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }

  dat <- read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )

  missing_columns <- setdiff(required_columns, names(dat))

  if (length(missing_columns) > 0) {
    stop(
      "The following required columns are missing from ",
      path,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  dat
}

as_logical_safe <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  tolower(trimws(as.character(x))) %in% c(
    "true", "t", "1", "yes", "y"
  )
}

change_results <- read_results(change_file) %>%
  mutate(
    fdr_significant = as_logical_safe(fdr_significant),
    bonferroni_significant = as_logical_safe(bonferroni_significant)
  )

followup_results <- read_results(followup_file) %>%
  mutate(
    fdr_significant = as_logical_safe(fdr_significant),
    bonferroni_significant = as_logical_safe(bonferroni_significant)
  )

# -----------------------------------------------------------------------------
# 3. Prepare all cross-system delta-change edges
# -----------------------------------------------------------------------------

cross_edges <- change_results %>%
  filter(
    predictor_system != outcome_system,
    edge_type == "cross_system",
    is.finite(beta),
    abs(beta) >= minimum_abs_beta,
    predictor_system %in% system_order,
    outcome_system %in% system_order
  ) %>%
  transmute(
    from = predictor_system,
    to = outcome_system,
    beta = beta,
    p = p,
    q = q_fdr_bh,
    fdr_significant = fdr_significant,
    bonferroni_significant = bonferroni_significant,

    sign = if_else(
      beta >= 0,
      "Positive coupling",
      "Inverse coupling"
    ),

    evidence = if_else(
      bonferroni_significant,
      "Bonferroni-significant",
      "Not Bonferroni-significant"
    ),

    edge_label = sprintf("%.3f", beta),
    edge_strength = abs(beta)
  )

if (nrow(cross_edges) == 0) {
  stop("No valid cross-system edges were found.")
}

expected_cross_edges <- length(system_order) * (length(system_order) - 1)

if (nrow(cross_edges) != expected_cross_edges) {
  warning(
    "Expected ",
    expected_cross_edges,
    " directed off-diagonal cross-system edges, but found ",
    nrow(cross_edges),
    "."
  )
}

# Scale cross-system edge widths. If all effects are equal, use the midpoint.
if (length(unique(cross_edges$edge_strength)) > 1) {
  cross_edges <- cross_edges %>%
    mutate(
      plot_width = rescale(
        edge_strength,
        to = cross_edge_width_range
      )
    )
} else {
  cross_edges <- cross_edges %>%
    mutate(
      plot_width = mean(cross_edge_width_range)
    )
}

# -----------------------------------------------------------------------------
# 4. Prepare within-system persistence self-loops
# -----------------------------------------------------------------------------

within_loops <- followup_results %>%
  filter(
    predictor_system == outcome_system,
    edge_type == "within_system",
    is.finite(beta),
    predictor_system %in% system_order
  ) %>%
  transmute(
    from = predictor_system,
    to = outcome_system,
    persistence_beta = beta,
    persistence_p = p,
    persistence_q = q_fdr_bh,
    persistence_fdr = fdr_significant,
    persistence_bonferroni = bonferroni_significant,

    loop_label = sprintf("%.2f", beta),
    loop_strength = abs(beta),

    evidence = if_else(
      bonferroni_significant,
      "Bonferroni-significant",
      "Not Bonferroni-significant"
    )
  )

if (nrow(within_loops) != length(system_order)) {
  warning(
    "Expected ",
    length(system_order),
    " within-system follow-up coefficients, but found ",
    nrow(within_loops),
    "."
  )
}

if (length(unique(within_loops$loop_strength)) > 1) {
  within_loops <- within_loops %>%
    mutate(
      loop_width = rescale(
        loop_strength,
        to = self_loop_width_range
      )
    )
} else {
  within_loops <- within_loops %>%
    mutate(
      loop_width = mean(self_loop_width_range)
    )
}

# -----------------------------------------------------------------------------
# 5. Prepare node table
# -----------------------------------------------------------------------------

nodes <- node_positions %>%
  left_join(
    within_loops %>%
      select(
        name = from,
        persistence_beta,
        persistence_p,
        persistence_q,
        persistence_fdr,
        persistence_bonferroni,
        loop_width
      ),
    by = "name"
  ) %>%
  mutate(
    node_label = name
  )

if (any(!is.finite(nodes$persistence_beta))) {
  stop(
    "At least one node is missing a valid within-system persistence coefficient."
  )
}

# -----------------------------------------------------------------------------
# 6. Construct graph for cross-system arrows
# -----------------------------------------------------------------------------

graph_object <- graph_from_data_frame(
  d = cross_edges,
  directed = TRUE,
  vertices = nodes
)

layout_matrix <- as.matrix(
  nodes[
    match(V(graph_object)$name, nodes$name),
    c("x", "y")
  ]
)

# -----------------------------------------------------------------------------
# 7. Create self-loop geometry
# -----------------------------------------------------------------------------

# Self-loops are drawn separately with geom_curve so that they appear as
# explicit circular/elliptical closed-arrow loops around each node.
#
# Loop placement is node-specific to avoid collisions:
#   Endocrine loop: above
#   Digestive loop: left
#   Hepatic loop: right
#   Immune loop: below

self_loop_geometry <- nodes %>%
  mutate(
    x_start = case_when(
      name == "Digestive" ~ x - 0.03,
      name == "Hepatic"   ~ x + 0.03,
      TRUE                ~ x - 0.27
    ),

    y_start = case_when(
      name == "Endocrine" ~ y + 0.12,
      name == "Immune"    ~ y - 0.12,
      TRUE                ~ y + 0.26
    ),

    x_end = case_when(
      name == "Digestive" ~ x - 0.03,
      name == "Hepatic"   ~ x + 0.03,
      TRUE                ~ x + 0.27
    ),

    y_end = case_when(
      name == "Endocrine" ~ y + 0.12,
      name == "Immune"    ~ y - 0.12,
      TRUE                ~ y - 0.26
    ),

    loop_curvature = case_when(
      name == "Endocrine" ~ -self_loop_strength,
      name == "Immune"    ~  self_loop_strength,
      name == "Digestive" ~  self_loop_strength,
      name == "Hepatic"   ~ -self_loop_strength,
      TRUE                ~  self_loop_strength
    ),

    label_x = case_when(
      name == "Endocrine" ~ x,
      name == "Immune"    ~ x,
      name == "Digestive" ~ x - 0.60,
      name == "Hepatic"   ~ x + 0.60,
      TRUE                ~ x
    ),

    label_y = case_when(
      name == "Endocrine" ~ y + 0.55,
      name == "Immune"    ~ y - 0.55,
      name == "Digestive" ~ y,
      name == "Hepatic"   ~ y,
      TRUE                ~ y
    ),

    loop_label = sprintf("%.2f", persistence_beta)
  )

# -----------------------------------------------------------------------------
# 8. Plot
# -----------------------------------------------------------------------------

network_plot <- ggraph(
  graph_object,
  layout = "manual",
  x = layout_matrix[, 1],
  y = layout_matrix[, 2]
) +

  # ---------------------------------------------------------------------------
  # Cross-system arrows: all off-diagonal results
  # ---------------------------------------------------------------------------

  geom_edge_fan(
    aes(
      width = plot_width,
      colour = sign,
      linetype = evidence,
      label = edge_label
    ),
    arrow = arrow(
      type = "closed",
      length = unit(4.0, "mm")
    ),
    start_cap = circle(node_radius_mm, "mm"),
    end_cap = circle(node_radius_mm + 0.7, "mm"),
    alpha = 0.78,
    angle_calc = "along",
    label_dodge = unit(2.0, "mm"),
    label_push = unit(1.3, "mm"),
    label_size = 3.0,
    show.legend = TRUE
  ) +

  # ---------------------------------------------------------------------------
  # Within-system closed-arrow self-loops
  # ---------------------------------------------------------------------------

  # geom_curve() requires one scalar curvature per layer.
  # Draw each self-loop in a separate layer so that every node can use
  # its own loop direction without triggering a length > 1 error.
  geom_curve(
    data = self_loop_geometry %>% filter(name == "Endocrine"),
    aes(
      x = x_start, y = y_start,
      xend = x_end, yend = y_end,
      linewidth = loop_width
    ),
    inherit.aes = FALSE,
    curvature = -self_loop_strength,
    colour = "grey20",
    linetype = "solid",
    alpha = 0.95,
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(4.0, "mm")),
    show.legend = FALSE
  ) +

  geom_curve(
    data = self_loop_geometry %>% filter(name == "Digestive"),
    aes(
      x = x_start, y = y_start,
      xend = x_end, yend = y_end,
      linewidth = loop_width
    ),
    inherit.aes = FALSE,
    curvature = self_loop_strength,
    colour = "grey20",
    linetype = "solid",
    alpha = 0.95,
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(4.0, "mm")),
    show.legend = FALSE
  ) +

  geom_curve(
    data = self_loop_geometry %>% filter(name == "Hepatic"),
    aes(
      x = x_start, y = y_start,
      xend = x_end, yend = y_end,
      linewidth = loop_width
    ),
    inherit.aes = FALSE,
    curvature = -self_loop_strength,
    colour = "grey20",
    linetype = "solid",
    alpha = 0.95,
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(4.0, "mm")),
    show.legend = FALSE
  ) +

  geom_curve(
    data = self_loop_geometry %>% filter(name == "Immune"),
    aes(
      x = x_start, y = y_start,
      xend = x_end, yend = y_end,
      linewidth = loop_width
    ),
    inherit.aes = FALSE,
    curvature = self_loop_strength,
    colour = "grey20",
    linetype = "solid",
    alpha = 0.95,
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(4.0, "mm")),
    show.legend = FALSE
  ) +

  # Self-loop coefficient labels.
  geom_text(
    data = self_loop_geometry,
    aes(
      x = label_x,
      y = label_y,
      label = loop_label
    ),
    inherit.aes = FALSE,
    size = 3.3,
    fontface = "bold",
    colour = "grey15"
  ) +

  # ---------------------------------------------------------------------------
  # Nodes
  # ---------------------------------------------------------------------------

  geom_node_point(
    size = node_size,
    shape = 21,
    fill = "grey96",
    colour = "black",
    stroke = 0.8,
    show.legend = FALSE
  ) +

  geom_node_text(
    aes(label = node_label),
    size = 4.0,
    fontface = "bold",
    lineheight = 0.95
  ) +

  # ---------------------------------------------------------------------------
  # Scales and legends
  # ---------------------------------------------------------------------------

  scale_edge_width_identity() +
  scale_linewidth_identity() +

  scale_edge_colour_manual(
    values = c(
      "Positive coupling" = "#B2182B",
      "Inverse coupling" = "#2166AC"
    ),
    breaks = c(
      "Positive coupling",
      "Inverse coupling"
    )
  ) +

  scale_edge_linetype_manual(
    values = c(
      "Bonferroni-significant" = "solid",
      "Not Bonferroni-significant" = "22"
    ),
    breaks = c(
      "Bonferroni-significant",
      "Not Bonferroni-significant"
    )
  ) +

  guides(
    edge_colour = guide_legend(
      title = "Cross-system direction",
      order = 1,
      override.aes = list(
        width = 2.5,
        alpha = 1
      )
    ),

    edge_linetype = guide_legend(
      title = "Cross-system evidence",
      order = 2,
      override.aes = list(
        colour = "grey25",
        width = 2.5,
        alpha = 1
      )
    )
  ) +

  coord_equal(
    xlim = c(-2.25, 2.25),
    ylim = c(-2.05, 2.05),
    clip = "off"
  ) +

  theme_void(base_size = 12) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 10.5,
      hjust = 0.5,
      margin = margin(b = 8)
    ),

    plot.caption = element_text(
      size = 9,
      hjust = 0,
      margin = margin(t = 10)
    ),

    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 30, 20, 30)
  ) +

  labs(
    title = paste0(
      "Within- and cross-system temporal coupling of ",
      "metabolomics mortality EPOCHs"
    ),

    subtitle = paste0(
      "Cross-system arrows show all baseline-to-annualized-change associations; ",
      "closed self-loops show within-system baseline-to-follow-up persistence"
    ),

    caption = paste0(
      "Cross-system arrow direction: baseline predictor system -> subsequent ",
      "annualized change in the outcome system. Cross-system arrow width is ",
      "proportional to |beta|; solid arrows pass Bonferroni correction and ",
      "dashed arrows do not. Red and blue indicate positive and inverse coupling, ",
      "respectively. Closed grey self-loops represent within-system longitudinal ",
      "persistence from the follow-up-level model; loop width is proportional to ",
      "|beta|. Labels report beta coefficients. Negative diagonal baseline-to-change ",
      "coefficients are intentionally excluded."
    )
  )

# -----------------------------------------------------------------------------
# 9. Save outputs
# -----------------------------------------------------------------------------

ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = network_plot,
  width = 10.8,
  height = 9.0,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = network_plot,
  width = 10.8,
  height = 9.0,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = paste0(out_prefix, ".svg"),
  plot = network_plot,
  width = 10.8,
  height = 9.0,
  units = "in",
  bg = "white"
)

# Save the exact plotted data for reproducibility.
write_tsv(
  cross_edges,
  paste0(out_prefix, "_plotted_cross_edges.tsv")
)

write_tsv(
  within_loops,
  paste0(out_prefix, "_plotted_within_system_loops.tsv")
)

write_tsv(
  self_loop_geometry,
  paste0(out_prefix, "_self_loop_geometry.tsv")
)

print(network_plot)

message("Saved:")
message("  ", paste0(out_prefix, ".pdf"))
message("  ", paste0(out_prefix, ".png"))
message("  ", paste0(out_prefix, ".svg"))
message("  ", paste0(out_prefix, "_plotted_cross_edges.tsv"))
message("  ", paste0(out_prefix, "_plotted_within_system_loops.tsv"))
message("  ", paste0(out_prefix, "_self_loop_geometry.tsv"))
