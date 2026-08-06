# ============================================================
# Brain proteomics mortality EPOCH: four interacting biological systems
#
# Purpose
#   Create a publication-ready two-panel figure summarizing how the major
#   FUMA-supported loci/genes converge on four interacting biological systems:
#     1) Neurodegenerative vulnerability
#     2) Vascular and thrombotic risk
#     3) Inflammaging and immune regulation
#     4) Metabolic and nutrient-sensing pathways
#
# Figure structure
#   Panel A: locus/gene-to-system network centered on brain proteomics mortality EPOCH
#   Panel B: evidence matrix showing which representative loci support each system
#
# Notes
#   - The mapping below is intentionally conservative and editable.
#   - A locus may connect to more than one system.
#   - This is a biological synthesis figure, not a formal enrichment analysis.
#   - Replace or expand the curated mapping after locus-specific colocalization,
#     fine-mapping, or experimental validation.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(patchwork)
  library(scales)
})

# ----------------------------
# 1. Input / output controls
# ----------------------------
out_dir <- "/Users/hao/Downloads"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_prefix <- file.path(
  out_dir,
  "brain_proteomics_mortality_EPOCH_four_interacting_systems"
)

# Font: use Times New Roman if available; otherwise R will fall back gracefully.
base_family <- "Times New Roman"

# Main figure controls
figure_width  <- 16
figure_height <- 10

# ----------------------------
# 2. Curated FUMA-supported locus table
# ----------------------------
# Evidence strength is a visual summary only:
#   3 = highly prominent / multi-layer support
#   2 = strong support
#   1 = supportive / hypothesis-generating
#
# Edit the genes and system assignments as your locus interpretation evolves.

locus_tbl <- tribble(
  ~locus_id,              ~lead_variant, ~chr, ~lead_p,      ~display_genes,                 ~locus_short,          ~evidence_note,
  "ICAM-TYK2 region",     "rs281439",    19,  5.74646e-148, "ICAM1/3/4/5; TYK2; DNMT1",    "ICAM-TYK2",          "Immune adhesion, cytokine signaling, endothelial inflammation",
  "ABO region",           "rs532436",    9,   1.31279e-32,  "ABO; ADAMTS13-linked eQTLs",  "ABO",                "Coagulation, thrombosis, vascular and blood-trait biology",
  "APOE cluster",         "rs429358",    19,  2.28e-26,     "APOE; TOMM40; APOC1/2/4",     "APOE-TOMM40",        "Neurodegeneration, lipid transport, atherosclerotic risk",
  "RPTOR-NPTX1 region",   "rs79155262",  17,  4.68062e-21,  "RPTOR; NPTX1; RNF213; SGSH",  "RPTOR-NPTX1",        "mTORC1 nutrient sensing and neuronal maintenance",
  "CNDP1-CNDP2 region",   "rs4891560",   18,  9.88761e-21,  "CNDP1; CNDP2",                "CNDP1/2",            "Carnosine metabolism, oxidative stress and metabolic resilience",
  "SORL1 region",         "rs2370794",   11,  3.57362e-16,  "SORL1; HSPA8; CRTAM; SC5D",   "SORL1",              "Endosomal trafficking, amyloid biology and immune regulation",
  "MLXIPL-BCL7B region",  "rs35107030",  7,   1.64217e-08,  "MLXIPL; BCL7B; NSUN5; BAZ1B", "MLXIPL-BCL7B",       "Carbohydrate-responsive lipogenesis and multi-tissue regulation",
  "C1QB-associated signal","rs6852182",   4,   2.39172e-10,  "C1QB; DLGAP3; RIMS3; IPO13",  "C1QB-related",       "Complement-linked immune signaling with neuronal connections"
) %>%
  mutate(
    minus_log10_p = -log10(lead_p),
    locus_label = paste0(locus_short, "\n", lead_variant),
    locus_label_matrix = paste0(locus_short, " (", lead_variant, ")")
  )

# ----------------------------
# 3. Four-system definitions
# ----------------------------
system_tbl <- tribble(
  ~system_id, ~system_label,                              ~system_short,               ~system_summary,
  "neuro",    "Neurodegenerative vulnerability",        "Neurodegeneration",        "Amyloid biology, neuronal maintenance, endosomal and protein trafficking",
  "vascular", "Vascular and thrombotic risk",            "Vascular / thrombosis",    "Endothelial dysfunction, coagulation, atherosclerosis and cerebrovascular injury",
  "immune",   "Inflammaging and immune regulation",      "Inflammaging / immunity",  "Chronic inflammation, complement, leukocyte adhesion and cytokine signaling",
  "metabolic","Metabolic and nutrient-sensing pathways", "Metabolic / nutrient sensing", "Lipid transport, mTOR signaling, oxidative stress and metabolic resilience"
)

# ----------------------------
# 4. Curated locus-to-system mapping
# ----------------------------
# Weight controls edge width and the matrix bubble size.
# Keep weights interpretable rather than overly precise.

mapping_tbl <- tribble(
  ~locus_id,               ~system_id,   ~weight, ~support_level,
  "ICAM-TYK2 region",      "vascular",   3,       "High",
  "ICAM-TYK2 region",      "immune",     3,       "High",

  "ABO region",            "vascular",   3,       "High",
  "ABO region",            "immune",     1,       "Supportive",

  "APOE cluster",          "neuro",      3,       "High",
  "APOE cluster",          "vascular",   3,       "High",
  "APOE cluster",          "metabolic",  3,       "High",

  "RPTOR-NPTX1 region",    "neuro",      2,       "Strong",
  "RPTOR-NPTX1 region",    "metabolic",  3,       "High",

  "CNDP1-CNDP2 region",    "vascular",   1,       "Supportive",
  "CNDP1-CNDP2 region",    "metabolic",  2,       "Strong",

  "SORL1 region",          "neuro",      3,       "High",
  "SORL1 region",          "immune",     1,       "Supportive",

  "MLXIPL-BCL7B region",   "metabolic",  3,       "High",
  "MLXIPL-BCL7B region",   "immune",     1,       "Supportive",

  "C1QB-associated signal","immune",     3,       "High",
  "C1QB-associated signal","neuro",      2,       "Strong"
)

# Optional system-to-system interactions for the conceptual network.
system_interactions <- tribble(
  ~from,       ~to,         ~interaction_label,
  "neuro",    "vascular",  "Neurovascular injury",
  "vascular", "immune",    "Endothelial inflammation",
  "immune",   "metabolic", "Immunometabolic signaling",
  "metabolic","neuro",     "Nutrient sensing and neuronal resilience",
  "vascular", "metabolic", "Lipid and cardiometabolic risk",
  "immune",   "neuro",     "Neuroinflammation"
)

# ----------------------------
# 5. Build network data
# ----------------------------
center_node <- tibble(
  name = "Brain proteomics mortality EPOCH",
  node_type = "center",
  label = "Brain proteomics\nmortality EPOCH",
  system_id = NA_character_,
  minus_log10_p = NA_real_
)

system_nodes <- system_tbl %>%
  transmute(
    name = system_id,
    node_type = "system",
    label = system_label,
    system_id = system_id,
    minus_log10_p = NA_real_
  )

locus_nodes <- locus_tbl %>%
  transmute(
    name = locus_id,
    node_type = "locus",
    label = locus_label,
    system_id = NA_character_,
    minus_log10_p = minus_log10_p
  )

nodes <- bind_rows(center_node, system_nodes, locus_nodes)

center_edges <- system_tbl %>%
  transmute(
    from = "Brain proteomics mortality EPOCH",
    to = system_id,
    edge_type = "center_to_system",
    weight = 2.8,
    label = NA_character_
  )

locus_edges <- mapping_tbl %>%
  transmute(
    from = system_id,
    to = locus_id,
    edge_type = "system_to_locus",
    weight = as.numeric(weight),
    label = NA_character_
  )

interaction_edges <- system_interactions %>%
  transmute(
    from,
    to,
    edge_type = "system_interaction",
    weight = 0.8,
    label = interaction_label
  )

edges <- bind_rows(center_edges, locus_edges, interaction_edges)

g <- graph_from_data_frame(
  d = edges,
  directed = FALSE,
  vertices = nodes
)

# ----------------------------
# 6. Fixed radial layout
# ----------------------------
# Fixed coordinates produce a stable, publication-ready arrangement.
# The four systems occupy the cardinal directions; loci are placed around them.

layout_tbl <- tribble(
  ~name,                               ~x,    ~y,
  "Brain proteomics mortality EPOCH",  0.0,   0.0,
  "neuro",                             0.0,   3.0,
  "vascular",                          3.8,   0.0,
  "immune",                            0.0,  -3.0,
  "metabolic",                        -3.8,   0.0,

  "APOE cluster",                      1.8,   4.7,
  "SORL1 region",                     -1.2,   5.2,
  "RPTOR-NPTX1 region",               -2.8,   4.0,
  "C1QB-associated signal",            2.6,  -4.2,
  "ICAM-TYK2 region",                  4.9,  -2.5,
  "ABO region",                        6.0,   0.8,
  "CNDP1-CNDP2 region",               -5.4,  -1.8,
  "MLXIPL-BCL7B region",              -5.8,   1.6
)

layout_mat <- layout_tbl %>%
  right_join(tibble(name = V(g)$name), by = "name") %>%
  arrange(match(name, V(g)$name)) %>%
  select(x, y) %>%
  as.matrix()

# ----------------------------
# 7. Visual settings
# ----------------------------
system_colors <- c(
  "neuro" = "#4C78A8",
  "vascular" = "#E45756",
  "immune" = "#72B7B2",
  "metabolic" = "#F2CF5B"
)

node_type_shapes <- c(
  "center" = 21,
  "system" = 22,
  "locus" = 21
)

# Assign each locus a display color based on its strongest mapped system.
locus_primary_system <- mapping_tbl %>%
  arrange(locus_id, desc(weight), system_id) %>%
  group_by(locus_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(locus_id, primary_system = system_id)

nodes_for_plot <- nodes %>%
  left_join(locus_primary_system, by = c("name" = "locus_id")) %>%
  mutate(
    fill_group = case_when(
      node_type == "center" ~ "center",
      node_type == "system" ~ system_id,
      node_type == "locus" ~ primary_system,
      TRUE ~ "center"
    ),
    node_size = case_when(
      node_type == "center" ~ 18,
      node_type == "system" ~ 13,
      node_type == "locus" ~ rescale(minus_log10_p, to = c(5.5, 9.5), from = range(locus_tbl$minus_log10_p)),
      TRUE ~ 5
    ),
    label_size = case_when(
      node_type == "center" ~ 5.0,
      node_type == "system" ~ 3.6,
      TRUE ~ 3.0
    ),
    fontface = case_when(
      node_type %in% c("center", "system") ~ "bold",
      TRUE ~ "plain"
    )
  )

V(g)$fill_group <- nodes_for_plot$fill_group[match(V(g)$name, nodes_for_plot$name)]
V(g)$node_type <- nodes_for_plot$node_type[match(V(g)$name, nodes_for_plot$name)]
V(g)$node_size <- nodes_for_plot$node_size[match(V(g)$name, nodes_for_plot$name)]
V(g)$label_size <- nodes_for_plot$label_size[match(V(g)$name, nodes_for_plot$name)]
V(g)$fontface <- nodes_for_plot$fontface[match(V(g)$name, nodes_for_plot$name)]
V(g)$label <- nodes_for_plot$label[match(V(g)$name, nodes_for_plot$name)]

E(g)$edge_type <- edges$edge_type
E(g)$weight <- edges$weight

fill_values <- c(
  "center" = "white",
  system_colors
)

# ----------------------------
# 8. Panel A: biological systems network
# ----------------------------
p_network <- ggraph(g, layout = "manual", x = layout_mat[, 1], y = layout_mat[, 2]) +

  # Conceptual interactions among systems.
  geom_edge_link(
    aes(
      filter = edge_type == "system_interaction",
      width = weight
    ),
    color = "grey60",
    linetype = "dashed",
    alpha = 0.55,
    show.legend = FALSE
  ) +

  # EPOCH-to-system links.
  geom_edge_link(
    aes(
      filter = edge_type == "center_to_system",
      width = weight
    ),
    color = "grey25",
    alpha = 0.75,
    show.legend = FALSE
  ) +

  # System-to-locus links, colored by system of origin.
  geom_edge_link(
    aes(
      filter = edge_type == "system_to_locus",
      width = weight,
      color = node1.fill_group
    ),
    alpha = 0.78,
    show.legend = FALSE
  ) +

  scale_edge_width(range = c(0.5, 2.0)) +
  scale_edge_color_manual(values = fill_values) +

  geom_node_point(
    aes(
      shape = node_type,
      size = node_size,
      fill = fill_group
    ),
    color = "grey15",
    stroke = 0.8
  ) +

  geom_node_text(
    aes(
      label = label,
      size = label_size,
      fontface = fontface
    ),
    family = base_family,
    lineheight = 0.95,
    repel = FALSE,
    color = "grey10"
  ) +

  scale_shape_manual(values = node_type_shapes, guide = "none") +
  scale_size_identity() +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_color_manual(values = fill_values, guide = "none") +
  coord_equal(clip = "off") +
  labs(
    tag = "A",
    title = "Brain proteomics mortality EPOCH converges on four interacting biological systems",
    subtitle = paste0(
      "Representative FUMA-supported loci are linked to one or more systems. ",
      "Locus-node size reflects -log10(P); edge width reflects curated support strength."
    )
  ) +
  theme_void(base_family = base_family) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "grey30", hjust = 0),
    plot.tag = element_text(face = "bold", size = 20),
    plot.margin = margin(10, 20, 10, 20)
  )

# ----------------------------
# 9. Panel B: locus-by-system evidence matrix
# ----------------------------
matrix_tbl <- mapping_tbl %>%
  left_join(
    locus_tbl %>%
      select(locus_id, locus_label_matrix, lead_p, minus_log10_p, display_genes),
    by = "locus_id"
  ) %>%
  left_join(
    system_tbl %>% select(system_id, system_short),
    by = "system_id"
  ) %>%
  mutate(
    locus_label_matrix = fct_reorder(locus_label_matrix, minus_log10_p),
    system_short = factor(
      system_short,
      levels = c(
        "Neurodegeneration",
        "Vascular / thrombosis",
        "Inflammaging / immunity",
        "Metabolic / nutrient sensing"
      )
    ),
    support_level = factor(
      support_level,
      levels = c("Supportive", "Strong", "High")
    )
  )

p_matrix <- ggplot(
  matrix_tbl,
  aes(x = system_short, y = locus_label_matrix)
) +
  geom_tile(
    fill = "white",
    color = "grey88",
    linewidth = 0.35
  ) +
  geom_point(
    aes(
      size = weight,
      fill = system_id,
      shape = support_level
    ),
    color = "grey15",
    stroke = 0.45,
    alpha = 0.96
  ) +
  scale_size_continuous(
    range = c(3.5, 9),
    breaks = c(1, 2, 3),
    labels = c("Supportive", "Strong", "High"),
    name = "Curated support"
  ) +
  scale_shape_manual(
    values = c("Supportive" = 21, "Strong" = 22, "High" = 24),
    name = "Support level"
  ) +
  scale_fill_manual(
    values = system_colors,
    guide = "none"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    tag = "B",
    title = "Representative loci show cross-system pleiotropy",
    subtitle = "Several loci connect the brain-derived mortality phenotype to more than one biological system.",
    x = NULL,
    y = NULL,
    caption = paste0(
      "This synthesis is based on FUMA positional mapping, eQTL mapping, chromatin interactions, ",
      "and prior GWAS-trait overlap. It is hypothesis-generating and does not establish causal genes."
    )
  ) +
  theme_minimal(base_family = base_family, base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      size = 10,
      angle = 20,
      hjust = 0
    ),
    axis.text.y = element_text(size = 9.5, color = "grey15"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey30"),
    plot.caption = element_text(size = 8.5, color = "grey35", hjust = 0),
    plot.tag = element_text(face = "bold", size = 20),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ----------------------------
# 10. Combine and save
# ----------------------------
combined_plot <- p_network / p_matrix +
  plot_layout(heights = c(1.35, 1)) +
  plot_annotation(
    title = "Genetic architecture of the brain proteomics mortality EPOCH",
    subtitle = paste0(
      "FUMA-supported loci converge on neurodegenerative, vascular-thrombotic, immune-inflammatory, ",
      "and metabolic-nutrient-sensing mechanisms relevant to lifespan and mortality."
    ),
    theme = theme(
      plot.title = element_text(
        family = base_family,
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        family = base_family,
        size = 11.5,
        hjust = 0.5,
        color = "grey25"
      ),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

print(combined_plot)

# Save data underlying the synthesis figure.
readr::write_tsv(locus_tbl, paste0(out_prefix, "_locus_table.tsv"))
readr::write_tsv(system_tbl, paste0(out_prefix, "_system_table.tsv"))
readr::write_tsv(mapping_tbl, paste0(out_prefix, "_locus_to_system_mapping.tsv"))
readr::write_tsv(system_interactions, paste0(out_prefix, "_system_interactions.tsv"))

# PDF
ggsave(
  filename = paste0(out_prefix, ".pdf"),
  plot = combined_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# PNG
ggsave(
  filename = paste0(out_prefix, ".png"),
  plot = combined_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = 500,
  bg = "white"
)

# SVG, when svglite is installed.
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    filename = paste0(out_prefix, ".svg"),
    plot = combined_plot,
    width = figure_width,
    height = figure_height,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
} else {
  message("Package 'svglite' is not installed; SVG output was skipped.")
}

message("Saved figure and supporting tables with prefix:")
message("  ", out_prefix)
