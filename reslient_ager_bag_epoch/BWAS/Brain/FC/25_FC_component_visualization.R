# ============================================================
# Network plot (UKB group-ICA d25 “good nodes” connectome-style)
# - Reads BWAS_result_Brain_FC.tsv (BAG x FC-edge results)
# - Keeps only P_value < 0.05/720/9
# - Maps IDP f_1..f_M to node pairs (standard lower-tri ordering)
# - Node size ~ number of significant BAG–IDP associations incident to node
# - Edge linewidth ~ number of significant BAGs for that FC edge
# - Labels top edges with the most significant BAG + f_*
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# -----------------------
# User inputs
# -----------------------
in_tsv  <- "/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Brain_FC.tsv"

alpha <- 0.05 / 720 / 9      # your threshold
label_top_edges <- 12        # number of edges to annotate

# -----------------------
# Helpers
# -----------------------
infer_n_from_m <- function(m) {
  n <- (1 + sqrt(1 + 8*m)) / 2
  if (abs(n - round(n)) > 1e-6) {
    stop(sprintf("Cannot infer integer n from m=%s (n=%s). Check IDP indexing.", m, n))
  }
  as.integer(round(n))
}

make_pair_table <- function(n) {
  k <- 0L
  out <- vector("list", n*(n-1)/2)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      k <- k + 1L
      out[[k]] <- tibble(edge_index = k, node_i = i, node_j = j)
    }
  }
  bind_rows(out)
}

# -----------------------
# Load + parse
# -----------------------
df <- readr::read_tsv(in_tsv, show_col_types = FALSE) %>%
  mutate(
    edge_index = as.integer(str_remove(IDP, "^f_")),
    Log_Odds   = as.numeric(Log_Odds),
    P_value    = as.numeric(P_value)
  )

m_max <- max(df$edge_index, na.rm = TRUE)
n_nodes <- infer_n_from_m(m_max)   # for 210 -> 21
pair_tbl <- make_pair_table(n_nodes)

df2 <- df %>%
  left_join(pair_tbl, by = "edge_index")

# -----------------------
# Filter significant rows
# -----------------------
sig <- df2 %>%
  filter(!is.na(P_value), P_value < alpha)

if (nrow(sig) == 0) {
  stop("No significant results at P < 0.05/720/9. Nothing to plot.")
}

# -----------------------
# Summarize at EDGE level
# -----------------------
edge_sum <- sig %>%
  group_by(edge_index, node_i, node_j) %>%
  summarise(
    n_sig   = n(),                    # number of significant BAG–IDP rows for this edge
    min_p   = min(P_value),
    top_bag = BAG[which.min(P_value)],
    top_beta = Log_Odds[which.min(P_value)],
    .groups = "drop"
  ) %>%
  mutate(
    top_dir = case_when(
      top_beta > 0 ~ "positive",
      top_beta < 0 ~ "negative",
      TRUE ~ "zero"
    )
  )

# -----------------------
# NODE counts (“component thickness”)
# -----------------------
node_counts <- sig %>%
  select(BAG, edge_index, node_i, node_j) %>%
  pivot_longer(cols = c(node_i, node_j), names_to = "end", values_to = "node") %>%
  group_by(node) %>%
  summarise(
    n_assoc = n(),               # total significant associations incident to node
    n_bag   = n_distinct(BAG),   # distinct BAGs hitting node
    .groups = "drop"
  )

node_df <- tibble(node = 1:n_nodes) %>%
  left_join(node_counts, by = "node") %>%
  mutate(
    n_assoc = replace_na(n_assoc, 0L),
    n_bag   = replace_na(n_bag, 0L)
  ) %>%
  mutate(
    theta = 2*pi*(node - 1)/n_nodes,
    x = cos(theta),
    y = sin(theta)
  )

# -----------------------
# Join coordinates onto edges
# -----------------------
edges_plot <- edge_sum %>%
  left_join(node_df %>% select(node, x, y), by = c("node_i" = "node")) %>%
  rename(x_i = x, y_i = y) %>%
  left_join(node_df %>% select(node, x, y), by = c("node_j" = "node")) %>%
  rename(x_j = x, y_j = y) %>%
  mutate(
    curvature = 0.25 * sign(x_i*y_j - y_i*x_j),
    mlog10p = -log10(min_p)
  )

# -----------------------
# Labels (compute constant n first; avoid n() in slice_head)
# -----------------------
n_label <- min(as.integer(label_top_edges), nrow(edges_plot))

lab_edges <- edges_plot %>%
  arrange(min_p) %>%
  slice_head(n = n_label) %>%
  mutate(
    x_mid = (x_i + x_j)/2,
    y_mid = (y_i + y_j)/2,
    label = paste0(top_bag, "\n", "f_", edge_index, "\nP=", format(min_p, digits=2, scientific=TRUE))
  )

# -----------------------
# Plot
# -----------------------
p <- ggplot() +
  geom_curve(
    data = edges_plot,
    aes(
      x = x_i, y = y_i, xend = x_j, yend = y_j,
      linewidth = n_sig,
      alpha = mlog10p,
      color = top_dir,
      curvature = curvature
    ),
    lineend = "round"
  ) +
  geom_point(
    data = node_df,
    aes(x = x, y = y, size = n_assoc),
    shape = 21, fill = "white", color = "black", stroke = 0.8
  ) +
  geom_text(
    data = node_df,
    aes(x = 1.08*x, y = 1.08*y, label = node),
    size = 3
  ) +
  geom_text_repel(
    data = lab_edges,
    aes(x = x_mid, y = y_mid, label = label),
    size = 3,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(positive = "firebrick3", negative = "royalblue3", zero = "grey50")) +
  scale_linewidth(range = c(0.3, 2.5)) +
  scale_alpha(range = c(0.15, 0.9)) +
  scale_size(range = c(2.5, 8)) +
  coord_equal() +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = "Significant Brain-FC edges (UKB group-ICA d25 good-nodes style)",
    subtitle = paste0("Threshold: P < 0.05/720/9 = ", signif(alpha, 3),
                      " | Nodes sized by # significant BAG–IDP associations",
                      " | Edges colored by sign of top BAG effect"),
    color = "Top effect sign",
    linewidth = "# sig BAGs\n(per edge)",
    size = "# sig assocs\n(per node)",
    alpha = "-log10(min P)\n(per edge)"
  )

print(p)
