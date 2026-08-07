#!/usr/bin/env Rscript

# ==============================================================================
# Bonferroni-significant EPOCH mediation map
# Revised version: indirect pathways only, solid lines, separate a/b annotations
#
# Scientific time order:
#
#   Baseline metabolomics-based disease EPOCH
#                -- a -->   MRI mortality EPOCH
#                -- b -->   Observed post-MRI mortality
#
# Visual encoding
# ---------------
#   COLOR = organ
#   SHAPE = disease endpoint for molecular exposure nodes
#   SOLID lines = indirect mediation pathway only
#   a annotation = exposure -> MRI mortality EPOCH association
#   b annotation = MRI mortality EPOCH -> observed mortality Cox association
#   Hollow nodes = pathway/mediator with PH violation for the MRI mediator
#
# All text is black.
# Plot/panel backgrounds are transparent.
#
# Local Mac root:
#   /Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock
#
# Required R packages:
#   ggplot2, dplyr, readr, stringr, scales, grid
#
# Optional:
#   svglite
# ==============================================================================


# ==============================================================================
# USER SETTINGS
# ==============================================================================

ROOT <- "/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock"

JOB_ID <- "16997583"

INPUT_FILE <- file.path(
  ROOT,
  "mediation_OLS_Cox_bootstrap_full_single_models",
  paste0("job_", JOB_ID),
  "collected_results",
  "OLS_Cox_bootstrap_mediation_Bonferroni_significant.tsv"
)

OUTPUT_DIR <- file.path(
  ROOT,
  "mediation_OLS_Cox_bootstrap_full_single_models",
  paste0("job_", JOB_ID),
  "collected_results",
  "figures"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

PNG_FILE <- file.path(
  OUTPUT_DIR,
  "Bonferroni_significant_indirect_mediation_network_precise.png"
)

PDF_FILE <- file.path(
  OUTPUT_DIR,
  "Bonferroni_significant_indirect_mediation_network_precise.pdf"
)

SVG_FILE <- file.path(
  OUTPUT_DIR,
  "Bonferroni_significant_indirect_mediation_network_precise.svg"
)


# ==============================================================================
# PACKAGES
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(scales)
  library(grid)
})


# ==============================================================================
# LOAD DATA
# ==============================================================================

if (!file.exists(INPUT_FILE)) {
  stop(
    paste0(
      "\nInput file not found:\n",
      INPUT_FILE,
      "\n\nPlease update ROOT or JOB_ID at the top of the script."
    )
  )
}

dat <- read_tsv(
  INPUT_FILE,
  show_col_types = FALSE,
  progress = FALSE
)

required_cols <- c(
  "model_id",
  "exposure_organ",
  "exposure_modality",
  "exposure_endpoint",
  "mediator_organ",
  "a_beta",
  "a_p_hc3",
  "b_log_hr",
  "b_hr",
  "b_p",
  "indirect_log_hr",
  "indirect_delta_p",
  "indirect_bonferroni_p",
  "boot_ci_low",
  "boot_ci_high"
)

missing_cols <- setdiff(
  required_cols,
  names(dat)
)

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns:",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}

# PH flag may already exist in the collected table.
if (
  !"ph_mediator_violation_p_lt_0_05" %in% names(dat)
  && !"ph_mediator_p" %in% names(dat)
) {
  warning(
    "Neither ph_mediator_violation_p_lt_0_05 nor ph_mediator_p is present. ",
    "PH-violation hollow-node encoding will be unavailable."
  )
}

dat <- dat %>%
  filter(
    is.finite(indirect_log_hr),
    is.finite(indirect_delta_p)
  ) %>%
  mutate(
    exposure_organ = as.character(exposure_organ),
    mediator_organ = as.character(mediator_organ),
    disease = as.character(exposure_endpoint),
    modality = as.character(exposure_modality),
    
    ph_mediator_violation = case_when(
      "ph_mediator_violation_p_lt_0_05" %in% names(.) ~
        as.logical(ph_mediator_violation_p_lt_0_05),
      
      "ph_mediator_p" %in% names(.) ~
        is.finite(ph_mediator_p) & ph_mediator_p < 0.05,
      
      TRUE ~ FALSE
    ),
    
    exposure_node_id = paste(
      exposure_organ,
      modality,
      disease,
      sep = "__"
    ),
    
    mediator_node_id = paste0(
      mediator_organ,
      "__mri_mortality"
    ),
    
    exposure_label = paste0(
      str_to_title(
        str_replace_all(
          exposure_organ,
          "_",
          " "
        )
      ),
      "\n",
      str_to_title(disease),
      " EPOCH"
    ),
    
    mediator_label = paste0(
      str_to_title(
        str_replace_all(
          mediator_organ,
          "_",
          " "
        )
      ),
      "\nMRI mortality EPOCH"
    ),
    
    indirect_abs = abs(indirect_log_hr)
  )

if (nrow(dat) == 0) {
  stop(
    "No finite Bonferroni-significant mediation results were found."
  )
}

non_metabolomics <- dat %>%
  filter(
    tolower(modality) != "metabolomics"
  )

if (nrow(non_metabolomics) > 0) {
  warning(
    paste0(
      nrow(non_metabolomics),
      " significant row(s) are not metabolomics. ",
      "The left-column heading is metabolomics-based because that is the ",
      "current significant result set."
    )
  )
}


# ==============================================================================
# VAN GOGH-INSPIRED ORGAN PALETTE
# ==============================================================================

organ_palette_master <- c(
  "brain"               = "#355C9A",
  "heart"               = "#D1495B",
  "spleen"              = "#E6B325",
  "digestive"           = "#2A9D8F",
  "endocrine"           = "#7251B5",
  "hepatic"             = "#D68C1F",
  "immune"              = "#4F7D4A",
  "metabolic"           = "#277DA1",
  "pulmonary"           = "#64B5CD",
  "reproductive_female" = "#C77DAB",
  "reproductive_male"   = "#8A6D46",
  "adipose"             = "#E9C46A",
  "kidney"              = "#668F80",
  "liver"               = "#B9770E",
  "pancreas"            = "#E76F51"
)

all_organs <- unique(
  c(
    dat$exposure_organ,
    dat$mediator_organ
  )
)

missing_palette_organs <- setdiff(
  all_organs,
  names(organ_palette_master)
)

if (length(missing_palette_organs) > 0) {
  
  extra_cols <- hue_pal(
    h.start = 20,
    c = 90,
    l = 55
  )(
    length(missing_palette_organs)
  )
  
  names(extra_cols) <- missing_palette_organs
  
  organ_palette_master <- c(
    organ_palette_master,
    extra_cols
  )
}

organ_palette <- organ_palette_master[
  all_organs
]


# ==============================================================================
# DISEASE SHAPES
# ==============================================================================

disease_shape_master <- c(
  "asthma"   = 21,
  "dementia" = 22,
  "mi"       = 24,
  "stroke"   = 23,
  "copd"     = 25
)

diseases <- unique(
  dat$disease
)

missing_diseases <- setdiff(
  diseases,
  names(disease_shape_master)
)

if (length(missing_diseases) > 0) {
  
  available_shapes <- c(
    0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12
  )
  
  extra_shapes <- available_shapes[
    seq_len(
      min(
        length(missing_diseases),
        length(available_shapes)
      )
    )
  ]
  
  names(extra_shapes) <- missing_diseases[
    seq_along(extra_shapes)
  ]
  
  disease_shape_master <- c(
    disease_shape_master,
    extra_shapes
  )
}

disease_shapes <- disease_shape_master[
  diseases
]


# ==============================================================================
# NODE LAYOUT
# ==============================================================================

# Exposure nodes are collapsed by organ + disease + modality.
exposure_nodes <- dat %>%
  distinct(
    exposure_node_id,
    exposure_organ,
    disease,
    modality,
    exposure_label
  ) %>%
  arrange(
    disease,
    exposure_organ,
    modality
  ) %>%
  mutate(
    x = 1.00,
    y = rev(
      seq_len(n())
    ) / (n() + 1)
  )

# Mark an exposure node hollow if ANY significant pathway from that exposure has
# a mediator PH violation. This highlights the spleen-mediated MI signals.
exposure_ph_flag <- dat %>%
  group_by(
    exposure_node_id
  ) %>%
  summarise(
    ph_violation_any = any(
      ph_mediator_violation,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

exposure_nodes <- exposure_nodes %>%
  left_join(
    exposure_ph_flag,
    by = "exposure_node_id"
  ) %>%
  mutate(
    ph_violation_any = ifelse(
      is.na(ph_violation_any),
      FALSE,
      ph_violation_any
    )
  )


mediator_nodes <- dat %>%
  distinct(
    mediator_node_id,
    mediator_organ,
    mediator_label
  ) %>%
  arrange(
    mediator_organ
  ) %>%
  mutate(
    x = 2.00,
    y = rev(
      seq_len(n())
    ) / (n() + 1)
  )

mediator_ph_flag <- dat %>%
  group_by(
    mediator_node_id
  ) %>%
  summarise(
    ph_violation_any = any(
      ph_mediator_violation,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

mediator_nodes <- mediator_nodes %>%
  left_join(
    mediator_ph_flag,
    by = "mediator_node_id"
  ) %>%
  mutate(
    ph_violation_any = ifelse(
      is.na(ph_violation_any),
      FALSE,
      ph_violation_any
    )
  )


mortality_node <- tibble(
  x = 3.00,
  y = 0.50
)


# ==============================================================================
# PATH-SPECIFIC EDGE LAYOUT
#
# We draw only the mediation path:
#
#   X --a--> M --b--> Y
#
# Every significant model has its own X->M edge.
# Every significant model also gets its own M->Y edge with tiny curvature offsets,
# so b annotations can be shown separately and precisely without geom_repel.
# ==============================================================================

path_dat <- dat %>%
  
  left_join(
    exposure_nodes %>%
      select(
        exposure_node_id,
        x_x = x,
        y_x = y
      ),
    by = "exposure_node_id"
  ) %>%
  
  left_join(
    mediator_nodes %>%
      select(
        mediator_node_id,
        x_m = x,
        y_m = y
      ),
    by = "mediator_node_id"
  ) %>%
  
  group_by(
    mediator_node_id
  ) %>%
  
  arrange(
    disease,
    exposure_organ,
    .by_group = TRUE
  ) %>%
  
  mutate(
    mediator_path_rank = row_number(),
    mediator_path_n = n(),
    
    # Small deterministic fan-out for M -> Y pathways.
    b_curve = ifelse(
      mediator_path_n == 1,
      0,
      seq(
        from = -0.16,
        to = 0.16,
        length.out = mediator_path_n
      )
    )
  ) %>%
  
  ungroup() %>%
  
  mutate(
    x_y = mortality_node$x,
    y_y = mortality_node$y,
    
    # -------------------------------------------------------------------------
    # PRECISE A-LABEL LOCATION
    #
    # Put a labels at 43% of the X->M segment, slightly above/below the line.
    # Offset is deterministic by model order to avoid identical coordinates.
    # -------------------------------------------------------------------------
    
    a_label_x = x_x + 0.43 * (x_m - x_x),
    
    a_base_y = y_x + 0.43 * (y_m - y_x),
    
    a_label_offset = case_when(
      row_number() %% 3 == 0 ~  0.018,
      row_number() %% 3 == 1 ~ -0.018,
      TRUE                  ~  0.000
    ),
    
    a_label_y = a_base_y + a_label_offset,
    
    a_label = paste0(
      "a = ",
      sprintf(
        "%.3f",
        a_beta
      )
    ),
    
    # -------------------------------------------------------------------------
    # PRECISE B-LABEL LOCATION
    #
    # Place b labels at 55% of the M->Y segment.
    # Add a controlled vertical offset based on pathway rank within each mediator.
    # -------------------------------------------------------------------------
    
    b_label_x = x_m + 0.55 * (x_y - x_m),
    
    b_base_y = y_m + 0.55 * (y_y - y_m)
  ) %>%
  
  group_by(
    mediator_node_id
  ) %>%
  
  mutate(
    b_label_offset = ifelse(
      n() == 1,
      0,
      seq(
        from = -0.085,
        to = 0.085,
        length.out = n()
      )
    ),
    
    b_label_y = b_base_y + b_label_offset,
    
    b_label = paste0(
      "b = ",
      sprintf(
        "%.3f",
        b_log_hr
      )
    )
  ) %>%
  
  ungroup()


# ==============================================================================
# COLUMN HEADINGS + DATA COLLECTION TIME
# ==============================================================================

column_headers <- tibble(
  x = c(
    1,
    2,
    3
  ),
  
  title_y = c(
    1.095,
    1.095,
    1.095
  ),
  
  time_y = c(
    1.020,
    1.020,
    1.020
  ),
  
  title = c(
    "METABOLOMICS-BASED\nDISEASE EPOCH",
    "MRI MORTALITY EPOCH",
    "OBSERVED MORTALITY"
  ),
  
  time = c(
    "Baseline molecular assessment\n(before MRI)",
    "Imaging visit 2\n(MRI assessment)",
    "Post-MRI follow-up\n(death or administrative censoring)"
  )
)


# ==============================================================================
# LEGEND DATA FOR PH VIOLATION
# ==============================================================================

ph_legend <- tibble(
  x = c(
    0.78,
    0.78
  ),
  
  y = c(
    -0.015,
    -0.055
  ),
  
  label = c(
    "Filled node: no mediator PH violation",
    "Hollow node: mediator PH violation (p < 0.05)"
  ),
  
  hollow = c(
    FALSE,
    TRUE
  )
)


# ==============================================================================
# PLOT
# ==============================================================================

p <- ggplot() +
  
  # ---------------------------------------------------------------------------
# COLUMN SEPARATORS
# ---------------------------------------------------------------------------

geom_vline(
  xintercept = c(
    1.5,
    2.5
  ),
  linewidth = 0.35,
  color = "black",
  alpha = 0.16,
  linetype = "dashed"
) +
  
  
  # ---------------------------------------------------------------------------
# INDIRECT PATH X -> M
# SOLID
# ---------------------------------------------------------------------------

geom_curve(
  data = path_dat,
  
  aes(
    x = x_x,
    y = y_x,
    xend = x_m,
    yend = y_m,
    color = exposure_organ,
    linewidth = indirect_abs,
    group = model_id
  ),
  
  curvature = 0.10,
  linetype = "solid",
  alpha = 0.85,
  lineend = "round",
  
  arrow = arrow(
    length = unit(
      2.5,
      "mm"
    ),
    type = "closed"
  )
) +
  
  
  # ---------------------------------------------------------------------------
# INDIRECT PATH M -> Y
# SOLID
#
# One model-specific curve per significant pathway.
# ---------------------------------------------------------------------------

geom_curve(
  data = path_dat,
  
  aes(
    x = x_m,
    y = y_m,
    xend = x_y,
    yend = y_y,
    color = mediator_organ,
    linewidth = indirect_abs,
    group = model_id
  ),
  
  curvature = 0.08,
  linetype = "solid",
  alpha = 0.78,
  lineend = "round",
  
  arrow = arrow(
    length = unit(
      2.5,
      "mm"
    ),
    type = "closed"
  )
) +
  
  
  # ---------------------------------------------------------------------------
# A-PATH ANNOTATIONS
# Exact fixed coordinates; no repel.
# ---------------------------------------------------------------------------

geom_label(
  data = path_dat,
  
  aes(
    x = a_label_x,
    y = a_label_y,
    label = a_label
  ),
  
  color = "black",
  fill = alpha(
    "white",
    0.82
  ),
  
  label.size = 0.18,
  label.r = unit(
    0.08,
    "lines"
  ),
  
  label.padding = unit(
    0.10,
    "lines"
  ),
  
  size = 2.55,
  fontface = "bold",
  show.legend = FALSE
) +
  
  
  # ---------------------------------------------------------------------------
# B-PATH ANNOTATIONS
# Exact fixed coordinates; no repel.
# ---------------------------------------------------------------------------

geom_label(
  data = path_dat,
  
  aes(
    x = b_label_x,
    y = b_label_y,
    label = b_label
  ),
  
  color = "black",
  fill = alpha(
    "white",
    0.82
  ),
  
  label.size = 0.18,
  label.r = unit(
    0.08,
    "lines"
  ),
  
  label.padding = unit(
    0.10,
    "lines"
  ),
  
  size = 2.55,
  fontface = "bold",
  show.legend = FALSE
) +
  
  
  # ---------------------------------------------------------------------------
# EXPOSURE NODES WITHOUT PH VIOLATION
# FILLED
# ---------------------------------------------------------------------------

geom_point(
  data = exposure_nodes %>%
    filter(
      !ph_violation_any
    ),
  
  aes(
    x = x,
    y = y,
    shape = disease,
    fill = exposure_organ
  ),
  
  color = "black",
  stroke = 1.05,
  size = 6.5
) +
  
  
  # ---------------------------------------------------------------------------
# EXPOSURE NODES WITH PH VIOLATION
# HOLLOW
#
# For the current result set these highlight the spleen-mediated MI signals.
# ---------------------------------------------------------------------------

geom_point(
  data = exposure_nodes %>%
    filter(
      ph_violation_any
    ),
  
  aes(
    x = x,
    y = y,
    shape = disease
  ),
  
  fill = "white",
  color = "black",
  stroke = 1.25,
  size = 6.5
) +
  
  
  # ---------------------------------------------------------------------------
# EXPOSURE LABELS
# ---------------------------------------------------------------------------

geom_text(
  data = exposure_nodes,
  
  aes(
    x = x - 0.055,
    y = y,
    label = exposure_label
  ),
  
  hjust = 1,
  color = "black",
  size = 3.35,
  lineheight = 0.90,
  fontface = "bold"
) +
  
  
  # ---------------------------------------------------------------------------
# MRI MEDIATOR NODES WITHOUT PH VIOLATION
# FILLED
# ---------------------------------------------------------------------------

geom_point(
  data = mediator_nodes %>%
    filter(
      !ph_violation_any
    ),
  
  aes(
    x = x,
    y = y,
    fill = mediator_organ
  ),
  
  shape = 22,
  color = "black",
  stroke = 1.15,
  size = 8.2
) +
  
  
  # ---------------------------------------------------------------------------
# MRI MEDIATOR NODES WITH PH VIOLATION
# HOLLOW
# ---------------------------------------------------------------------------

geom_point(
  data = mediator_nodes %>%
    filter(
      ph_violation_any
    ),
  
  aes(
    x = x,
    y = y
  ),
  
  shape = 22,
  fill = "white",
  color = "black",
  stroke = 1.35,
  size = 8.2
) +
  
  
  # ---------------------------------------------------------------------------
# MRI LABELS
# ---------------------------------------------------------------------------

geom_text(
  data = mediator_nodes,
  
  aes(
    x = x,
    y = y - 0.060,
    label = mediator_label
  ),
  
  color = "black",
  size = 3.50,
  lineheight = 0.90,
  fontface = "bold",
  vjust = 1
) +
  
  
  # ---------------------------------------------------------------------------
# OBSERVED MORTALITY NODE
# ---------------------------------------------------------------------------

geom_point(
  data = mortality_node,
  
  aes(
    x = x,
    y = y
  ),
  
  shape = 8,
  color = "#E6B325",
  size = 10,
  stroke = 1.5
) +
  
  geom_point(
    data = mortality_node,
    
    aes(
      x = x,
      y = y
    ),
    
    shape = 8,
    color = "black",
    size = 7.2,
    stroke = 1.2
  ) +
  
  geom_text(
    data = mortality_node,
    
    aes(
      x = x,
      y = y - 0.075
    ),
    
    label = "Observed\nall-cause mortality",
    color = "black",
    size = 4.05,
    fontface = "bold",
    lineheight = 0.90,
    vjust = 1
  ) +
  
  
  # ---------------------------------------------------------------------------
# COLUMN TITLES
# ---------------------------------------------------------------------------

geom_text(
  data = column_headers,
  
  aes(
    x = x,
    y = title_y,
    label = title
  ),
  
  color = "black",
  size = 4.15,
  fontface = "bold",
  lineheight = 0.92
) +
  
  
  # ---------------------------------------------------------------------------
# DATA COLLECTION TIME
# ---------------------------------------------------------------------------

geom_text(
  data = column_headers,
  
  aes(
    x = x,
    y = time_y,
    label = time
  ),
  
  color = "black",
  size = 3.0,
  fontface = "italic",
  lineheight = 0.95
) +
  
  
  # ---------------------------------------------------------------------------
# SMALL PH-VIOLATION EXPLANATION
# ---------------------------------------------------------------------------

geom_point(
  data = ph_legend %>%
    filter(
      !hollow
    ),
  
  aes(
    x = x,
    y = y
  ),
  
  shape = 22,
  fill = "#E6B325",
  color = "black",
  size = 4.3,
  stroke = 1
) +
  
  geom_point(
    data = ph_legend %>%
      filter(
        hollow
      ),
    
    aes(
      x = x,
      y = y
    ),
    
    shape = 22,
    fill = "white",
    color = "black",
    size = 4.3,
    stroke = 1.1
  ) +
  
  geom_text(
    data = ph_legend,
    
    aes(
      x = x + 0.035,
      y = y,
      label = label
    ),
    
    hjust = 0,
    color = "black",
    size = 2.75
  ) +
  
  
  # ---------------------------------------------------------------------------
# SCALES
# ---------------------------------------------------------------------------

scale_color_manual(
  values = organ_palette,
  name = "Organ"
) +
  
  scale_fill_manual(
    values = organ_palette,
    name = "Organ"
  ) +
  
  scale_shape_manual(
    values = disease_shapes,
    name = "Disease"
  ) +
  
  scale_linewidth_continuous(
    range = c(
      0.8,
      3.7
    ),
    name = "|Indirect effect|"
  ) +
  
  
  # ---------------------------------------------------------------------------
# PLOT RANGE
# ---------------------------------------------------------------------------

coord_cartesian(
  xlim = c(
    0.62,
    3.28
  ),
  
  ylim = c(
    -0.09,
    1.15
  ),
  
  clip = "off"
) +
  
  
  # ---------------------------------------------------------------------------
# TITLE / CAPTION
# ---------------------------------------------------------------------------

labs(
  title = "Metabolomics-based disease EPOCHs converge on MRI mortality pathways",
  
  subtitle = paste0(
    "Bonferroni-significant indirect associations (n = ",
    nrow(dat),
    ")  |  solid arrows = mediation path  |  color = organ  |  shape = disease"
  ),
  
  caption = paste0(
    "a = association of baseline metabolomics-based disease EPOCH with later MRI mortality EPOCH; ",
    "b = Cox log-hazard coefficient for the MRI mortality EPOCH conditional on exposure and covariates. ",
    "Hollow nodes indicate pathways involving an MRI mediator with proportional-hazards p < 0.05. ",
    "Statistical mediation does not establish a causal natural indirect effect."
  )
) +
  
  
  # ---------------------------------------------------------------------------
# LEGENDS
# ---------------------------------------------------------------------------

guides(
  fill = guide_legend(
    order = 1,
    override.aes = list(
      shape = 21,
      size = 5,
      color = "black"
    )
  ),
  
  color = "none",
  
  shape = guide_legend(
    order = 2,
    override.aes = list(
      fill = "white",
      color = "black",
      size = 5
    )
  ),
  
  linewidth = guide_legend(
    order = 3
  )
) +
  
  
  # ---------------------------------------------------------------------------
# TRANSPARENT THEME / BLACK TEXT
# ---------------------------------------------------------------------------

theme_void(
  base_size = 12
) +
  
  theme(
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    
    legend.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    
    legend.box.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    
    legend.key = element_rect(
      fill = "transparent",
      color = NA
    ),
    
    plot.title = element_text(
      color = "black",
      size = 20,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 6
      )
    ),
    
    plot.subtitle = element_text(
      color = "black",
      size = 11.5,
      hjust = 0.5,
      margin = margin(
        b = 12
      )
    ),
    
    plot.caption = element_text(
      color = "black",
      size = 8.5,
      hjust = 0,
      lineheight = 1.15,
      margin = margin(
        t = 12
      )
    ),
    
    legend.position = "bottom",
    legend.box = "vertical",
    
    legend.title = element_text(
      color = "black",
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      color = "black",
      size = 9
    ),
    
    text = element_text(
      color = "black"
    ),
    
    plot.margin = margin(
      t = 25,
      r = 60,
      b = 35,
      l = 155
    )
  )


# ==============================================================================
# SAVE
# ==============================================================================

ggsave(
  filename = PNG_FILE,
  plot = p,
  width = 16.5,
  height = 11.0,
  units = "in",
  dpi = 400,
  bg = "transparent"
)

ggsave(
  filename = PDF_FILE,
  plot = p,
  width = 16.5,
  height = 11.0,
  units = "in",
  device = cairo_pdf,
  bg = "transparent"
)

if (
  requireNamespace(
    "svglite",
    quietly = TRUE
  )
) {
  
  ggsave(
    filename = SVG_FILE,
    plot = p,
    width = 16.5,
    height = 11.0,
    units = "in",
    device = svglite::svglite,
    bg = "transparent"
  )
}


# ==============================================================================
# DISPLAY + TERMINAL SUMMARY
# ==============================================================================

print(p)

cat("\n")
cat("============================================================\n")
cat("Indirect mediation figure completed\n")
cat("============================================================\n")
cat("Input: ", INPUT_FILE, "\n", sep = "")
cat(
  "Bonferroni-significant models plotted: ",
  nrow(dat),
  "\n",
  sep = ""
)
cat(
  "Models with mediator PH violation: ",
  sum(
    dat$ph_mediator_violation,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)
cat("PNG:   ", PNG_FILE, "\n", sep = "")
cat("PDF:   ", PDF_FILE, "\n", sep = "")

if (
  requireNamespace(
    "svglite",
    quietly = TRUE
  )
) {
  cat("SVG:   ", SVG_FILE, "\n", sep = "")
}

cat("============================================================\n")