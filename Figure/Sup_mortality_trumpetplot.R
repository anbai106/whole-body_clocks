.libPaths('/home/hao/R/x86_64-pc-linux-gnu-library/4.3')
library('TrumpetPlots')
library(ggplot2)
library(tidyverse)
library(ggridges)
library(data.table)
library(gridExtra)
rm(list = ls()) ## Clear global environment

# Define organ colors
organ_color <- c("brain"="#619C58",
                 "adipose"="#387AA4",
                 "heart" = "#4A9078",
                 "kidney" ="#EB5F2C",
                 "liver" ="#CEBB30",
                 "pancreas" ="#5A2657",
                 "spleen" = "#D4AEBF"
)

# List of organs and corresponding files
organs <- names(organ_color)

# Initialize empty list to store plots
plot_list <- list()

# Generate trumpet plots for each organ
for (i in seq_along(organs)) {
  organ <- organs[i]
  file <- paste0('/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/Result/TrumpetPlots_input_', organ, '.tsv')
  
  # Check if file exists
  if (!file.exists(file)) {
    warning(paste("File not found for organ:", organ))
    next
  }
  
  # Read data
  df <- fread(file)
  
  # Verify color
  organ_col <- organ_color[[organ]]
  if (is.null(organ_col)) {
    warning(paste("Color not found for organ:", organ))
    next
  }
  
  # Create trumpet plot for the current organ
  p <- plot_trumpets(
    dataset = df,
    calculate_power = FALSE,
    show_power_curves = FALSE,
    analysis_color_palette = c(organ_col)
  ) +
    ggtitle(organ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Add plot to the list
  plot_list[[organ]] <- p
}

# Combine all plots into a single figure with 11 subfigures
if (length(plot_list) > 0) {
  grid.arrange(grobs = plot_list, ncol = 3)
} else {
  warning("No plots generated due to missing files or color issues.")
}

print('Stop')