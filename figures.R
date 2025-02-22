# Assemble panels into figures, using saved raster images and RDS ggplots

library(tidyverse)
library(patchwork)
library(grid)
library(png)
library(glue)
library(cowplot)
theme_set(theme_cowplot())

# Load my circuit diagram image as a raster grob
circuit_diagram_path <- 'da_circuit.png'
circuit_diagram_grob <- rasterGrob(readPNG(circuit_diagram_path))
circuit_diagram_patchwork_element <- wrap_elements(circuit_diagram_grob)
# Just testing if it loaded properly
ggsave('circuit_diagram_patchwork_element.png', circuit_diagram_patchwork_element, width = 4, height = 4)

# Load my electronegativity plots that I serialized
category_structure_pair <- '3-5:zincblende'
electronegativity_comparison_plot_path <- glue('{category_structure_pair}_electronegativity_comparison_plot.rds')
electronegativity_comparison_plot <- readr::read_rds(electronegativity_comparison_plot_path)
hardness_regression_plot_path <- glue('{category_structure_pair}_hardness_regression_plot.rds')
hardness_regression_plot <- readr::read_rds(hardness_regression_plot_path)

# Combine in this layout:
# circuit diagram | empty space
# hardness regression (full row)
# electronegativity comparison | empty space
panels <- (circuit_diagram_patchwork_element | plot_spacer()) / 
  hardness_regression_plot / 
  (electronegativity_comparison_plot | plot_spacer())
ggsave('figure_transparent.png', panels, width = 8, height = 6)
# Use ImageMagick to make the background white
system('convert figure_transparent.png -background white -alpha remove -alpha off figure.png')
