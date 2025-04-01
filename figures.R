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
hardness_regression_plot <- readr::read_rds(hardness_regression_plot_path) +
    # Change the color of the hlines and vlines to black
    geom_hline(yintercept = 0, color = 'black') +
    geom_vline(xintercept = 0, color = 'black') +
    # Rename the legend title for color to "cation"
    labs(color = 'Cation')

# Save the electronegativity comparison plot as a figure that covers half a
# powerpoint slide
ggsave('electronegativity_comparison_plot_powerpoint.png', electronegativity_comparison_plot,
       width = unit(5.68, 'in'), height = unit(4.76, 'in'))
# Save the hardness regression plot as a figure for powerpoint slide
ggsave('hardness_regression_plot_powerpoint.png', hardness_regression_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# List of energy curve plots
energy_curve_plot_list <- readr::read_rds('3-5:zincblende:0_energy_with_nofield_plots.rds')
# Example for a figure
example_energy_curve_plot <- energy_curve_plot_list[['GaP']] +
            geom_line(color = 'black') +
            ggtitle('GaP energy (in nuclear potential)')
ggsave('example_energy_curve_plot.png', example_energy_curve_plot,
       width = unit(5.68, 'in'), height = unit(4.76, 'in'))

# Same thing for energy derivative plots
electronegativity_plot_list <- readr::read_rds('3-5:zincblende:0_lam_plots.rds')
example_electronegativity_plot <- electronegativity_plot_list[['GaP']] +
            ggtitle('GaP electronegativity difference') +
            theme(legend.position = 'bottom')
ggsave('example_electronegativity_plot.png', example_electronegativity_plot,
       width = unit(5.68, 'in'), height = unit(4.76, 'in'))


# Plot layout design so that the hardness regression plot takes up two "slots"
design <- "
    ABC
    DDE
"
panels <- circuit_diagram_patchwork_element + example_energy_curve_plot + example_electronegativity_plot +
    hardness_regression_plot + electronegativity_comparison_plot +
    patchwork::plot_layout(design = design) +
    # Label each panel with capital letters
    patchwork::plot_annotation(tag_levels = 'A')

ggsave('figure_transparent.png', panels, width = unit(10, 'in'), height = unit(7.5, 'in'))
# Use ImageMagick to make the background white
system('convert figure_transparent.png -background white -alpha remove -alpha off figure.png')
