# Make a scatter plot of pairwise distances between the atoms in rocksalt vs
# zincblende scatterplot
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggdark)
library(tidymodels)
theme_set(dark_mode(theme_cowplot()))

pairwise_distance_table <- read_csv('pairwise_distance_table.csv', col_types = cols(
  structure_id = col_character(),
  formula = col_character(),
  crystal_structure = col_character(),
  distance = col_double()
))

pairwise_distance_plot <- pairwise_distance_table |>
    pivot_wider(id_cols = formula, names_from = crystal_structure, values_from = distance) |>
    ggplot(aes(y = rocksalt, x = zincblende, label = formula)) +
    geom_point() +
    geom_text_repel() +
    coord_obs_pred() +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    # Lines showing 10% changes
    geom_abline(intercept = 0, slope = 1.1, linetype = 'dashed') +
    geom_abline(intercept = 0, slope = 0.9, linetype = 'dashed') +
    ggtitle('Distance between nearest\nneighbor heteroatoms')

ggsave('pairwise_distance_plot.png', pairwise_distance_plot, 
       height = unit(4.76, 'in'), width = unit(5.67, 'in'))
