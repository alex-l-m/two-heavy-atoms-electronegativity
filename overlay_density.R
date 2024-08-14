# Make a csv file of a density with periodic boundary conditions, give in the
# density sampled on a grid without periodic boundary conditions
library(tidyverse)

# Output density path
output_path <- commandArgs(trailingOnly = TRUE)[1]

# Input density
density_path <- 'density.txt'
density <- read_table(density_path, skip = 1,
                      col_names = c('x', 'y', 'z', 'density'),
                      col_types = cols(
                          x = col_double(),
                          y = col_double(),
                          z = col_double(),
                          density = col_double()
                      ))

# Indices of grid points (non-unique, since they represent indices in the unit
# cell)
grid <- read_csv('grid.csv', col_types = cols(
    i = col_integer(),
    j = col_integer(),
    k = col_integer()
))

# These were output by the same script with rows in the same order. Therefore I
# should be able to just concatenate them
density_with_grid <- bind_cols(grid, density)

# Sum up the densities with the same index in the unit cell
density_summed <- density_with_grid |>
    group_by(i, j, k) |>
    summarize(density = sum(density), .groups = 'drop')

write_csv(density_summed, output_path)
