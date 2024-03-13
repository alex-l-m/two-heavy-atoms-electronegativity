# Combine all densities into separate compressed files for each molecule,
# containing density differences from the no-field condition
library(tidyverse)
library(glue)

simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_integer(),
    field_value = col_double(),
    molecule_charge = col_integer(),
    gamess_input_file = col_character()
))

# Load and save each file separately in a loop, because if I combine them all
# into one file, I run out of memory or something when splitting it
# Make a list of simulation tables for each individual molecule (labelled by a
# formula, since there should be no isomers) to use as the basis of the loop
simulation_list <- simulation_table |>
    filter(molecule_charge == 0) |>
    group_by(formula) |>
    group_split()

header <- c('x', 'y', 'z', 'coordinate_1', 'coordinate_2', 'density')
for (this_simulation_table in simulation_list) {
    this_formula <- this_simulation_table$formula[1]
    density_plane <- this_simulation_table |>
        mutate(
            inpath = glue('densities_plane/{combination_id}.txt'),
            combination_id = str_extract(inpath, '^densities_plane/(.*).txt$', 1)
        ) |>
        # Not all simulations succeeded. Therefore, filter for simulations where
        # the inpath exists
        filter(file.exists(inpath)) |>
        group_by(combination_id, field_value) |>
        reframe(read_table(inpath,
                           col_names = header,
                           col_types = cols(.default = col_double())))
    density_plane_nofield <- density_plane |>
        filter(field_value == 0) |>
        select(x, y, z, density)
    density_difference <- density_plane |>
        left_join(density_plane_nofield, by = c('x', 'y', 'z'), suffix = c('', '_nofield')) |>
        mutate(density_difference = density - density_nofield)
    outpath <- glue('density_difference_tables/{this_formula}.csv.gz')
    write_csv(density_difference, outpath)
}
