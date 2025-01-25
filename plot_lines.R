# Pot density along a line between the atoms

library(tidyverse)
library(cowplot)
library(ggdark)
theme_set(dark_mode(theme_cowplot()))
library(glue)

# Atomic numbers, for ordering
atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_integer()
))

simulations <- read_csv('simulations.csv', col_types = cols(
    simulation_id = col_character(),
    potential = col_character(),
    structure_id = col_character(),
    cation = col_character(),
    anion = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    log_file_path = col_character(),
    cube_file_path = col_character(),
    pot_file_path = col_character(),
    hartree_pot_path = col_character(),
    acceptor_proatom_path = col_character(),
    donor_proatom_path = col_character(),
    promolecule_path = col_character()
))

# Functions plotted on a line between atoms
line_tbl <- read_csv('lines.csv.gz', col_types = cols(
    cube_id = col_character(),
    i = col_double(),
    j = col_double(),
    k = col_double(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    density = col_double(),
    squarednorm = col_double(),
    projection_squarednorm = col_double(),
    simulation_id = col_character(),
    field_number = col_double(),
    structure_id = col_character(),
    cube_file_type = col_character(),
    cube_file_path = col_character()
))

# Ordered anion and cation symbols
cation_levels <- simulations |>
    select(cation) |>
    distinct() |>
    left_join(atomic_numbers, by = c('cation' = 'symbol')) |>
    arrange(atomic_number) |>
    pull(cation)
anion_levels <- simulations |>
    select(anion) |>
    distinct() |>
    left_join(atomic_numbers, by = c('anion' = 'symbol')) |>
    arrange(atomic_number) |>
    pull(anion)

# How to label the cube file types for plotting
cube_type_labels <- c(acceptor_proatom_path = 'Acceptor atom ref',
                      donor_proatom_path = 'Donor atom ref',
                      cube_file_path = 'Actual density')

# Create the output directory
line_density_dir <- 'line_density_plots'
dir.create(line_density_dir, showWarnings = FALSE)

# Make a plot for each field number
for (this_field_number in unique(line_tbl$field_number))
{
    # Plot density
    line_plt <- line_tbl |>
        filter(field_number == this_field_number) |>
        # Relabel and order the cube types
        mutate(cube_type_label_character = cube_type_labels[cube_file_type],
               cube_type_label = factor(cube_type_label_character,
                                        levels = cube_type_labels)) |>
        # Join with the simulation table so I can use the extra annotations
        left_join(simulations,
                  by = c('simulation_id', 'field_number', 'structure_id')) |>
        # Order the anion and cations
        mutate(anion = factor(anion, levels = anion_levels),
               cation = factor(cation, levels = cation_levels)) |>
        ggplot(aes(x = sqrt(projection_squarednorm), y = density, color = cube_type_label)) +
        # Rows for cations, columns for anions
        facet_grid(cation ~ anion) +
        geom_line() +
        # Plot density on a log scale so I can see what's happening in the bonds
        scale_y_log10(limits = c(.01, 10)) +
        xlab('Distance along line between atoms (Å)') +
        theme(legend.position = 'bottom',
              # Add a background grid, thin and lightly coloured
              panel.grid.major = element_line(size = 0.1, colour = 'grey80'))
    
    ggsave(glue('{line_density_dir}/line_density_{this_field_number}.png'), line_plt, height = unit(4.76, 'in'), width = unit(11.5, 'in'), dpi = 96)
}
