# Plot a vector field representing a gradient that I took
library(tidyverse)
library(glue)
library(ggrepel)
library(cowplot)
library(ggdark)
theme_set(dark_mode(theme_nothing()))
library(gganimate)

grad_slices <- read_csv('slices.csv.gz', col_types = cols(
    cube_id = col_character(),
    i = col_integer(),
    j = col_integer(),
    k = col_integer(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    density = col_double(),
    simulation_id = col_character(),
    field_number = col_integer(),
    structure_id = col_character(),
    cube_file_type = col_character(),
    cube_file_path = col_character()
)) |>
    # Filter for only the gradient components
    filter(cube_file_type %in% c('dhartree_grad_x', 'dhartree_grad_y', 'dhartree_grad_z')) |>
    # Temporary: filter for only gallium arsenide, until the script is modified
    # so it loops over structures
    filter(structure_id == 'zincblende_Ga_As')

gradient <- grad_slices |>
    select(simulation_id, field_number, x, y, z, cube_file_type, density) |>
    pivot_wider(names_from = cube_file_type, values_from = density) |>
    # Calculate the norm of the gradient. I might use this for filtering later
    mutate(grad_norm = sqrt(dhartree_grad_x^2 + dhartree_grad_y^2 + dhartree_grad_z^2))

# Read the atoms
atoms <- read_csv('atoms.csv', col_types = cols(
    atom_id = col_integer(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    symbol = col_character()
))

text_table <- gradient |>
    distinct(simulation_id, 
             field_number)
# Recipe:
# https://r-graphics.org/recipe-miscgraph-vectorfield
# Vector field plot
# Value to scale gradients by
s <- 100
vfield <- ggplot() +
    geom_segment(data = gradient, aes(x = x, y = y, xend = x + dhartree_grad_x/s, yend = y + dhartree_grad_y/s),
                 arrow = arrow(length = unit(0.1/3, 'cm'))) +
    # Label with the element symbol
    geom_text(data = atoms, aes(x = x, y = y, label = symbol), box.padding = 0.5, color = 'red') +
    # Animate
    # Manual frames with no tweening
    transition_manual(field_number) +
    # Include the simulation id as text
    geom_text(data = text_table, aes(label = simulation_id), x = 0, y = -.1, color = 'red')

#ggsave('vfield.png', vfield, height = unit(4.76, 'in'), width = unit(5.67, 'in'))
anim_save('gradient_field.gif',
          vfield, duration = 5, width = 500)
