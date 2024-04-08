library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
theme_set(dark_mode(theme_cowplot()))

smoothed_energy <- read_csv('smoothed_energy.csv.gz', col_types = cols(
    formula = col_character(),
    symbol_donor = col_character(),
    symbol_acceptor = col_character(),
    charge_transfer = col_double(),
    energy = col_double(),
    derivative = col_double(),
    .default = col_double()
))

energy_derivatives <- read_csv('energy_derivatives.csv.gz', col_types = cols(
    combination_id = col_character(),
    charge_transfer = col_double(),
    derivative = col_double()
))

# Energy derivatives table is minimalistic, just has the x and y values of a
# function. I need the formula for grouping a plot, and I want to filter for
# the no field molecule, so I need the simulation data table
simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer()
))

nofield_derivatives <- simulation_table |>
    filter(molecule_charge == 0, field_value == 0) |>
    inner_join(energy_derivatives, by = 'combination_id') |>
    select(combination_id, formula, charge_transfer, derivative)

# Do the same for ordinary energies
energy_charge <- read_csv('energy_charge.csv.gz', col_types = cols(
    combination_id = col_character(),
    energy = col_double(),
    formula = col_character(),
    charge_acceptor = col_double(),
    charge_donor = col_double(),
    symbol_acceptor = col_character(),
    symbol_donor = col_character()
))

nofield_energies <- energy_charge |>
    select(combination_id, energy, charge_acceptor) |>
    inner_join(simulation_table, by = 'combination_id') |>
    filter(molecule_charge == 0, field_value == 0) |>
    select(combination_id, formula, charge_acceptor, energy) |>
    rename(charge_transfer = charge_acceptor)

# Also make a table of all energies with the same conventions for naming
# variables, for judging the fit of the smoothed curve
energies <- energy_charge |>
    select(combination_id, energy, charge_acceptor) |>
    inner_join(simulation_table, by = 'combination_id') |>
    filter(molecule_charge == 0) |>
    select(combination_id, formula, charge_acceptor, energy) |>
    rename(charge_transfer = charge_acceptor)

this_theme <- 
    theme(
        # x axis text is too crowded, rotate it
        axis.text.x = element_text(angle = 90)
    )

# Make a plot judging the quality of fit of the smoothed energy function to the
# original function
energy_smoothing_validation_plot <- smoothed_energy |>
    ggplot(mapping = aes(x = charge_transfer, y = energy)) +
    facet_wrap(~ formula, scales = 'free', nrow = 2) +
    geom_line() +
    geom_point(mapping = aes(x = charge_transfer, y = energy),
               data = energies) +
    this_theme
ggsave('energy_smoothing_validation.png', energy_smoothing_validation_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

energy_with_nofield <- smoothed_energy |>
    ggplot(mapping = aes(x = charge_transfer, y = energy)) +
    facet_wrap(~ formula, scales = 'free', nrow = 2) +
    geom_smooth(method = lmrob, formula = y ~ x + I(x^2), se = FALSE) +
    geom_line() +
    geom_point(mapping = aes(x = charge_transfer, y = energy),
               data = nofield_energies) +
    this_theme
ggsave('energy_with_nofield.png', energy_with_nofield, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

energy_derivatives_with_nofield <- smoothed_energy |>
    ggplot(aes(x = charge_transfer, y = derivative)) +
    facet_wrap(vars(formula), scales = 'free', nrow = 2) +
    geom_line() +
    geom_point(mapping = aes(x = charge_transfer, y = derivative),
               data = nofield_derivatives) +
    this_theme
ggsave('energy_derivatives_with_nofield.png', energy_derivatives_with_nofield, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Zoom in to look for linearity near neutral molecule
energy_derivatives_zoomed <- smoothed_energy |>
    ggplot(aes(x = charge_transfer, y = derivative)) +
    facet_wrap(vars(formula), nrow = 2) +
    # Add a horizontal line so we can see how far the no-field condition is
    # from zero
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # Mark equalized charge with a vertical line
    geom_vline(xintercept = 0, linetype = 'dashed') +
    # Add a linear fit to check how close to linear they are
    geom_smooth(method = lmrob, formula = y ~ x, se = FALSE) +
    geom_line() +
    geom_point(mapping = aes(x = charge_transfer, y = derivative), data = nofield_derivatives) +
    ylim(c(-5,5)) +
    this_theme
ggsave('energy_derivatives_zoomed.png', energy_derivatives_zoomed, width = unit(11.5, 'in'), height = unit(4.76, 'in'))
