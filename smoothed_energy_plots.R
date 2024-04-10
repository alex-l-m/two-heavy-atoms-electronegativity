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

nofield_simulation_table <- simulation_table |>
    filter(molecule_charge == 0, is.na(field_value))

nofield_derivatives <- nofield_simulation_table |>
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
    inner_join(nofield_simulation_table, by = 'combination_id') |>
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
    geom_point(data = energies) +
    this_theme
ggsave('energy_smoothing_validation.png', energy_smoothing_validation_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

energy_with_nofield <- smoothed_energy |>
    ggplot(mapping = aes(x = charge_transfer, y = energy)) +
    facet_wrap(~ formula, scales = 'free', nrow = 2) +
    geom_smooth(method = lmrob, formula = y ~ x + I(x^2), se = FALSE) +
    geom_line() +
    geom_point(data = nofield_energies) +
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

# Make a version of the plot with the theoretical lines overlaid
ipea_fukui <- read_csv('ipea_fukui.csv.gz', col_types = cols(
  formula = col_character(),
  ip = col_double(),
  ea = col_double(),
  electronegativity = col_double(),
  hardness = col_double(),
  lower_fukui_acceptor = col_double(),
  lower_fukui_donor = col_double(),
  upper_fukui_acceptor = col_double(),
  upper_fukui_donor = col_double()
))
slopes <- ipea_fukui |>
    group_by(formula) |>
    transmute(
        slope_1 = hardness * (1/lower_fukui_acceptor + 1/upper_fukui_donor),
        slope_2 = hardness * (1/upper_fukui_acceptor + 1/lower_fukui_donor)
    ) |>
    ungroup()

theoretical_lines <- slopes |>
    inner_join(nofield_derivatives, by = 'formula') |>
    mutate(
        intercept_1 = derivative - slope_1 * charge_transfer,
        intercept_2 = derivative - slope_2 * charge_transfer
    )

energy_derivatives_with_nofield_lines <- energy_derivatives_with_nofield +
    geom_abline(aes(slope = slope_1, intercept = intercept_1), data = theoretical_lines) +
    geom_abline(aes(slope = slope_2, intercept = intercept_2), data = theoretical_lines)

ggsave('energy_derivatives_with_nofield_lines.png', energy_derivatives_with_nofield_lines, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

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
    geom_point(data = nofield_derivatives) +
    ylim(c(-5,5)) +
    this_theme
ggsave('energy_derivatives_zoomed.png', energy_derivatives_zoomed, width = unit(11.5, 'in'), height = unit(4.76, 'in'))


# Read the 'Lam' values
lam <- read_csv('lam_values.csv.gz', col_types = cols(
    combination_id = col_character(),
    lam = col_double()
))
# Doesn't have the formula column, or the charge transfer. Join with the table
# containing energies and charges, just for the charges
lam_charges <- lam |>
    inner_join(energy_charge, by = 'combination_id') |>
    select(formula, combination_id, lam, charge_donor, charge_acceptor)

# Considering it as estimates of the derivative
lam_derivs <- lam_charges |>
    mutate(charge_transfer = charge_acceptor) |>
    select(formula, combination_id, lam, charge_transfer) |>
    # Just putting the negative sign because that makes it match, I don't get
    # why. Factor of two is because my linear combination has total weight two,
    # though I thought that means I need to divide by two, but multiplying
    # makes it match
    mutate(derivative = - 2 * lam)

# Table containing both the lambda values and the energy derivatives, for
# comparison
smoothed_energy_deriv_only <- smoothed_energy |>
    select(formula, charge_transfer, derivative)
# Each Lam value comes from a specific simulation, but each energy derivative
# comes from a smoothing process. So only Lam has a combination id, but this
# won't be relevant in the comparison
lam_comparison <- bind_rows(
        `2 * Lam` = select(lam_derivs, -combination_id),
        `dE/dq` = smoothed_energy_deriv_only, .id = 'computation') |>
    # Put Lam first, so it's the same color if I plot just that in another plot
    mutate(computation = factor(computation, levels = c('2 * Lam', 'dE/dq')))

# Add the Lam estimates of derivative, from the simulations where no field was
# applied, to the table of no field derivatives for comparison
# Start with the table of no field simulations, and add the Lam estimates of
# derivatives
# These aren't in the Lam values table because CDFT was not done, so there is
# no "Lam" value reported in the log. However, assuming it's a Lagrange
# multiplier, it can be assumed to be zero
nofield_charges <- energy_charge |>
    inner_join(nofield_simulation_table, by = c('combination_id', 'formula')) |>
    select(combination_id, formula, charge_acceptor, energy) |>
    rename(charge_transfer = charge_acceptor)
lam_derivs_nofield <- nofield_charges |>
    mutate(derivative = 0) |>
    select(formula, charge_transfer, derivative)
# Bind rows to make the comparison table
nofield_deriv_comparison <- bind_rows(
        `2 * Lam` = lam_derivs_nofield,
        `dE/dq` = nofield_derivatives, .id = 'computation') |>
    # Factor with same levels as the table with all field values
    mutate(computation = factor(computation, levels = c('2 * Lam', 'dE/dq')))

lam_plot <- lam_comparison |>
    ggplot(aes(x = charge_transfer, y = derivative, color = computation)) +
    facet_wrap(~ formula, scales = 'free', nrow = 2) +
    xlim(-1, 1) +
    theme(legend.position = 'bottom') +
    # Put a vertical line to indicate 0
    geom_vline(xintercept = 0, linetype = 'dashed') +
    # Put a horizontal line to indicate zero
    geom_hline(yintercept = 0, linetype = 'dashed')+
    xlab('charge of electron acceptor') +
    ylab('electronegativity difference') +
    # Remove the title from the legend
    guides(color = guide_legend(title = NULL)) +
    geom_line() +
    # Dots for the no-field values
    geom_point(data = nofield_deriv_comparison) +
    ylim(-10, 10)
ggsave('lam_comparison.png', lam_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Maybe I should reduce duplication by saving a lot of the settings to a
# variable?
lam_plot_zoomed <- lam_comparison |>
    # Only positive values of electronegativity, only negative values of charge
    # transfer
    filter(derivative > 0, charge_transfer < 0) |>
    ggplot(aes(x = charge_transfer, y = derivative, color = computation)) +
    # Scales shouldn't be free
    facet_wrap(~ formula, nrow = 2) +
    theme(legend.position = 'bottom') +
    # Put a vertical line to indicate 0
    geom_vline(xintercept = 0, linetype = 'dashed') +
    # Put a horizontal line to indicate zero
    geom_hline(yintercept = 0, linetype = 'dashed')+
    xlab('charge of electron acceptor') +
    ylab('electronegativity difference') +
    # Remove the title from the legend
    guides(color = guide_legend(title = NULL)) +
    geom_line() +
    # Dots for the no-field values
    geom_point(data = nofield_deriv_comparison)
ggsave('lam_comparison_zoomed.png', lam_plot_zoomed,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))
