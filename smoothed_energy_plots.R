library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
theme_set(dark_mode(theme_cowplot()))
library(glue)
library(ggrepel)
library(broom)

# Smoothed energies and derivatives
smoothed_energy <- read_csv('smoothed_energy.csv.gz', col_types = cols(
    formula = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    group_id = col_character(),
    charge = col_double(),
    total_energy = col_double()
))

smoothed_energy_derivative <- read_csv('smoothed_energy_derivatives.csv.gz',
                                       col_types = cols(
    formula = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    group_id = col_character(),
    charge = col_double(),
    total_energy = col_double()
))

# Raw energies and derivatives
charge_energy <- read_csv('charge_energy.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    charge = col_double(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    total_energy = col_double(),
    cdft = col_logical(),
    group_id = col_character()
))

energy_derivatives <- read_csv('energy_derivatives.csv.gz', col_types = cols(
  combination_id = col_character(),
  symbol = col_character(),
  other_symbol = col_character(),
  donor_or_acceptor = col_character(),
  group_id = col_character(),
  charge = col_double(),
  formula = col_character(),
  cdft = col_logical(),
  field_number = col_double(),
  field_value = col_double(),
  total_energy = col_double()
))

# Energies for the simulations where no field was applied, so I can add them as
# a dot on the plots. I expect them to be at the minimum of the parabola
nofield_energies <- charge_energy |>
    filter(! cdft)
# Energy derivatives where no field was applied. I expect these to be zero
nofield_derivatives <- energy_derivatives |>
    filter(! cdft)

this_theme <- 
    theme(
        # x axis text is too crowded, rotate it
        axis.text.x = element_text(angle = 90),
        legend.position = 'bottom'
    )

# Axis label for charge of the acceptor group, which I use in multiple plots
acceptor_charge_label <- xlab('Charge of electron acceptor group')

# Make a plot judging the quality of fit of the smoothed energy function to the
# original function
energy_smoothing_validation_plot <- smoothed_energy |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(mapping = aes(x = charge, y = total_energy)) +
    facet_wrap(~ formula, scales = 'free', ncol = 3) +
    geom_line() +
    geom_point(data =
               filter(charge_energy, donor_or_acceptor == 'acceptor'),
               size = 0, color = 'red') +
    acceptor_charge_label +
    this_theme
ggsave('energy_smoothing_validation.png', energy_smoothing_validation_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Plotting the smoothed energy, with points to indicate the no-field simulations
energy_with_nofield <- smoothed_energy |>
    left_join(nofield_energies, by = c('formula', 'donor_or_acceptor'),
              suffix = c('', '_nofield'), relationship = 'many-to-one') |>
    mutate(energy_above_nofield = total_energy - total_energy_nofield) |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(mapping = aes(x = charge, y = energy_above_nofield)) +
    facet_wrap(~ formula, scales = 'free', ncol = 3) +
    geom_smooth(method = lmrob, formula = y ~ x + I(x^2), se = FALSE) +
    geom_line() +
    geom_point(data = filter(mutate(nofield_energies, energy_above_nofield = 0), donor_or_acceptor == 'acceptor')) +
    acceptor_charge_label +
    this_theme +
    ylab('Energy above no field (eV)')

ggsave('energy_with_nofield.png', energy_with_nofield, width = unit(11.5, 'in'), height = unit(4.76, 'in'))
computation_levels <- c('dE/dq', '2 * Lam')

energy_derivatives_with_nofield <- smoothed_energy_derivative |>
    filter(donor_or_acceptor == 'acceptor') |>
    mutate(computation = factor('dE/dq', levels = computation_levels)) |>
    ggplot(aes(x = charge, y = total_energy, color = computation)) +
    facet_wrap(vars(formula), scales = 'free', ncol = 3) +
    geom_line() +
#    geom_point(data = filter(nofield_energies, donor_or_acceptor == 'acceptor')) +
    acceptor_charge_label +
    ylab('electronegativity difference (V)') +
    # Put a vertical line to indicate 0
    geom_vline(xintercept = 0, linetype = 'dashed') +
    # Put a horizontal line to indicate zero
    geom_hline(yintercept = 0, linetype = 'dashed')+
    this_theme +
    guides(color = guide_legend(title = NULL))
ggsave('energy_derivatives_with_nofield.png', energy_derivatives_with_nofield, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Zoom in to look for linearity near neutral molecule
energy_derivatives_zoomed <- smoothed_energy_derivative |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(aes(x = charge, y = total_energy)) +
    facet_wrap(vars(formula), ncol = 3) +
    # Add a horizontal line so we can see how far the no-field condition is
    # from zero
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # Mark equalized charge with a vertical line
    geom_vline(xintercept = 0, linetype = 'dashed') +
    # Add a linear fit to check how close to linear they are
    geom_smooth(method = lmrob, formula = y ~ x, se = FALSE) +
    geom_line() +
    geom_point(data = filter(nofield_energies, donor_or_acceptor == 'acceptor')) +
    ylim(c(-5,5)) +
    acceptor_charge_label +
    this_theme
ggsave('energy_derivatives_zoomed.png', energy_derivatives_zoomed, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

lam_values <- energy_derivatives |>
    filter(donor_or_acceptor == 'acceptor') |>
    mutate(`2 * Lam` = 2 * field_value,
           `dE/dq` = total_energy) |>
    pivot_longer(c(`2 * Lam`, `dE/dq`), names_to = 'computation', values_to = 'electronegativity')
lam_plot <- lam_values |>
    ggplot(aes(x = charge, y = electronegativity, color = computation)) +
    facet_wrap(~ formula, ncol = 3) +
    # Put a vertical line to indicate 0
    geom_vline(xintercept = 0, linetype = 'dashed') +
    # Put a horizontal line to indicate zero
    geom_hline(yintercept = 0, linetype = 'dashed')+
    xlab('charge of electron acceptor') +
    ylab('electronegativity difference (V)') +
    # Remove the title from the legend
    guides(color = guide_legend(title = NULL)) +
    geom_line() +
    acceptor_charge_label +
    this_theme
ggsave('lam_comparison.png', lam_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))
