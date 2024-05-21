library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
theme_set(dark_mode(theme_cowplot()))
library(glue)
library(ggrepel)

# Smoothed energies and derivatives
smoothed_energy <- read_csv('smoothed_energy.csv.gz', col_types = cols(
    formula = col_character(),
    molecule_charge = col_double(),
    symbol = col_character(),
    donor_or_acceptor = col_character(),
    group_bader_charge = col_double(),
    iqa_group_energy = col_double(),
    iqa_interaction_energy = col_double(),
    total_energy = col_double(),
    other_symbol = col_character(),
    group_id = col_character()
))

smoothed_energy_derivative <- read_csv('smoothed_energy_derivatives.csv.gz',
                                       col_types = cols(
    formula = col_character(),
    molecule_charge = col_double(),
    symbol = col_character(),
    donor_or_acceptor = col_character(),
    group_bader_charge = col_double(),
    iqa_group_energy = col_double(),
    iqa_interaction_energy = col_double(),
    total_energy = col_double(),
    other_symbol = col_character(),
    group_id = col_character()
))

# Raw energies and derivatives
charge_energy <- read_csv('charge_energy.csv.gz', col_types = cols(
    formula = col_character(),
    symbol = col_character(),
    donor_or_acceptor = col_character(),
    combination_id = col_character(),
    cdft = col_logical(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer(),
    group_bader_charge = col_double(),
    total_energy = col_double(),
    iqa_group_energy = col_double(),
    iqa_interaction_energy = col_double(),
    other_symbol = col_character(),
    group_id = col_character()
))

energy_derivatives <- read_csv('energy_derivatives.csv.gz', col_types = cols(
    formula = col_character(),
    molecule_charge = col_double(),
    symbol = col_character(),
    donor_or_acceptor = col_character(),
    combination_id = col_character(),
    cdft = col_logical(),
    field_number = col_double(),
    field_value = col_double(),
    group_bader_charge = col_double(),
    iqa_group_energy = col_double(),
    iqa_interaction_energy = col_double(),
    total_energy = col_double()
))

# Energies for the simulations where no field was applied, so I can add them as
# a dot on the plots. I expect them to be at the minimum of the parabola
nofield_energies <- charge_energy |>
    filter(! cdft & molecule_charge == 0)
# Energy derivatives where no field was applied. I expect these to be zero
nofield_derivatives <- energy_derivatives |>
    filter(! cdft & molecule_charge == 0)

this_theme <- 
    theme(
        # x axis text is too crowded, rotate it
        axis.text.x = element_text(angle = 90)
    )

# Axis label for charge of the acceptor group, which I use in multiple plots
acceptor_charge_label <- xlab('Bader charge of electron acceptor group')

# Make a plot judging the quality of fit of the smoothed energy function to the
# original function
energy_smoothing_validation_plot <- smoothed_energy |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(mapping = aes(x = group_bader_charge, y = total_energy)) +
    facet_wrap(~ formula, scales = 'free', nrow = 3) +
    geom_line() +
    geom_point(data =
               filter(charge_energy, cdft, donor_or_acceptor == 'acceptor'),
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
    ggplot(mapping = aes(x = group_bader_charge, y = energy_above_nofield)) +
    facet_wrap(~ formula, scales = 'free', nrow = 3) +
    geom_smooth(method = lmrob, formula = y ~ x + I(x^2), se = FALSE) +
    geom_line() +
    geom_point(data = filter(mutate(nofield_energies, energy_above_nofield = 0), donor_or_acceptor == 'acceptor')) +
    acceptor_charge_label +
    this_theme
ggsave('energy_with_nofield.png', energy_with_nofield, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

energy_derivatives_with_nofield <- smoothed_energy_derivative |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(aes(x = group_bader_charge, y = total_energy)) +
    facet_wrap(vars(formula), scales = 'free', nrow = 3) +
    geom_line() +
#    geom_point(data = filter(nofield_energies, donor_or_acceptor == 'acceptor')) +
    acceptor_charge_label +
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
    filter(donor_or_acceptor == 'acceptor') |>
    mutate(
        intercept_1 = total_energy - slope_1 * group_bader_charge,
        intercept_2 = total_energy - slope_2 * group_bader_charge
    )

energy_derivatives_with_nofield_lines <- energy_derivatives_with_nofield +
    geom_abline(aes(slope = slope_1, intercept = intercept_1), data = theoretical_lines, color = 'orange') +
    geom_abline(aes(slope = slope_2, intercept = intercept_2), data = theoretical_lines, color = 'orange')

ggsave('energy_derivatives_with_nofield_lines.png', energy_derivatives_with_nofield_lines, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Zoom in to look for linearity near neutral molecule
energy_derivatives_zoomed <- smoothed_energy_derivative |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(aes(x = group_bader_charge, y = total_energy)) +
    facet_wrap(vars(formula), nrow = 3) +
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

# Read the 'Lam' values
lam <- read_csv('lam_values.csv.gz', col_types = cols(
    combination_id = col_character(),
    lam = col_double()
))
# Doesn't have the formula column, or the charge transfer. Join with the table
# containing energies and charges, just for the charges
# This is a messy table, an ugly hack. Really lam should be added during
# smoothing... or something, I dont' know
lam_charges <- lam |>
    inner_join(charge_energy,
               by = 'combination_id')

# Considering it as estimates of the derivative
lam_derivs <- lam_charges |>
    # Just putting the negative sign because that makes it match, I don't get
    # why. Factor of two is because my linear combination has total weight two,
    # though I thought that means I need to divide by two, but multiplying
    # makes it match
    mutate(total_energy = ifelse(donor_or_acceptor == 'acceptor', -1, 1) * 2 * lam)

# Table containing both the lambda values and the energy derivatives, for
# comparison
# Each Lam value comes from a specific simulation, but each energy derivative
# comes from a smoothing process. So only Lam has a combination id, but this
# won't be relevant in the comparison
lam_comparison <- bind_rows(
        `2 * Lam` = select(lam_derivs, -combination_id),
        `dE/dq` = energy_derivatives,
        .id = 'computation') |>
    # Put Lam first, so it's the same color if I plot just that in another plot
    mutate(computation = factor(computation, levels = c('2 * Lam', 'dE/dq')))

lam_plot <- lam_comparison |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(aes(x = group_bader_charge, y = total_energy, color = computation)) +
    facet_wrap(~ formula, scales = 'free', nrow = 3) +
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
    acceptor_charge_label +
    this_theme
ggsave('lam_comparison.png', lam_plot,
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Plot of the IQA group energies
iqa_energy_plot <- smoothed_energy |>
    # Charges between -1 and 1
    filter(-1 < group_bader_charge & group_bader_charge < 1) |>
    # I don't see why this should be necessary
    filter(! is.na(other_symbol)) |>
    ggplot(aes(x = group_bader_charge, y = iqa_group_energy,
               group = group_id, color = other_symbol)) +
    facet_wrap(vars(symbol), scales = 'free_y', ncol = 4) +
    geom_line() +
#    geom_point(data = nofield_energies) +
    theme(axis.text.x = element_text(angle = 90))
ggsave('iqa_energy.png', iqa_energy_plot, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Plot of the IQA group energy derivatives
iqa_energy_derivative_plot <- smoothed_energy_derivative |>
    filter(-1 < group_bader_charge & group_bader_charge < 1) |>
    # I don't see why this should be necessary
    filter(! is.na(other_symbol)) |>
    ggplot(aes(x = group_bader_charge,
               y = iqa_group_energy,
               group = group_id, color = other_symbol)) +
    facet_wrap(vars(symbol), ncol = 4) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90)) +
    # Vertical line to indicate the location of zero charge
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_point(data = nofield_derivatives) +
    # No variable name label on the legend
    guides(color = guide_legend(title = NULL))
ggsave('iqa_energy_derivative.png', iqa_energy_derivative_plot, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

# Plot of the IQA interaction energies
interaction_plot <- smoothed_energy |>
    left_join(nofield_energies, by = c('formula', 'donor_or_acceptor'),
              suffix = c('', '_nofield'), relationship = 'many-to-one') |>
    group_by(symbol, other_symbol) |>
    mutate(energy_above_nofield =
           iqa_interaction_energy - iqa_interaction_energy_nofield) |>
    ungroup() |>
    filter(-1 < group_bader_charge & group_bader_charge < 1) |>
    ggplot(aes(x = group_bader_charge, y = energy_above_nofield)) +
    facet_grid(symbol ~ other_symbol) +
    geom_point(data = mutate(nofield_energies, energy_above_nofield = 0)) +
    geom_line()
ggsave('interaction_plot.png', interaction_plot, 
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))

