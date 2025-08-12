library(tidyverse)
library(cowplot)
library(robustbase)
theme_set(theme_cowplot() + theme(plot.background = element_rect(fill = 'white')))
library(glue)
library(ggrepel)
library(broom)

TO_EV <- 27.211386246

# Atomic numbers, for ordering
atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_integer()
))

# Smoothed energies and derivatives
smoothed_energy <- read_csv('smoothed_energy.csv.gz', col_types = cols(
    crystal_structure = col_character(),
    structure_id = col_character(),
    category = col_character(),
    symbol_anion = col_character(),
    symbol_cation = col_character(),
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
    crystal_structure = col_character(),
    structure_id = col_character(),
    category = col_character(),
    symbol_anion = col_character(),
    symbol_cation = col_character(),
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
    bader_charge = col_double(),
    cp2k_hirshfeld_charge = col_double(),
    weight_function_derivative = col_double(),
    structure_id = col_character(),
    field_number = col_integer(),
    field_value = col_double(),
    total_energy = col_double(),
    cdft = col_logical(),
    group_id = col_character(),
    category = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character(),
    crystal_structure = col_character(),
    scale_number = col_integer(),
    scale = col_double(),
    unscaled_structure_id = col_character(),
    formula = col_character(),
    first_iteration_charge = col_double(),
    last_iteration_charge = col_double(),
    lagged_first_iteration_charge = col_double(),
    lagged_last_iteration_charge = col_double(),
    electronegativity_field_discrete = col_double(),
    electronegativity_field_analytic = col_double(),
    unscaled_structure_id = col_character(),
    scale_number = col_integer(),
    scale = col_double()
))

energy_derivatives <- read_csv('energy_derivatives.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    group_id = col_character(),
    charge = col_double(),
    crystal_structure = col_character(),
    structure_id = col_character(),
    category = col_character(),
    symbol_anion = col_character(),
    symbol_cation = col_character(),
    formula = col_character(),
    cdft = col_logical(),
    field_number = col_integer(),
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

# Make separate plots for each category and crystal structure
# Also the scale, although I'm not going to bother changing the variable name
# to reflect that. So these are really category structure scale triples
# This requires changing a shocking number of lines, so if I add anything else
# to that I should refactor somehow
category_structure_pairs <- charge_energy |>
    mutate(category_structure_pair = glue('{category}:{crystal_structure}:{scale_number}')) |>
    distinct(category_structure_pair) |>
    pull(category_structure_pair)
for (category_structure_pair in category_structure_pairs)
{
    # Make a plot judging the quality of fit of the smoothed energy function to the
    # original function
    formula_order <- charge_energy |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
        select(formula, symbol_cation, symbol_anion) |>
        # Order based on atomic numbers of the elements
        left_join(atomic_numbers, by = c('symbol_cation' = 'symbol')) |>
        left_join(atomic_numbers, by = c('symbol_anion' = 'symbol'),
                                  suffix = c('_cation', '_anion')) |>
        arrange(atomic_number_cation, atomic_number_anion) |>
        distinct(formula) |>
        pull(formula)
    energy_smoothing_validation_tbl <- smoothed_energy |>
        filter(donor_or_acceptor == 'acceptor') |>
        mutate(formula = factor(formula, levels = formula_order)) |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair)

    energy_smoothing_validation_plot <- energy_smoothing_validation_tbl |>
        ggplot(mapping = aes(x = charge, y = total_energy)) +
        facet_wrap(~ formula, scales = 'free', ncol = 3) +
        geom_line() +
        geom_point(data =
                   filter(charge_energy, donor_or_acceptor == 'acceptor', glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair),
                   size = 0, color = 'red') +
        acceptor_charge_label +
        this_theme
    ggsave(glue('{category_structure_pair}_energy_smoothing_validation.png'), energy_smoothing_validation_plot,
           width = unit(11.5, 'in'), height = unit(4.76, 'in'))

    # Plotting the smoothed energy, with points to indicate the no-field simulations
    energy_with_nofield_tbl <- smoothed_energy |>
        left_join(nofield_energies, by =
                      c('crystal_structure', 'category', 'structure_id',
                        'symbol_anion', 'symbol_cation', 'formula',
                        'symbol', 'other_symbol', 'donor_or_acceptor', 'group_id'),
                  suffix = c('', '_nofield'), relationship = 'many-to-one') |>
        mutate(energy_above_nofield = total_energy - total_energy_nofield) |>
        filter(donor_or_acceptor == 'acceptor') |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair)
    energy_with_nofield_facet_plot <- energy_with_nofield_tbl |>
        ggplot(mapping = aes(x = charge, y = energy_above_nofield)) +
        facet_wrap(~ formula, scales = 'free', ncol = 3) +
        geom_smooth(mapping = aes(color = 'Quadratic fit'),
                    lwd = 2, lty = 'dashed',
                    formula = y ~ x + I(x^2),
                    se = FALSE, method = lm) +
        # Don't show the name of the guide in the color legend (that is,
        # display it, but without the label "color")
        scale_color_discrete(guide = guide_legend(title = NULL)) +
        geom_line(mapping = aes(color = 'Energy')) +
        geom_point(data = filter(mutate(nofield_energies, energy_above_nofield = 0), donor_or_acceptor == 'acceptor' & glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair)) +
        acceptor_charge_label +
        this_theme +
        ylab('Energy above no field (eV)')

    ggsave(glue('{category_structure_pair}_energy_with_nofield.png'), energy_with_nofield_facet_plot, width = unit(11.5, 'in'), height = unit(4.76, 'in'))
    computation_levels <- c('dE/dq', 'correction * Lam')

    # Same plots, but as a list of individual plots rather than as a facet
    energy_with_nofield_plot_tbl <- energy_with_nofield_tbl |>
        group_by(formula) |>
        nest() |>
        mutate(plot = lapply(data, function(data) {
        data |>
            ggplot(mapping = aes(x = charge, y = energy_above_nofield)) +
            geom_smooth(mapping = aes(color = 'Quadratic fit'),
                        lwd = 2, lty = 'dashed',
                        formula = y ~ x + I(x^2),
                        se = FALSE, method = lm) +
            geom_line(mapping = aes(color = 'Energy')) +
            # Don't show the name of the guide in the color legend (that is,
            # display it, but without the label "color")
            scale_color_discrete(guide = guide_legend(title = NULL)) +
            acceptor_charge_label +
            this_theme +
            ylab('Energy above no field (eV)')
        }))
    energy_with_nofield_plot_list <- energy_with_nofield_plot_tbl$plot
    names(energy_with_nofield_plot_list) <- energy_with_nofield_plot_tbl$formula
    # Save as rds using readr
    readr::write_rds(energy_with_nofield_plot_list, glue('{category_structure_pair}_energy_with_nofield_plots.rds'))

    
    energy_derivatives_with_nofield <- smoothed_energy_derivative |>
        filter(donor_or_acceptor == 'acceptor') |>
        mutate(computation = factor('dE/dq', levels = computation_levels)) |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
        ggplot(aes(x = charge, y = total_energy, color = computation)) +
        facet_wrap(vars(formula), scales = 'free', ncol = 3) +
        geom_line() +
        geom_point(data = mutate(filter(nofield_derivatives,
                                        donor_or_acceptor == 'acceptor' & glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair), computation = factor('dE/dq', levels = computation_levels))) +
        acceptor_charge_label +
        ylab('Δelectronegativity (V)') +
        # Put a vertical line to indicate 0
        geom_vline(xintercept = 0, linetype = 'dashed') +
        # Put a horizontal line to indicate zero
        geom_hline(yintercept = 0, linetype = 'dashed')+
        this_theme +
        guides(color = guide_legend(title = NULL))
    ggsave(glue('{category_structure_pair}_energy_derivatives_with_nofield.png'), energy_derivatives_with_nofield, width = unit(11.5, 'in'), height = unit(4.76, 'in'))
    
    # Zoom in to look for linearity near neutral molecule
    energy_derivatives_zoomed <- smoothed_energy_derivative |>
        filter(donor_or_acceptor == 'acceptor') |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
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
        geom_point(data = mutate(filter(nofield_derivatives,
                                        donor_or_acceptor == 'acceptor' & glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair), computation = factor('dE/dq', levels = computation_levels))) +
        ylim(c(-5,5)) +
        acceptor_charge_label +
        this_theme
    ggsave(glue('{category_structure_pair}_energy_derivatives_zoomed.png'), energy_derivatives_zoomed, width = unit(11.5, 'in'), height = unit(4.76, 'in'))

    # Verify the analytic and discrete electronegativities come out the same
    electronegativity_differentiation_check <- charge_energy |>
        filter(donor_or_acceptor == 'acceptor') |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
        ggplot(aes(x = electronegativity_field_discrete, y = electronegativity_field_analytic)) +
        facet_wrap(vars(formula), ncol = 3) +
        geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
        geom_point() +
        xlab('Discrete electronegativity') +
        ylab('Analytic electronegativity') +
        this_theme
    ggsave(glue('{category_structure_pair}_differentiation_check.png'),
           electronegativity_differentiation_check,
           width = unit(11.5, 'in'), height = unit(4.76, 'in'))
    
    electronegativity_difference_from_field <- charge_energy |>
        select(combination_id, field_number, donor_or_acceptor, electronegativity_field_analytic) |>
        pivot_wider(names_from = donor_or_acceptor, values_from = electronegativity_field_analytic) |>
        group_by(combination_id, field_number) |>
        transmute(electronegativity_difference_from_field = acceptor + donor) |>
        ungroup()
    lam_values <- energy_derivatives |>
        filter(donor_or_acceptor == 'acceptor') |>
        # Join correction factors to correct the potential
        left_join(electronegativity_difference_from_field,
                  by = c('combination_id', 'field_number')) |>
        mutate(`correction * Lam` = electronegativity_difference_from_field,
               `dE/dq` = total_energy) |>
        pivot_longer(c(`correction * Lam`, `dE/dq`), names_to = 'computation', values_to = 'electronegativity')
    lam_facet_plot_tbl <- lam_values |>
        mutate(formula = factor(formula, levels = formula_order)) |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair)
    lam_facet_plot <- lam_facet_plot_tbl |>
        ggplot(aes(x = charge, y = electronegativity, color = computation)) +
        facet_wrap(~ formula, ncol = 3) +
        # Put a vertical line to indicate 0
        geom_vline(xintercept = 0, linetype = 'dashed') +
        # Put a horizontal line to indicate zero
        geom_hline(yintercept = 0, linetype = 'dashed')+
        ylab('Δelectronegativity (V)') +
        # Remove the title from the legend
        guides(color = guide_legend(title = NULL)) +
        geom_line() +
        acceptor_charge_label +
        this_theme
    ggsave(glue('{category_structure_pair}_lam_comparison.png'), lam_facet_plot,
           width = unit(11.5, 'in'), height = unit(4.76, 'in'))
    
    # Same thing but a list of plots instead of facets
    lam_plot_tbl <- lam_facet_plot_tbl |>
        group_by(formula) |>
        nest() |>
        mutate(plot = lapply(data, function(data) {
            data |>
                ggplot(aes(x = charge, y = electronegativity, color = computation)) +
                # Put a vertical line to indicate 0
                geom_vline(xintercept = 0, linetype = 'dashed') +
                # Put a horizontal line to indicate zero
                geom_hline(yintercept = 0, linetype = 'dashed')+
                ylab('Δelectronegativity (V)') +
                # Remove the title from the legend
                guides(color = guide_legend(title = NULL)) +
                geom_line() +
                acceptor_charge_label
        }))
    lam_plot_list <- lam_plot_tbl$plot
    names(lam_plot_list) <- lam_plot_tbl$formula
    # Save as rds using readr
    readr::write_rds(lam_plot_list, glue('{category_structure_pair}_lam_plots.rds'))

    # Make a plot like the lam comparison plot but instead of comparing lam to
    # derivative, show just derivative, but with various kinds of charge
    charge_comparison_tbl <- energy_derivatives |>
        inner_join(select(charge_energy,
                         # Identify the simulation
                         combination_id, scale_number,
                         # Identify the atom
                         donor_or_acceptor,
                         # The charge columns that I want
                         cp2k_hirshfeld_charge, bader_charge),
                  by = c('combination_id', 'scale_number', 'donor_or_acceptor'),
                  relationship = 'one-to-one') |>
        rename(electronegativity = total_energy) |>
        select(combination_id, scale_number,
               category, crystal_structure, formula,
               donor_or_acceptor, symbol, other_symbol,
               electronegativity,
               charge, cp2k_hirshfeld_charge, bader_charge) |>
        # Pivot to create a single charge column so that I can compare kinds of
        # charge
        pivot_longer(c(charge, cp2k_hirshfeld_charge, bader_charge), names_to = 'charge_definition', values_to = 'charge') |>
        filter(donor_or_acceptor == 'acceptor')
    charge_comparison_plot <- charge_comparison_tbl |>
        mutate(formula = factor(formula, levels = formula_order)) |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
        ggplot(aes(x = charge, y = electronegativity, color = charge_definition)) +
        facet_wrap(~ formula, ncol = 3) +
        # Put a vertical line to indicate 0
        geom_vline(xintercept = 0, linetype = 'dashed') +
        # Put a horizontal line to indicate zero
        geom_hline(yintercept = 0, linetype = 'dashed')+
        ylab('dE/dq (V)') +
        # Remove the title from the legend
        guides(color = guide_legend(title = NULL)) +
        geom_line() +
        acceptor_charge_label +
        this_theme
    ggsave(glue('{category_structure_pair}_charge_comparison_electronegativity.png'), charge_comparison_plot,
           width = unit(11.5, 'in'), height = unit(4.76, 'in'))
    
    # Fit lines to the lambda values, from the x axis and y axis
    lam_line_data <- lam_values |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
        filter(computation == 'correction * Lam') |>
        # Compute rankings, so I can select the bottom and top charges from each
        # material
        group_by(formula) |>
        mutate(bottom_charge_rank = rank(charge, ties.method = 'min'),
               top_charge_rank = rank(-charge, ties.method = 'min')) |>
        ungroup() |>
        # Label as belonging to the line from the ground state or the line from the
        # neutral state
        mutate(line = if_else(bottom_charge_rank <= 5, 'ground',
                              if_else(top_charge_rank <= 5, 'neutral',
                                      'neither'))) |>
        # Filter out observations belonging to neither
        filter(line != 'neither') |>
        # Compute slope and intercept of the lines
        group_by(formula, line) |>
        summarize(slope = coef(lm(electronegativity ~ charge))['charge'],
                  intercept = coef(lm(electronegativity ~ charge))['(Intercept)'],
                  .groups = 'drop')
    
    lam_lines_plot <- lam_values |>
        filter(computation == 'correction * Lam') |>
        filter(glue('{category}:{crystal_structure}:{scale_number}') == category_structure_pair) |>
        ggplot(aes(x = charge, y = electronegativity, color = computation)) +
        facet_wrap(~ formula, ncol = 3) +
        # Put a vertical line to indicate 0
        geom_vline(xintercept = 0, linetype = 'dashed') +
        # Put a horizontal line to indicate zero
        geom_hline(yintercept = 0, linetype = 'dashed')+
        ylab('Δelectronegativity (V)') +
        # Remove the title from the legend
        guides(color = guide_legend(title = NULL)) +
        geom_line() +
        # Don't need to filter lam_line_data for the right category and crystal
        # structure because that was already done during creation
        geom_abline(data = lam_line_data,
                    aes(intercept = intercept, slope = slope, color = line)) +
        acceptor_charge_label +
        this_theme +
        theme(legend.position = 'bottom')
    
    ggsave(glue('{category_structure_pair}_lam_lines_comparison.png'), lam_lines_plot,
           width = unit(11.5, 'in'), height = unit(4.76, 'in'))
}
