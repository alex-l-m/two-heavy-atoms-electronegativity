library(tidyverse)
library(cowplot)
library(robustbase)
library(glue)
library(ggrepel)
library(tune)

this_theme <- 
    theme(
        # Usually x axis text is too crowded
        # But it depends on the font size I'm using and rotating loses vertical space
        axis.text.x = element_text(angle = 90),
        plot.background = element_rect(fill = 'white')
    )

theme_set(theme_cowplot(font_size = 24) + this_theme)

atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_double()
))
# Ordered vector to use as factor levels
element_levels <- atomic_numbers |>
    arrange(atomic_number) |>
    pull(symbol)

pauling_electronegativity <- read_csv('pauling_electronegativity.csv', col_types = cols(
    symbol = col_character(),
    pauling_electronegativity = col_double()
)) |>
    mutate(symbol = factor(symbol, levels = element_levels))

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
    electronegativity_field_analytic = col_double()
)) |>
    mutate(
        symbol = factor(symbol, levels = element_levels),
        other_symbol = factor(other_symbol, levels = element_levels)
    )

unique_scale_numbers <- charge_energy |>
    distinct(scale_number)

scale_numbers <- charge_energy |>
    select(combination_id, scale_number) |>
    distinct()

elements <- charge_energy |>
    select(combination_id, symbol, other_symbol) |>
    distinct()

# All formulas, with "symbol" as the anion and "other_symbol" as the cation
formula_symbols <- charge_energy |>
    filter(donor_or_acceptor == 'acceptor') |>
    distinct(formula, symbol, other_symbol)

# Charge of the acceptor atoms in ground state structures
acceptor_charge <- charge_energy |>
    filter(donor_or_acceptor == 'acceptor' & scale_number == 0 &
           crystal_structure == 'zincblende' & !cdft) |>
    select(formula, charge) |>
    rename(acceptor_charge = charge)

category_structure_pairs <- charge_energy |>
    mutate(category_structure_pair = glue('{category}:{crystal_structure}')) |>
    distinct(category_structure_pair) |>
    pull(category_structure_pair)
for (category_structure_pair in category_structure_pairs)
{
    loo <- read_csv(glue('{category_structure_pair}_loo.csv.gz'),
                    col_types = cols(
        formula = col_character(),
        electronegativity_donor = col_double(),
        electronegativity_acceptor = col_double(),
        hardness_donor = col_double(),
        hardness_acceptor = col_double(),
        eeq_acceptor_charge = col_double()
    ))
    # Lines indicating the regression fits from leave one out
    # y intercept is the difference between electronegativities, slope is the
    # sum of the hardnesses
    loo_lines <- loo |>
        select(formula, electronegativity_donor, electronegativity_acceptor,
               hardness_donor, hardness_acceptor) |>
        mutate(
            # The slope is the sum of the hardnesses
            slope = hardness_donor + hardness_acceptor,
            # The y intercept is the difference between the electronegativities
            y_intercept = electronegativity_acceptor - electronegativity_donor) |>
        select(-electronegativity_donor, -electronegativity_acceptor,
               -hardness_donor, -hardness_acceptor) |>
        left_join(formula_symbols, by = 'formula',
                  relationship = 'one-to-one')
    loo_charge_pred <- loo |>
        select(formula, eeq_acceptor_charge)
    loo_charge_comparison_table <- loo_charge_pred |>
        left_join(acceptor_charge, by = 'formula')

    loo_charge_comparison_plot <- loo_charge_comparison_table |>
        ggplot(aes(x = acceptor_charge, y = eeq_acceptor_charge,
                   label = formula)) +
        geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
        geom_point() +
        geom_label_repel() +
        coord_obs_pred()
    ggsave(glue('{category_structure_pair}_loo_charge_comparison_plot.png'),
           height = unit(4.76, 'in'), width = unit(5.67, 'in'))
    write_rds(loo_charge_comparison_plot,
              glue('{category_structure_pair}_loo_charge_comparison_plot.rds'))

    ranked_elements <-
        read_csv(glue('{category_structure_pair}_ranked_elements.csv.gz'), col_types = cols(
            symbol = col_character(),
            atomic_number = col_double(),
            element_rank = col_double()
        )) |>
        mutate(symbol = factor(symbol, levels = element_levels))

    # Read the previously saved parameter estimates from regression
    estimates <- read_csv(glue('{category_structure_pair}_regression_estimates.csv.gz'),
             col_types = cols(
                 term = col_character(),
                 estimate = col_double(),
                 std.error = col_double(),
                 statistic = col_double(),
                 p.value = col_double()
             )) |>
        select(term, estimate)

    electronegativity_comparison_tbl <- ranked_elements |>
        select(symbol) |>
        left_join(pauling_electronegativity, by = 'symbol',
                  relationship = 'one-to-one') |>
        cross_join(unique_scale_numbers) |>
        # Create a column with the name of the term corresponding to the
        # regression electronegativity of the element
        mutate(term = glue('electronegativity_{symbol}_S{scale_number}')) |>
        # Retrieve the regression electronegativity by joining with the
        # parameter table
        left_join(estimates, by = 'term', relationship = 'one-to-one') |>
        # If a symbol is missing that means there wasn't a parameter estimate
        # for it which means I implicitly set it to zero by leaving it out of
        # the regression model. So replace missing values with zero
        mutate(regression_electronegativity = replace_na(estimate, 0)) |>
        select(-term, -estimate)

    electronegativity_comparison_plot <- electronegativity_comparison_tbl |>
        ggplot(aes(x = pauling_electronegativity, y = regression_electronegativity,
                   label = symbol)) +
        # Divide it up by scale
        facet_wrap(~ scale_number) +
        geom_point() +
        geom_smooth(method = lmrob, se = FALSE) +
        geom_label_repel() +
        ylab('Regression electronegativity (V)') +
        xlab('Pauling electronegativity')

    ggsave(glue('{category_structure_pair}_electronegativity_comparison_plot.png'),
           electronegativity_comparison_plot,
           height = unit(4.76, 'in'), width = unit(5.67, 'in'))

    # Same plot, but for just the equilibrium geometry
    electronegativity_comparison_plot_eq <- electronegativity_comparison_tbl |>
        filter(scale_number == 0) |>
        ggplot(aes(x = pauling_electronegativity, y = regression_electronegativity,
                   label = symbol)) +
        geom_point() +
        geom_smooth(method = lmrob, se = FALSE) +
        geom_label_repel() +
        ylab('Regression electronegativity (V)') +
        xlab('Pauling electronegativity')
    
    electronegativity_comparison_base <- glue('{category_structure_pair}_electronegativity_comparison_plot_eq')
    ggsave(glue('{electronegativity_comparison_base}.png'),
           electronegativity_comparison_plot_eq,
           height = unit(4.76, 'in'), width = unit(5.67, 'in'))
    write_rds(electronegativity_comparison_plot_eq,
             glue('{electronegativity_comparison_base}.rds'))
    
    hardness_terms <- read_csv(glue('{category_structure_pair}_hardness_terms.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            category = col_character(),
            variable_contribution = col_double()
        ))
    # Read the table of electronegativity differences
    electronegativity_differences <- read_csv(glue('{category_structure_pair}_electronegativity_differences.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            variable_contribution = col_double()
        ))
    regression_plot_table <- charge_energy |>
        select(combination_id, symbol, other_symbol, charge,
               scale_number, structure_id, donor_or_acceptor) |>
        filter(donor_or_acceptor == 'acceptor') |>
        left_join(rename(electronegativity_differences,
                         electronegativity_difference = variable_contribution),
                  by = 'combination_id', relationship = 'one-to-one')

    hardness_regression_plot <- regression_plot_table |>
        # Filter so that I'm only making plots where "symbol" corresponds to
        # the symbol of the acceptor
        filter(donor_or_acceptor == 'acceptor') |>
        ggplot(aes(x = charge, y = electronegativity_difference, color = other_symbol, group = structure_id)) +
        facet_wrap(~ symbol) +
        geom_hline(yintercept = 0, linetype = 'dotted') +
        geom_vline(xintercept = 0, linetype = 'dotted') +
        geom_line() +
        # I want this plot to also include the theoretical lines
        geom_abline(data = loo_lines,
                    mapping = aes(slope = slope, intercept = y_intercept,
                                  color = other_symbol),
                    linetype = 'dashed') +
        scale_x_continuous(breaks = seq(-1, 1, 0.5), name = 'Charge of acceptor group') +
        scale_y_continuous(name = 'Δelectronegativity (V)') +
        theme(legend.position = 'bottom',
              # Settings for facet panels so they're sufficiently well spaced
              # https://stackoverflow.com/questions/28652284/how-to-change-color-of-facet-borders-when-using-facet-grid
              panel.spacing = unit(.1, 'lines'),
              panel.border = element_rect(color = 'grey', fill = NA, size = 1)
        )
    hardness_regression_base <- glue('{category_structure_pair}_hardness_regression_plot')
    ggsave(glue('{hardness_regression_base}.png'), hardness_regression_plot, height = unit(4.76, 'in'), width = unit(11.5, 'in'))
    write_rds(hardness_regression_plot, glue('{hardness_regression_base}.rds'))
}
