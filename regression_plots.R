library(tidyverse)
library(cowplot)
library(robustbase)
library(glue)
library(ggrepel)

this_theme <- 
    theme(
        # Usually x axis text is too crowded
        # But it depends on the font size I'm using and rotating loses vertical space
        axis.text.x = element_text(angle = 90),
        plot.background = element_rect(fill = 'white')
    )

theme_set(theme_cowplot(font_size = 24) + this_theme)

pauling_electronegativity <- read_csv('pauling_electronegativity.csv', col_types = cols(
    symbol = col_character(),
    pauling_electronegativity = col_double()
))

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
))

unique_scale_numbers <- charge_energy |>
    distinct(scale_number)

scale_numbers <- charge_energy |>
    select(combination_id, scale_number) |>
    distinct()

elements <- charge_energy |>
    select(combination_id, symbol, other_symbol) |>
    distinct()

category_structure_pairs <- charge_energy |>
    mutate(category_structure_pair = glue('{category}:{crystal_structure}')) |>
    distinct(category_structure_pair) |>
    pull(category_structure_pair)
for (category_structure_pair in category_structure_pairs)
{
    ranked_elements <-
        read_csv(glue('{category_structure_pair}_ranked_elements.csv.gz'), col_types = cols(
            symbol = col_character(),
            atomic_number = col_double(),
            element_rank = col_double()
        ))
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
    regression_plot_table <- hardness_terms |>
        mutate(
               # Extract the symbol from the category
               symbol = str_match(category, '([^_]*)_S\\d+')[, 2],
               this_charge = variable_contribution) |>
        select(-category, -variable_contribution) |>
        left_join(elements, by = c('combination_id', 'symbol'), relationship = 'many-to-one') |>
        # Join the scale numbers. I originally joined this because I was
        # intending to use it as a grouping variable. It makes more sense, I
        # realized, to just use the structure id. But there's no harm in having
        # an extra variable in the table
        left_join(scale_numbers, by = 'combination_id',
                  relationship = 'many-to-one') |>
        left_join(rename(electronegativity_differences,
                         electronegativity_difference = variable_contribution),
                  by = 'combination_id', relationship = 'many-to-one') |>
        # Restore the information about structure id, for grouping on the plot
        left_join(distinct(charge_energy, combination_id, structure_id),
                  by = 'combination_id', relationship = 'many-to-one') |>
        # Restore information about which is donor and which is acceptor
        # This won't work for 4-4's
        left_join(distinct(charge_energy, combination_id, symbol,
                           donor_or_acceptor),
                  by = c('combination_id', 'symbol'),
                  relationship = 'many-to-one')

    # Make plots showing the data underlying the regression
    # Make a vector of elements ordered by their rank to use as factor levels
    ordered_selected_elements <- ranked_elements |>
        arrange(element_rank) |>
        pull(symbol)
    hardness_regression_plot <- regression_plot_table |>
        # Filter so that I'm only making plots where "symbol" corresponds to
        # the symbol of the acceptor
        filter(donor_or_acceptor == 'acceptor') |>
        # Turn the symbols into factors based on the ordering implied by the
        # atomic numbers
        mutate(symbol = factor(symbol, levels = ordered_selected_elements),
               other_symbol = factor(other_symbol, levels = ordered_selected_elements)) |>
        ggplot(aes(x = this_charge, y = electronegativity_difference, color = other_symbol, group = structure_id)) +
        facet_wrap(~ symbol) +
        geom_hline(yintercept = 0, linetype = 'dotted') +
        geom_vline(xintercept = 0, linetype = 'dotted') +
        geom_line() +
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
