library(broom)
library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
library(glue)
library(ggrepel)

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
    structure_id = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    total_energy = col_double(),
    cdft = col_logical(),
    group_id = col_character(),
    category = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character(),
    crystal_structure = col_character(),
    formula = col_character()
))

elements <- charge_energy |>
    select(combination_id, symbol, other_symbol) |>
    distinct()

this_theme <- 
    theme(
        # Usually x axis text is too crowded
        # But it depends on the font size I'm using and rotating loses vertical space
        axis.text.x = element_text(angle = 90)
    )

theme_set(dark_mode(theme_cowplot(font_size = 24)) + this_theme)

category_structure_pairs <- charge_energy |>
    mutate(category_structure_pair = glue('{category}:{crystal_structure}')) |>
    distinct(category_structure_pair) |>
    pull(category_structure_pair)
for (category_structure_pair in category_structure_pairs)
{
    # I used to filter for only those where a charge was applied, and those
    # were computed with constrained DFT. I don't do that anymore but the name
    # remains
    cdft_charges <- charge_energy |>
        filter(glue('{category}:{crystal_structure}') == category_structure_pair) |>
        select(combination_id, formula, symbol, other_symbol, donor_or_acceptor,
               charge)
    
    electronegativity_terms <- cdft_charges |>
        mutate(category = symbol,
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1)) |>
        select(combination_id, category, variable_contribution)
    
    hardness_terms <- cdft_charges |>
        mutate(category = symbol,
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution)
    
    interaction_terms <- cdft_charges |>
        mutate(category = formula, 
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution)
    
    electronegativity_differences <- cdft_charges |>
        # Would be "conceptually" more accurate to subtract the values of the donor
        # and acceptor after negating one of them, but the result will be the same
        # as doubling
        filter(donor_or_acceptor == 'acceptor') |>
        mutate(variable_contribution = 2 * field_value) |>
        select(combination_id, variable_contribution)
    
    regression_terms <- bind_rows(
            electronegativity_difference = electronegativity_differences,
            electronegativity = electronegativity_terms,
            hardness = hardness_terms,
            interaction = interaction_terms,
            .id = 'variable_type') |>
        mutate(variable_name = ifelse(is.na(category), variable_type,
                                      glue('{variable_type}_{category}')))
    
    pre_regression_table <- regression_terms |>
        # Temporary: no interaction terms
        filter(variable_type != 'interaction') |>
        group_by(combination_id, variable_name) |>
        summarize(variable_value = sum(variable_contribution), .groups = 'drop') |>
        pivot_wider(names_from = variable_name, values_from = variable_value,
                    values_fill = 0)
    
    regression_table <- pre_regression_table |>
        # Remove electronegativity of nitrogen to use it as a reference
        select(-electronegativity_N) |>
        select(-combination_id) |>
        as.data.frame()
    rownames(regression_table) <- pre_regression_table$combination_id
    
    linear_model <- lmrob(electronegativity_difference ~ . + 0,
                          data = regression_table)
    estimates <- tidy(linear_model)
    write_csv(estimates, glue('{category_structure_pair}_electronegativity_regression_estimates.csv'))
    
    regression_electronegativity <- estimates |>
        mutate(
            term_type = str_match(term, '([^_]*)_([^_]*)')[,2],
            symbol = str_match(term, '([^_]*)_([^_]*)')[,3]
        ) |>
        filter(term_type == 'electronegativity') |>
        select(symbol, estimate) |>
        rename(regression_electronegativity = estimate) |>
        bind_rows(tibble(symbol = 'N', regression_electronegativity = 0))
    
    electronegativity_comparison_tbl <- left_join(regression_electronegativity,
                                                  pauling_electronegativity,
                                                  by = 'symbol')
    
    electronegativity_comparison_plot <- electronegativity_comparison_tbl |>
        ggplot(aes(x = pauling_electronegativity, y = regression_electronegativity,
                   label = symbol)) +
        geom_point() +
        geom_smooth(method = lmrob, se = FALSE) +
        geom_label_repel() +
        ylab('Regression electronegativity (V, relative to N)') +
        xlab('Pauling electronegativity')
    
    ggsave(glue('{category_structure_pair}_electronegativity_comparison_plot.png'),
           electronegativity_comparison_plot,
           height = unit(4.76, 'in'), width = unit(5.67, 'in'))
    
    regression_plot_table <- hardness_terms |>
        rename(symbol = category, this_charge = variable_contribution) |>
        left_join(elements, by = c('combination_id', 'symbol'), relationship = 'many-to-one') |>
        left_join(rename(electronegativity_differences,
                         electronegativity_difference = variable_contribution),
                  by = 'combination_id', relationship = 'many-to-one')
    
    # Make plots showing the data underlying the regression
    # I don't like having to rename here, needs refactoring upstream
    hardness_regression_plot <- regression_plot_table |>
        ggplot(aes(x = this_charge, y = electronegativity_difference, color = other_symbol)) +
        facet_wrap(~ symbol) +
        geom_hline(yintercept = 0, linetype = 'dotted') +
        geom_vline(xintercept = 0, linetype = 'dotted') +
        geom_line() +
        xlab('Charge of acceptor group') +
        ylab('Δpotential (V)')
    ggsave(glue('{category_structure_pair}_hardness_regression_plot.png'), hardness_regression_plot, height = unit(4.76, 'in'), width = unit(11.5, 'in'))
}
