library(broom)
library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
library(glue)
library(ggrepel)

# Atomic numbers, for ordering
atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_double()
))

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
    formula = col_character(),
    electronegativity_field_discrete = col_double()
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
               charge, electronegativity_field_discrete)
    
    # Order elements by atomic number, to select the least for use as a
    # reference for electronegativity and hardness
    ordered_selected_elements <- cdft_charges |>
        select(symbol) |>
        distinct() |>
        left_join(atomic_numbers, by = 'symbol') |>
        arrange(atomic_number) |>
        pull(symbol)
    reference_element <- ordered_selected_elements[1]
    write_csv(tibble(symbol = reference_element),
              glue('{category_structure_pair}_reference_element.csv'))

    # For each element, the formula with the minimum other element, to use as a
    # reference for the interaction term
    minimum_formulas <- cdft_charges |>
        select(formula, symbol, other_symbol) |>
        distinct() |>
        left_join(atomic_numbers, by = c('symbol' = 'symbol')) |>
        left_join(atomic_numbers, by = c('other_symbol' = 'symbol'),
                                  suffix = c('', '_other')) |>
        group_by(symbol) |>
        summarize(minimum_formula = formula[which.min(atomic_number_other)],
                  .groups = 'drop')
    reference_formulas <- minimum_formulas |>
        rename(formula = minimum_formula) |>
        distinct(formula)
    write_csv(reference_formulas, glue('{category_structure_pair}_reference_formulas.csv'))

    electronegativity_terms <- cdft_charges |>
        mutate(category = symbol,
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1)) |>
        select(combination_id, category, variable_contribution) |>
        # Remove electronegativity of an arbitrary element to use it as a
        # reference
        filter(category != reference_element)
    
    hardness_terms <- cdft_charges |>
        mutate(category = symbol,
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution) |>
        # Remove hardness of an arbitrary element to use it as a
        # reference
        # This is harder to explain, but since all that matters is total
        # hardness, and there's always a cation and an anion, you can "move"
        # hardness from all cations to all anions without changing any of the
        # totals. So we decide in advance how much to "move": whatever makes
        # the hardness of the reference element zero
        filter(category != reference_element)
    
    interaction_terms <- cdft_charges |>
        # There's linear constraint on the interaction terms as well. Every
        # atom can be thought of as having an "adjusted hardness", which is the
        # hardness minus the interaction term. So the "true" hardness is never
        # realized in practice. So you could just consider the smallest
        # interaction term to be zero, and the hardness to be the highest
        # achievable hardness
        anti_join(reference_formulas, by = 'formula') |>
        mutate(category = formula, 
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution)
    
    electronegativity_differences <- cdft_charges |>
        # Sum the electronegativities for each atom to get the derivative of
        # energy as one atom gains charge and the other loses it
        group_by(combination_id) |>
        summarize(variable_contribution = sum(electronegativity_field_discrete),
                  .groups = 'drop') |>
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
        group_by(combination_id, variable_name) |>
        summarize(variable_value = sum(variable_contribution), .groups = 'drop') |>
        pivot_wider(names_from = variable_name, values_from = variable_value,
                    values_fill = 0)
    
    regression_table <- pre_regression_table |>
        select(-combination_id) |>
        as.data.frame()
    rownames(regression_table) <- pre_regression_table$combination_id
    
    linear_model <- lm(electronegativity_difference ~ . + 0,
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
        bind_rows(tibble(symbol = reference_element,
                         regression_electronegativity = 0))
    
    electronegativity_comparison_tbl <- left_join(regression_electronegativity,
                                                  pauling_electronegativity,
                                                  by = 'symbol')
    
    electronegativity_comparison_plot <- electronegativity_comparison_tbl |>
        ggplot(aes(x = pauling_electronegativity, y = regression_electronegativity,
                   label = symbol)) +
        geom_point() +
        geom_smooth(method = lmrob, se = FALSE) +
        geom_label_repel() +
        ylab('Regression electronegativity (V)') +
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
        # Turn the symbols into factors based on the ordering implied by the
        # atomic numbers
        mutate(symbol = factor(symbol, levels = ordered_selected_elements),
               other_symbol = factor(other_symbol, levels = ordered_selected_elements)) |>
        ggplot(aes(x = this_charge, y = electronegativity_difference, color = other_symbol)) +
        facet_wrap(~ symbol) +
        geom_hline(yintercept = 0, linetype = 'dotted') +
        geom_vline(xintercept = 0, linetype = 'dotted') +
        geom_line() +
        xlab('Charge of acceptor group') +
        ylab('Δpotential (V)')
    ggsave(glue('{category_structure_pair}_hardness_regression_plot.png'), hardness_regression_plot, height = unit(4.76, 'in'), width = unit(11.5, 'in'))
}
