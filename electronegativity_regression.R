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
    bader_charge = col_double(),
    cp2k_hirshfeld_charge = col_double(),
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
    scale_number = col_integer(),
    scale = col_double(),
    unscaled_structure_id = col_character(),
    formula = col_character(),
    first_iteration_charge = col_double(),
    last_iteration_charge = col_double(),
    lagged_first_iteration_charge = col_double(),
    lagged_last_iteration_charge = col_double(),
    electronegativity_field_discrete = col_double(),
    unscaled_structure_id = col_character(),
    scale_number = col_integer(),
    scale = col_double()
))

unique_scale_numbers <- charge_energy |>
    distinct(scale_number) |>
    pull(scale_number)

elements <- charge_energy |>
    select(combination_id, symbol, other_symbol) |>
    distinct()

scale_numbers <- charge_energy |>
    select(combination_id, scale_number) |>
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
               charge, electronegativity_field_discrete,
               scale_number)
    
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

    # As a reference for the interaction terms, all formulas with the "minimum"
    # (by atomic number) cation or anion
    reference_formulas <- cdft_charges |>
        distinct(formula, symbol, other_symbol, donor_or_acceptor) |>
        left_join(atomic_numbers, by = 'symbol') |>
        group_by(donor_or_acceptor) |>
        filter(atomic_number == min(atomic_number)) |>
        distinct(formula) |>
        # Also add a scale number that I'm going to use as a reference. Only
        # the interaction terms with that scale number need to be removed
        mutate(scale_number = 1)
    write_csv(reference_formulas, glue('{category_structure_pair}_reference_formulas.csv'))

    electronegativity_terms <- cdft_charges |>
        # Electronegativity terms correspond to an element
        # Also make them dependent on scale
        mutate(category = glue('{symbol}_S{scale_number}'),
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1)) |>
        select(combination_id, category, variable_contribution) |>
        # Remove electronegativity of an arbitrary element to use it as a
        # reference
        filter(category != reference_element)
    
    hardness_terms <- cdft_charges |>
        # Same as electronegativity: each term corresponds to an element, also
        # dependent on scale
        mutate(category = glue('{symbol}_S{scale_number}'),
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
        # I think this only needs to be done for one scale though, so the table
        # actually also contains scale numbers, even though I haven't updated
        # the name from "reference formulas"
        anti_join(reference_formulas, by = c('formula', 'scale_number')) |>
        mutate(category = glue('{formula}_S{scale_number}'),
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
        # The categories are all {formula}_S{scale_number} or {symbol}_S{scale_number}
        # Therefore variable names will consist of three parts: type, group, scale number
        # (where group is either formula or scale number)
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
    
    electronegativity_param_regex <- '([^_]*)_([^_]*)_S(\\d+)'
    regression_electronegativity <- estimates |>
        mutate(
            term_type = str_match(term, electronegativity_param_regex)[, 2],
            symbol = str_match(term, electronegativity_param_regex)[, 3],
            scale_number = as.integer(str_match(term, electronegativity_param_regex)[, 4])
        ) |>
        filter(term_type == 'electronegativity') |>
        select(symbol, scale_number, estimate) |>
        rename(regression_electronegativity = estimate) |>
        bind_rows(tibble(symbol = reference_element,
                         regression_electronegativity = 0,
                         scale_number = unique_scale_numbers))
    
    electronegativity_comparison_tbl <- left_join(regression_electronegativity,
                                                  pauling_electronegativity,
                                                  by = 'symbol')
    
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
    # I don't like having to rename here, needs refactoring upstream
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
        xlab('Charge of acceptor group') +
        ylab('Δpotential (V)')
    ggsave(glue('{category_structure_pair}_hardness_regression_plot.png'), hardness_regression_plot, height = unit(4.76, 'in'), width = unit(11.5, 'in'))
}
