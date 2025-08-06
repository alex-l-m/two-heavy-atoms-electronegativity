library(tidyverse)
library(glue)

# Atomic numbers, for ordering
atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_integer()
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
    ranked_elements <- cdft_charges |>
        select(symbol) |>
        distinct() |>
        left_join(atomic_numbers, by = 'symbol') |>
        # Ranks
        # Choosing a deterministic ties method with integer output, but there
        # shouldn't be any ties, so it doesn't matter
        mutate(element_rank = rank(atomic_number, ties.method = 'first'))

    # As a reference for the interaction terms, all formulas with the "minimum"
    # (by atomic number) cation or anion
    ranked_formulas <- cdft_charges |>
        distinct(formula, symbol, donor_or_acceptor) |>
        left_join(atomic_numbers, by = 'symbol') |>
        group_by(donor_or_acceptor) |>
        # Ranks
        # This time, multiple formulas can be the minimum, I need a method
        # for handling ties that doesn't give unique ranks
        # Therefore, using the 'min' method, which gives output like:
        # > rank(c(15, 15, 20), ties.method = 'min')
        # [1] 1 1 3
        mutate(formula_rank = rank(atomic_number, ties.method = 'min')) |>
        select(formula, donor_or_acceptor, formula_rank) |>
        # For each formula there should be two distinct rank columns, one for
        # the donor and one for the acceptor
        pivot_wider(
            names_from = donor_or_acceptor,
            values_from = formula_rank,
            names_prefix = 'formula_rank_'
        )

    electronegativity_terms <- cdft_charges |>
        # Add the rank column for referencing
        left_join(ranked_elements, by = 'symbol') |>
        rename(electronegativity_term_rank = element_rank) |>
        # Electronegativity terms correspond to an element
        # Also make them dependent on scale
        mutate(category = glue('{symbol}_S{scale_number}'),
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1)) |>
        select(combination_id, category, variable_contribution, 
               electronegativity_term_rank)

    write_csv(electronegativity_terms,
              glue('{category_structure_pair}_electronegativity_terms.csv.gz'))
    
    hardness_terms <- cdft_charges |>
        # Hardness for one arbitrary element has to be removed to use it as a
        # reference
        # This is harder to explain, but since all that matters is total
        # hardness, and there's always a cation and an anion, you can "move"
        # hardness from all cations to all anions without changing any of the
        # totals. So we decide in advance how much to "move": whatever makes
        # the hardness of the reference element zero
        # Add the rank column for later removal
        left_join(ranked_elements, by = 'symbol') |>
        rename(hardness_term_rank = element_rank) |>
        # Same as electronegativity: each term corresponds to an element, also
        # dependent on scale
        mutate(category = glue('{symbol}_S{scale_number}'),
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution,
               hardness_term_rank)

    write_csv(hardness_terms, glue('{category_structure_pair}_hardness_terms.csv.gz'))
    
    interaction_terms <- cdft_charges |>
        # There's linear constraint on the interaction terms as well. Every
        # atom can be thought of as having an "adjusted hardness", which is the
        # hardness minus the interaction term. So the "true" hardness is never
        # realized in practice. So you could just consider the smallest
        # interaction term to be zero, and the hardness to be the highest
        # achievable hardness
        left_join(ranked_formulas, by = 'formula') |>
        rename(interaction_term_donor_rank = formula_rank_donor,
               interaction_term_acceptor_rank = formula_rank_acceptor) |>
        mutate(category = glue('{formula}_S{scale_number}'),
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution,
               interaction_term_donor_rank, interaction_term_acceptor_rank)

    write_csv(interaction_terms, glue('{category_structure_pair}_interaction_terms.csv.gz'))
    
    electronegativity_differences <- cdft_charges |>
        # Sum the electronegativities for each atom to get the derivative of
        # energy as one atom gains charge and the other loses it
        group_by(combination_id) |>
        summarize(variable_contribution = sum(electronegativity_field_discrete),
                  .groups = 'drop') |>
        select(combination_id, variable_contribution)
    write_csv(electronegativity_differences,
              glue('{category_structure_pair}_electronegativity_differences.csv.gz'))
}
