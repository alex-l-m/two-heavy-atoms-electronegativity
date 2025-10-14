library(tidyverse)
library(glue)

cell_sizes <- read_csv('cell_sizes.csv', col_types = cols(
    structure_id = col_character(),
    cell_size = col_double()
))

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
    # Charges which will be used in this iteration of the loop
    these_charges <- charge_energy |>
        filter(glue('{category}:{crystal_structure}') == category_structure_pair) |>
        select(combination_id, formula, symbol, other_symbol, donor_or_acceptor,
               charge, electronegativity_field_analytic,
               scale_number,
               # Also needed structure_id for joining cell sizes
               # In the future might just want cell sizes in charge_energy
               structure_id)
    
    # Order elements by atomic number, to select the least for use as a
    # reference for electronegativity and hardness
    ranked_elements <- these_charges |>
        select(symbol) |>
        distinct() |>
        left_join(atomic_numbers, by = 'symbol') |>
        # Ranks
        # Choosing a deterministic ties method with integer output, but there
        # shouldn't be any ties, so it doesn't matter
        mutate(element_rank = rank(atomic_number, ties.method = 'first'))
    write_csv(ranked_elements,
              glue('{category_structure_pair}_ranked_elements.csv.gz'))

    # As a reference for the interaction terms, all formulas with the "minimum"
    # (by atomic number) cation or anion
    ranked_formulas <- these_charges |>
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
    write_csv(ranked_formulas,
              glue('{category_structure_pair}_ranked_formulas.csv.gz'))

    electronegativity_terms <- these_charges |>
        # Add the rank column for referencing
        left_join(ranked_elements, by = 'symbol') |>
        rename(electronegativity_term_rank = element_rank) |>
        # Electronegativity terms correspond to an element
        # Also making them dependent on scale, since I have enough data to fit
        # scale-dependence, but may change this decision if I ever actually
        # include multiple scales in the analysis
        # Also need to reverse sign for a donor; this is because we're
        # constructing the regression equation for the electronegativity
        # difference, that is, the electronegativity of the acceptor minus the
        # electronegativity of the donor, so donor terms are negative
        mutate(category = glue('{symbol}_S{scale_number}'),
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1)) |>
        select(combination_id, category, variable_contribution, 
               electronegativity_term_rank)

    write_csv(electronegativity_terms,
              glue('{category_structure_pair}_electronegativity_terms.csv.gz'))
    
    hardness_terms <- these_charges |>
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
               # Note that the contributions from both atoms end up the same
               # sign; one is added and the other is subtracted, but one has
               # positive charge and the other negative
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * charge) |>
        select(combination_id, category, variable_contribution,
               hardness_term_rank)

    write_csv(hardness_terms, glue('{category_structure_pair}_hardness_terms.csv.gz'))
    
    # The interaction terms the way I originally did it, where they are
    # arbitrary parameters that depend on the material. The resulting model
    # can't be used for prediction. Originally I did it this way just because I
    # didn't want to take a chance with any specific model of the interaction,
    # except the implicit assumption that it doesn't depend on the charges
    # (which would be violated by an overlap correction that depends on the
    # sizes of the Hirshfeld atoms). I still think this analysis could be
    # useful for something, maybe just for comparing these terms to models
    # based on stronger assumptions, so I'm keeping this code.
    interaction_terms <- these_charges |>
        # There's linear constraint on the interaction terms as well. Every
        # atom can be thought of as having an "adjusted hardness", which is the
        # hardness minus the interaction term. So the "true" hardness is never
        # realized in practice. So you could just consider the smallest
        # interaction term to be zero, and the hardness to be the highest
        # achievable hardness
        # Another issue to consider is sign. For an interaction term, the
        # relevant charge is the charge of the other atom. Using the assumption
        # of zero total charge, that's a sign flip. This will have to change if
        # I ever have non-constant total charge.
        # Contrasting to the hardness, this shows that the interaction
        # "softens" the atoms.
        left_join(ranked_formulas, by = 'formula') |>
        rename(interaction_term_donor_rank = formula_rank_donor,
               interaction_term_acceptor_rank = formula_rank_acceptor) |>
        mutate(category = glue('{formula}_S{scale_number}'),
               variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * (-charge)) |>
        select(combination_id, category, variable_contribution,
               interaction_term_donor_rank, interaction_term_acceptor_rank)

    write_csv(interaction_terms, glue('{category_structure_pair}_interaction_terms.csv.gz'))

    # Naive Coulomb interaction term, but with a constant learned by
    # regression. In the model where the material is a lattice of point
    # charges, the constant reflects the universal Coulomb constant and the
    # crystal-structure specific Madelung constant. And, depending on what
    # measure of cell size I use, some adjustment for that. It should really
    # also depend on the dielectric constant and on the overlap. I'm leaving
    # these out of this model, but as a result the constant that is estimated
    # from the data cannot be expected to match the theoretical value for a
    # lattice of point charges.
    point_charge_coulomb_terms <- these_charges |>
        left_join(cell_sizes, by = 'structure_id') |>
        # The category here is tricky. It will have to get more specific if I
        # ever include multiple space groups in the same regression, because
        # each will have its own Madelung constant. If I include multiple
        # scales, the same one should be used across scales, because I take
        # into account changes in the lattice constant with the cell size. So
        # for now, I won't have any category.
        mutate(variable_contribution = ifelse(donor_or_acceptor == 'acceptor',
                                              1, -1) * (-charge) / cell_size) |>
        select(combination_id, variable_contribution)
    write_csv(point_charge_coulomb_terms,
              glue('{category_structure_pair}_point_charge_coulomb_terms.csv.gz'))

    # Sum the electronegativities for each atom to get the derivative of energy
    # as one atom gains charge and the other loses it
    # Theoretically I should be able to use the discrete or analytic
    # electronegativities here. I'm not even sure which is more accurate since
    # the "analytic" one relies on a numerical integration. Currently the
    # discrete one is undefined for the first frame (though that depends on how
    # I take the numerical derivative, so it may change), so using analytic
    electronegativity_differences <- these_charges |>
        group_by(combination_id) |>
        summarize(variable_contribution = sum(electronegativity_field_analytic),
                  .groups = 'drop') |>
        select(combination_id, variable_contribution)
    write_csv(electronegativity_differences,
              glue('{category_structure_pair}_electronegativity_differences.csv.gz'))
}
