library(broom)
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

# Make a table that maps combination ids to formulas, since these will be used
# to decide which rows to leave out during cross validation
id2formula <- charge_energy |>
    select(combination_id, formula) |>
    distinct()

# Also just a list of unique formulas, which will serve as the list of formula
# to leave out
formula <- id2formula |>
    distinct(formula)

# Symbols of the donor and acceptor in each formula
donor_acceptor_symbols <- charge_energy |>
    distinct(formula, symbol, donor_or_acceptor) |>
    pivot_wider(names_from = donor_or_acceptor, values_from = symbol)

param_regex <- '([^_]*)_([^_]*)_S(\\d+)'
# Function for filtering on one of the reference ranking columns. It's a
# predicate that means "above the minimum rank", but missing values are
# also counted as true, since that means it isn't the right term for this
# kind of referencing
not_reference <- function(x)
{
    minimum_rank <- min(x, na.rm = TRUE)
    is_minimum_rank <- sapply(x == minimum_rank, isTRUE)
    return(!is_minimum_rank)
}

for (category_structure_pair in category_structure_pairs)
{
    electronegativity_terms <- read_csv(glue('{category_structure_pair}_electronegativity_terms.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            category = col_character(),
            variable_contribution = col_double(),
            electronegativity_term_rank = col_integer()
        ))
    hardness_terms <- read_csv(glue('{category_structure_pair}_hardness_terms.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            category = col_character(),
            variable_contribution = col_double(),
            hardness_term_rank = col_integer()
        ))

    electronegativity_differences <- read_csv(glue('{category_structure_pair}_electronegativity_differences.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            variable_contribution = col_double()
        ))

    regression_terms <- bind_rows(
            electronegativity_difference = electronegativity_differences,
            electronegativity = electronegativity_terms,
            hardness = hardness_terms,
            .id = 'variable_type') |>
        # The categories are all {formula}_S{scale_number} or {symbol}_S{scale_number}
        # Therefore variable names will consist of three parts: type, group, scale number
        # (where group is either formula or scale number)
        mutate(variable_name = ifelse(is.na(category), variable_type,
                                      glue('{variable_type}_{category}'))) |>
        # Add formula which will be used to set up the cross-validation
        left_join(id2formula, by = 'combination_id')
    
    # Estimate parameters for the elements in the left out formula, and predict
    # charge of the acceptor
    loo_pre_regression_table <- formula |>
        rename(formula_left_out = formula) |>
        bind_rows(tibble(formula_left_out = 'none')) |>
        cross_join(regression_terms) |>
        filter(formula != formula_left_out) |>
        group_by(formula_left_out) |>
        # Filter out the lowest ranked terms according to each of the reference
        # rank columns and then remove those columns
        filter(not_reference(electronegativity_term_rank) &
               not_reference(hardness_term_rank)) |>
        ungroup() |>
        select(-electronegativity_term_rank,
               -hardness_term_rank) |>
        group_by(combination_id, formula, formula_left_out, variable_name) |>
        summarize(variable_value = sum(variable_contribution), .groups = 'drop') |>
        # Remove the combination id so won't end up being used as a predictor
        # Convert from a tidy table to a design matrix
        pivot_wider(names_from = variable_name, values_from = variable_value,
                    values_fill = 0)

    loo_models <- loo_pre_regression_table |>
        select(-combination_id, -formula) |>
        group_by(formula_left_out) |>
        # Make a list column of regression models
        summarize(model = list(
                lm(electronegativity_difference ~ . + 0,
                   data = pick(everything()))
        ))

    # For each formula, make predictions using the model where that formula was
    # left out
    loo_eq_pred <- loo_pre_regression_table |>
        filter(formula_left_out == 'none') |>
        select(-formula_left_out) |>
        group_by(formula) |>
        nest() |>
        left_join(loo_models, by = c('formula' = 'formula_left_out')) |>
        reframe(
            tibble(combination_id = data[[1]]$combination_id,
                   predicted_electronegativity_difference =
                       predict(model[[1]], newdata = data[[1]]))
        )
    write_csv(loo_eq_pred, glue('{category_structure_pair}_loo_eq_pred.csv.gz'))

    loo_parameters <- loo_models |>
        group_by(formula_left_out) |>
        reframe(tidy(model[[1]]))
    write_csv(loo_parameters, glue('{category_structure_pair}_loo_parameters.csv.gz'))

    loo <- loo_parameters |>
        ungroup() |>
        # Join the information on the symbols of the donor and acceptor
        left_join(donor_acceptor_symbols,
            by = c('formula_left_out' = 'formula')) |>
        # Extract information about what kind of parameter this is and what
        # element it describes
        mutate(
            param_type = str_extract(term, param_regex, group = 1),
            param_symbol = str_extract(term, param_regex, group = 2)
        ) |>
        # Filter for only the parameters of elements that are contained in the
        # formula
        filter(param_symbol == donor | 
               param_symbol == acceptor) |>
        mutate(param = ifelse(param_symbol == donor,
                              glue('{param_type}_donor'), 
                              glue('{param_type}_acceptor'))) |>
        select(formula_left_out, param, estimate) |>
        pivot_wider(names_from = param, values_from = estimate) |>
        # If a symbol is missing that means there wasn't a parameter estimate
        # for it which means I implicitly set it to zero by leaving it out of
        # the regression model. So replace missing values with zero
        mutate(electronegativity_donor = replace_na(electronegativity_donor, 0),
               electronegativity_acceptor = replace_na(electronegativity_acceptor, 0),
               hardness_donor = replace_na(hardness_donor, 0),
               hardness_acceptor = replace_na(hardness_acceptor, 0)) |>
        # Use the standard electronegativity equalization formula to make
        # predictions of charge
        # Electronegativity of the acceptor is higher, but the acceptor charge
        # is negative, so we need a negative sign
        mutate(eeq_acceptor_charge = 
                   -(electronegativity_acceptor - electronegativity_donor) /
                   (hardness_donor + hardness_acceptor)
               ) |>
        rename(formula = formula_left_out)
        
        write_csv(loo, glue('{category_structure_pair}_loo.csv.gz'))
}
