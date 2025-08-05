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
for (category_structure_pair in category_structure_pairs)
{
    electronegativity_terms <- read_csv(glue('{category_structure_pair}_electronegativity_terms.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            category = col_character(),
            variable_contribution = col_double()
        ))
    hardness_terms <- read_csv(glue('{category_structure_pair}_hardness_terms.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            category = col_character(),
            variable_contribution = col_double()
        ))
    interaction_terms <- read_csv(glue('{category_structure_pair}_interaction_terms.csv.gz'),
        col_types = cols(
            combination_id = col_character(),
            category = col_character(),
            variable_contribution = col_double()
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
    write_csv(estimates, glue('{category_structure_pair}_electronegativity_regression_estimates.csv.gz'))
}
