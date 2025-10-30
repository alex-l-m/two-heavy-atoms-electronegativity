library(tidyverse)

# Make the supplementary table showing the maximum deviation between the two
# methods of computing the electronegativity difference
electronegativity_comparison_tbl <- read_csv('3-5:zincblende:0_lam_comparison_table.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    group_id = col_character(),
    charge = col_double(),
    crystal_structure = col_character(),
    category = col_character(),
    structure_id = col_character(),
    symbol_anion = col_character(),
    symbol_cation = col_character(),
    formula = col_character(),
    scale_number = col_integer(),
    cdft = col_logical(),
    field_number = col_double(),
    field_value = col_double(),
    total_energy = col_double(),
    electronegativity_difference_from_field = col_double(),
    computation = col_character(),
    electronegativity = col_double()
))

electronegativity_max_deviations <- electronegativity_comparison_tbl |>
    select(combination_id, structure_id, computation, electronegativity) |>
    pivot_wider(id_cols = c(structure_id, combination_id), names_from = computation, values_from = electronegativity) |>
    mutate(absolute_error = abs(`correction * Lam` - `dE/dq`)) |>
    group_by(structure_id) |>
    summarize(maximum_absolute_error = max(absolute_error, na.rm = TRUE),
              .groups = 'drop')

write_csv(electronegativity_max_deviations, 'manuscript/supplementary_table_electronegativity_max_deviations.csv')
