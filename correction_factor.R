# Using the table of charges from integration, calculate the change in charge
# at the last iteration, as well as the change in charge at the first iteration
library(tidyverse)

simulations <- read_csv('simulations.csv', col_types = cols(
    simulation_id = col_character(),
    potential = col_character(),
    structure_id = col_character(),
    cation = col_character(),
    anion = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    log_file_path = col_character(),
    cube_file_path = col_character()
))

charges_from_integration <- read_csv('charges_from_integration.csv', col_types = cols(
    simulation_id = col_character(),
    symbol = col_character(),
    charge = col_double(),
    iteration = col_integer()
))

# Table containing the charge at the first iteration, the last iteration, and
# lagged versions of these columns
lagged_first_last_iter <- simulations |>
    select(simulation_id, structure_id, field_number) |>
    inner_join(charges_from_integration, by = 'simulation_id') |>
    group_by(structure_id, symbol, field_number) |>
    # The iterations are 0 indexed, but iteration 0 really means no update has
    # been done, so I want iteration 1, after the first update
    summarize(first_iteration_charge = charge[iteration == 1],
              last_iteration_charge = charge[iteration == max(iteration)],
              .groups = 'drop') |>
    arrange(structure_id, symbol, field_number) |>
    mutate(lagged_first_iteration_charge = lag(first_iteration_charge),
           lagged_last_iteration_charge = lag(last_iteration_charge)) |>
    # Use these numbers to calculate the correction factor
    mutate(correction_factor_contribution =
           (first_iteration_charge - lagged_last_iteration_charge) /
           (last_iteration_charge - lagged_last_iteration_charge))

write_csv(lagged_first_last_iter, 'correction_factor_contribution.csv.gz')

correction_factor <- lagged_first_last_iter |>
    # Sum across atoms
    group_by(structure_id, field_number) |>
    summarize(correction_factor = sum(correction_factor_contribution),
              .groups = 'drop')
    
write_csv(correction_factor, 'correction_factor.csv.gz')

