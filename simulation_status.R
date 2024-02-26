# Make a table showing the status of the simulations and integrations

library(tidyverse)

# Combination id of simulations where GAMESS completed, as indicated by
# existence of an output file
gamess_output_paths <- Sys.glob('gamess_output/*.dat')
gamess_output_filenames <- basename(gamess_output_paths)
gamess_complete_combinations <-
    str_extract(gamess_output_filenames, '^(.*)\\.dat$', 1)

# Load premade list of combination ids of unconverged simulations
unconverged_list_filename <- 'unconverged_combination_id.txt'
unconverged_combinations <- read_lines(unconverged_list_filename)

# Combination id of simulations where integration of atomic regions with AIMAll
# has completed, as indicated by existence of an output file
integration_results_paths <- Sys.glob('aimall/*.sum')
integration_results_filenames <- basename(integration_results_paths)
integration_complete_combinations <- 
    str_extract(integration_results_filenames, '^(.*)\\.sum$', 1)

# Load table of simulations
simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer(),
    gamess_input_file = col_character()
))

# Table of simulation status
simulation_status <- simulation_table |>
    select(combination_id) |>
    mutate(
        gamess_complete = combination_id %in% gamess_complete_combinations,
        converged = ! ifelse(gamess_complete, combination_id %in% unconverged_combinations, NA),
        integration_complete = combination_id %in% integration_complete_combinations
    )

write_csv(simulation_status, 'simulation_status.csv.gz')
