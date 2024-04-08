# Make a table showing the status of the simulations and integrations

library(tidyverse)

# Combination id of simulations where Q-Chem completed, as indicated by
# existence of the .wfn file produced
qchem_output_paths <- Sys.glob('aimall/*.wfn')
qchem_output_filenames <- basename(qchem_output_paths)
qchem_complete_combinations <-
    str_extract(qchem_output_filenames, '^(.*)\\.wfn$', 1)

# Combination id of simulations where integration of atomic regions with AIMAll
# has completed, as indicated by existence of an output file
integration_results_paths <- Sys.glob('aimall/*.sum')
integration_results_filenames <- basename(integration_results_paths)
integration_complete_combinations <- 
    str_extract(integration_results_filenames, '^(.*)\\.sum$', 1)

# Combination id of simulations where integration didn't work
bad_integration_list_filename <- 'bad_integration_combination_id.txt'
bad_integration_combinations <- read_lines(bad_integration_list_filename)

# Load table of simulations
simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer()
))

# Table of simulation status
simulation_status <- simulation_table |>
    select(combination_id) |>
    mutate(
        qchem_complete = combination_id %in% qchem_complete_combinations,
        integration_complete = combination_id %in% integration_complete_combinations,
        accurate_integration = ! ifelse(integration_complete,
                                        combination_id %in% bad_integration_combinations,
                                        NA),
    )

write_csv(simulation_status, 'simulation_status.csv.gz')
