library(tidyverse)
library(glue)

TO_EV <- 27.211386246

simulations <- read_csv('simulations.csv', col_types = cols(
    simulation_id = col_character(),
    potential = col_character(),
    structure_id = col_character(),
    cation = col_character(),
    anion = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    log_file_path = col_character(),
    cube_file_path = col_character(),
    pot_file_path = col_character(),
    hartree_pot_path = col_character()
))

# Function for extracting the first-iteration energy from a log file
# I don't like that I have to change this whenever I change the update method
energy_line_regex <- ' *1 +P_Mix/Diag. +[0-9.E+-]+ +[0-9.E+-]+ +[0-9.E+-]+ +([0-9.E+-]+) +[0-9.E+-]+'
extract_first_energy <- function(inpath)
{
    # Check if the file exists
    if (!file.exists(inpath))
    {
        return(NA)
    }
    intext <- read_file(inpath)
    
    first_energy_str <- str_match(intext, energy_line_regex)[1,2]
    first_energy_au <- as.numeric(first_energy_str)
    first_energy <- first_energy_au * TO_EV

    return(first_energy)
}

# Get the initial energies from each log file
initial_energies <- simulations |>
    mutate(total_energy = map_dbl(log_file_path, extract_first_energy))
# The energy I'm interested in is the energy with no field applied
# Get the no field energy and assign a combination id
energies <- initial_energies |>
    filter(potential == 'nuclei') |>
    mutate(combination_id = glue('{structure_id}_F{field_number}')) |>
    # This is also where I will keep the simulation annotations that should
    # depend only on the combination
    select(combination_id, structure_id, field_number, field_value, total_energy)

# Write the table of energies
write_csv(energies, 'energies.csv.gz')
