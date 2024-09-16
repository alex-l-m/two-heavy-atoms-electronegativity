library(tidyverse)
library(glue)

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
    charge = col_double()
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
    first_energy <- as.numeric(first_energy_str)

    return(first_energy)
}

# Get the charges from the simulations with a field applied
charges <- simulations |>
    # Filter should be redundant since charges should only be calculated with
    # the field
    filter(potential == 'field') |>
    inner_join(charges_from_integration, by = 'simulation_id') |>
    # Assign a combination id
    mutate(combination_id = glue('{structure_id}_F{field_number}')) |>
    # Decide donor or acceptor status based on whether the element symbol
    # matches the cation or the anion
    mutate(donor_or_acceptor = if_else(symbol == cation, 'donor',
                                    if_else(symbol == anion, 'acceptor',
                                            NA))) |>
    # Retrieve the other symbol, either cation or anion
    mutate(other_symbol = if_else(donor_or_acceptor == 'donor', anion,
                                  if_else(donor_or_acceptor == 'acceptor', cation,
                                          NA))) |>
    select(combination_id, symbol, other_symbol, donor_or_acceptor, charge)

# Get the initial energies from each log file
initial_energies <- simulations |>
    mutate(total_energy = map_dbl(log_file_path, extract_first_energy))
# The energy I'm interested in is the energy with no field applied
# Get the no field energy and assign a combination idea
energies <- initial_energies |>
    filter(potential == 'nuclei') |>
    mutate(combination_id = glue('{structure_id}_F{field_number}')) |>
    # This is also where I will keep the simulation annotations that should
    # depend only on the combination
    select(combination_id, structure_id, field_number, field_value, total_energy)

# Make the table of charges and energies
charge_energy <- charges |>
    inner_join(energies, by = 'combination_id') |>
    # Column that indicates whether a field was applied
    # Called "cdft" for backwards compatibility reasons for now
    mutate(cdft = field_number > 0) |>
    # Assign a "group_id" by combining information fom other columns
    mutate(group_id = glue('{structure_id}:{symbol}:{donor_or_acceptor}'))

# Join with metadata for structures, including crystal structure, for easy plotting later
structure_metadata <- read_csv('selected_structure_files.csv', col_types = cols(
    category = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character(),
    crystal_structure = col_character(),
    structure_file_path = col_character(),
    structure_id = col_character()
)) |>
    select(-structure_file_path) |>
    mutate(formula = glue('{symbol_cation}{symbol_anion}'))
charge_energy_annotated <- charge_energy |>
    left_join(structure_metadata, by = 'structure_id')

write_csv(charge_energy_annotated, 'charge_energy.csv.gz')

