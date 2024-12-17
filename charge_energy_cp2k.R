library(tidyverse)
library(glue)

TO_EV <- 27.211386246

simulations <- read_csv('simulations.csv', col_types = cols(
    simulation_id = col_character(),
    potential = col_character(),
    structure_id = col_character(),
    cation = col_character(),
    anion = col_character(),
    field_number = col_integer(),
    field_value = col_double(),
    log_file_path = col_character(),
    cube_file_path = col_character(),
    pot_file_path = col_character(),
    hartree_pot_path = col_character()
))

# Read table of energies parsed from the log files
energies <- read_csv('energies.csv.gz', col_types = cols(
    combination_id = col_character(),
    structure_id = col_character(),
    field_number = col_integer(),
    field_value = col_double(),
    total_energy = col_double()
))

# Reading the Hirshfeld-I charges that I calculated by integrating against the
# cube file
# Why doesn't this have donor_or_acceptor?
# Won't there be duplicates in 4-4 materials?
charges_from_integration_all_iterations <- read_csv('charges_from_integration.csv', col_types = cols(
    simulation_id = col_character(),
    symbol = col_character(),
    charge = col_double(),
    iteration = col_integer()
))
charges_from_integration <- charges_from_integration_all_iterations |>
    # Select only the last iteration, which is the Hirshfeld-I charge
    group_by(simulation_id, symbol) |>
    filter(iteration == max(iteration)) |>
    select(-iteration)

# Reading the Bader charges
bader_charges_raw <- read_csv('bader_charges_raw.csv.gz', col_types = cols(
    Id = col_double(),
    cp = col_double(),
    ncp = col_double(),
    Name = col_character(),
    Z = col_double(),
    mult = col_character(),
    Volume = col_double(),
    Pop = col_double(),
    Lap = col_double(),
    simulation_id = col_character()
))

# Number of valence electrons. Needed to convert population into charge. though
# that will be later, since the element is not yet available to join on
n_valence_electrons <- read_csv('n_valence_electrons.csv', col_types = cols(
    symbol = col_character(),
    valence_electrons = col_integer()
))
bader_charges <- bader_charges_raw |>
    group_by(simulation_id) |>
    transmute(donor_or_acceptor = ifelse(Id == 1, 'donor', ifelse(Id == 2, 'acceptor', NA)),
              bader_population = Pop) |>
    ungroup()

# Reading the Hirshfeld charges calculated by CP2K
cp2k_hirshfeld_charges <- read_csv('cp2k_hirshfeld_charges.csv.gz', col_types = cols(
    simulation_id = col_character(),
    symbol = col_character(),
    donor_or_acceptor = col_character(),
    charge = col_double()
)) |>
    rename(cp2k_hirshfeld_charge = charge)

# Get the charges from the simulations with a field applied
charges <- simulations |>
    # Filter should be redundant since charges should only be calculated with
    # the field
    filter(potential == 'field') |>
    inner_join(charges_from_integration, by = 'simulation_id',
    # Every simulation should correspond to two charges, one for each atom, so
    # this relationship should be one-to-many
               relationship='one-to-many') |>
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
    # Join the Bader charges
    left_join(bader_charges, by = c('simulation_id', 'donor_or_acceptor'),
              # Join columns identify an atom, so it should be one to one
              relationship = 'one-to-one') |>
    # This actually gives Bader populations. Join number of valence electrons
    # and calculate charges
    left_join(n_valence_electrons, by = 'symbol',
              relationship = 'many-to-one') |>
    mutate(bader_charge = valence_electrons - bader_population) |>
    # Join CP2K's Hirshfeld charges
    # This is on the level of atoms, so one-to-one
    left_join(cp2k_hirshfeld_charges,
              by = c('simulation_id', 'symbol', 'donor_or_acceptor'),
              # This is on the level of atoms, so one-to-one
              relationship = 'one-to-one') |>
    select(combination_id, symbol, other_symbol, donor_or_acceptor, charge, bader_charge, cp2k_hirshfeld_charge)

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
    structure_id = col_character(),
    scale_number = col_integer(),
    scale = col_double(),
    unscaled_structure_id = col_character()
)) |>
    select(-structure_file_path) |>
    mutate(formula = glue('{symbol_cation}{symbol_anion}'))
charge_energy_annotated <- charge_energy |>
    left_join(structure_metadata,
              # Different lattice constants are distinct stuctures with their
              # own structure id that includes the scale, so structure_id is
              # sufficient to join on
              by = 'structure_id',
              # Since a charge corresponds to a simulation, and there are
              # multiple simulations of the same structure, this relationship
              # is expected to be many-to-one
              relationship = 'many-to-one')

# Make a table of the ground state charges for each structure
ground_state_charges <- charge_energy_annotated |>
    filter(donor_or_acceptor == 'acceptor' & field_number == 0) |>
    select(structure_id, donor_or_acceptor, charge) |>
    rename(ground_state_charge = charge)

# Hack: Remove all simulations that go past an integer value of the
# corresponding structure
combinations_to_keep <- charge_energy_annotated |>
    inner_join(ground_state_charges,
               by = c('structure_id', 'donor_or_acceptor'),
    # All simulations of the same structure are being compared to the same
    # ground state, so this relationship should be many-to-one
               relationship = 'many-to-one') |>
    filter(charge < ceiling(ground_state_charge)) |>
    select(combination_id)
charge_energy_annotated <- charge_energy_annotated |>
    inner_join(combinations_to_keep, by = 'combination_id')

# Calculation of correction factor relating field values to electronegativity
# Right now the calcultaion is based on discrete differences, but I'll probably
# add one based on integrating the density later
# Table containing the charge at the first iteration, the last iteration, and
# lagged versions of these columns
electronegativity_field_discrete <- simulations |>
    select(simulation_id, structure_id, field_number) |>
    inner_join(charges_from_integration_all_iterations,
               by = 'simulation_id') |>
    group_by(structure_id, symbol, field_number) |>
    # The iterations are 0 indexed, but iteration 0 really means no update has
    # been done, so I want iteration 1, after the first update
    # Next step will be a discrete difference across field numbers, so only
    # drop that group
    summarize(first_iteration_charge = charge[iteration == 1],
              last_iteration_charge = charge[iteration == max(iteration)],
              .groups = 'drop_last') |>
    arrange(structure_id, symbol, field_number) |>
    mutate(lagged_first_iteration_charge = lag(first_iteration_charge),
           lagged_last_iteration_charge = lag(last_iteration_charge)) |>
    # Use these numbers to calculate the correction factor
    mutate(correction_factor_contribution =
           (first_iteration_charge - lagged_last_iteration_charge) /
           (last_iteration_charge - lagged_last_iteration_charge)) |>
    ungroup()
    
charge_energy_electronegativity <- charge_energy_annotated |>
    left_join(electronegativity_field_discrete,
              # This join looks wrong. Symbol doesn't identify an atom in the
              # 4-4's. Need to change this if I include the 4-4's again?
              by = c('structure_id', 'symbol', 'field_number'),
              relationship = 'one-to-one') |>
    mutate(electronegativity_field_discrete =
           field_value * TO_EV * correction_factor_contribution) |>
    # Actually I don't need the "correction factor contribution" anymore since
    # I'm not summing across atoms
    select(-correction_factor_contribution)

write_csv(charge_energy_electronegativity, 'charge_energy.csv.gz')
