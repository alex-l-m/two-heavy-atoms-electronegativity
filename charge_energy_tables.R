library(readr)
library(dplyr)
library(tidyr)
library(stringr)

one_tbl <- tibble(path = Sys.glob('aimall_tbl/*_oneatom.csv')) |>
    mutate(combination_id = str_match(path, 'aimall_tbl/(.*)_oneatom.csv')[,2]) |>
    group_by(combination_id) |>
    reframe(read_csv(path, col_names = c('section', 'atom_id', 'property', 'value'), col_types = cols(.default = col_character(), value = col_double()))) |>
    # AIMAll seems to maintain atom order but rename them, one-indexed in all
    # caps. To match my own naming scheme I need to capitalize only the first
    # letter
    mutate(atom_id = str_to_title(atom_id))

write_csv(one_tbl, 'combined_one_tbl.csv.gz')

# Coordinates of each atom
coordinates <- read_csv('coordinates.csv.gz', col_types = cols(
    formula = col_character(),
    atom_id = col_character(),
    symbol = col_character(),
    atomic_number = col_integer(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    donor_or_acceptor = col_character()
))

simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer(),
    gamess_input_file = col_character()
))

# Bader charge of each atom
bader_charge_atom <- one_tbl |>
    filter(section == 'Some Atomic Properties' & property == 'q(A)') |>
    group_by(combination_id, atom_id) |>
    transmute(bader_charge = value) |>
    ungroup()

write_csv(bader_charge_atom, 'bader_charge_atom.csv.gz')

# Bader charge of each group
bader_charge_group <- bader_charge_atom |>
    left_join(simulation_table, by = 'combination_id') |>
    left_join(coordinates, by = c('formula', 'atom_id'),
              relationship = 'many-to-one') |>
    group_by(combination_id, formula, donor_or_acceptor) |>
    summarize(
        heavy_atom_symbol = symbol[symbol != 'H'][1],
        total_bader_charge = sum(bader_charge),
        .groups = 'drop'
    ) |>
    rename(symbol = heavy_atom_symbol)

write_csv(bader_charge_group, 'bader_charge_group.csv.gz')

# The total energies of each molecule. These don't seem to be available
# directly, so I obtain them by summing the atomic energies
TO_EV <- 27.211386246
combined_atom_energies <- tibble(path = Sys.glob('aimall_tbl/*_oneatom.csv')) |>
    mutate(combination_id = str_match(path, 'aimall_tbl/(.*)_oneatom.csv')[,2]) |>
    group_by(combination_id) |>
    # Read each table, treating the value column as numeric
    reframe(read_csv(path,
        col_names = c('section', 'atom', 'property', 'value'),
        col_types = cols(.default = col_character(), value = col_double()))) |>
    # Filter for the IQA energy property. I just sum this to get total energy
    filter(section == 'IQA Atomic Energy Components', property == 'E_IQA(A)') |>
    # Remove the property column and rename the energy column
    select(-property) |>
    rename(energy = value)
total_atom_energies <- combined_atom_energies |>
    group_by(combination_id) |>
    summarize(energy = TO_EV * sum(energy), .groups = 'drop')

write_csv(total_atom_energies, 'total_atom_energies.csv.gz')
