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

# Detect combination id's for simulations where integration with AIMAll found
# non-nuclear attractors. Non-nuclear attractors are labeled with 'NNA' in
# AIMAll, but since I converted to title case, they will be labeled 'Nna'
# Table counting non-nuclear attractors
nna_number <- one_tbl |>
    # An arbitrary table of atoms, to look at the names
    filter(section == 'Nuclear Charges and Cartesian Coordinates' &
           property == 'Charge') |>
    # Count the number of non-nuclear attractors
    group_by(combination_id) |>
    mutate(is_nna = str_detect(atom_id, 'Nna')) |>
    summarize(nna_number = sum(is_nna), .groups = 'drop') |>
    select(combination_id, nna_number) |>
    # Join with field information, since really I want to bound the field,
    # not just filter out non-nuclear attractors
    inner_join(simulation_table, by = 'combination_id') |>
    select(combination_id, formula, field_number, field_value, nna_number)
write_csv(nna_number, 'nna_number.csv.gz')

# Table containing combination id's of simulations without non-nuclear
# attractors, for filtering with an inner join. Should also exclude simulations
# without non-nuclear attractors, if the field is more extreme than a
# simulation with a non-nuclear attractor
# This filtering technique doesn't work if there's NNA's without a field
directional_field_nna <- nna_number |>
    # Calculate bounds for the field
    filter(field_value != 0) |>
    mutate(direction = sign(field_value),
           strength = abs(field_value))
nna_at_field_or_weaker <- directional_field_nna |>
    left_join(directional_field_nna, by = c('formula', 'direction'), suffix = c('', '_other'), relationship = 'many-to-many') |>
    # What's the strongest field under the constraint that it has no NNA's,
    # and no weaker field has any NNA?
    filter(strength >= strength_other) |>
    group_by(combination_id) |>
    summarize(nna_at_field_or_weaker = sum(nna_number_other != 0),
              .groups = 'drop')
write_csv(nna_at_field_or_weaker, 'nna_at_field_or_weaker.csv.gz')

no_nna <- nna_at_field_or_weaker |>
    filter(nna_at_field_or_weaker == 0) |>
    select(combination_id)

# Bader charge of each atom
bader_charge_atom <- one_tbl |>
    filter(section == 'Some Atomic Properties' & property == 'q(A)') |>
    group_by(combination_id, atom_id) |>
    transmute(bader_charge = value) |>
    ungroup()

write_csv(bader_charge_atom, 'bader_charge_atom.csv.gz')

# Bader charge of each group
bader_charge_group <- bader_charge_atom |>
    inner_join(no_nna, by = 'combination_id') |>
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
    inner_join(no_nna, by = 'combination_id') |>
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
