library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)

# Make a table for filtering based on whether integration was accurate. If I
# don't filter these out, I get simulations of neutral molecules where the
# Bader charges don't add up to zero.
simulation_status <- read_csv('simulation_status.csv.gz', col_types = cols(
    combination_id = col_character(),
    integration_complete = col_logical(),
    accurate_integration = col_logical()
))
accurate_integration <- simulation_status |>
    filter(accurate_integration) |>
    select(combination_id)

one_tbl <- tibble(path = Sys.glob('aimall_tbl/*_oneatom.csv')) |>
    mutate(combination_id = str_match(path, 'aimall_tbl/(.*)_oneatom.csv')[,2]) |>
    group_by(combination_id) |>
    reframe(read_csv(path,
            col_names = c('section', 'atom_id', 'property', 'value'),
            col_types = cols(.default = col_character(), value = col_double()))) |>
    # AIMAll seems to maintain atom order but rename them, one-indexed in all
    # caps. To match my own naming scheme I need to capitalize only the first
    # letter
    mutate(atom_id = str_to_title(atom_id)) |>
    # Filter out simulations with inaccurate integration
    inner_join(accurate_integration, by = 'combination_id')

write_csv(one_tbl, 'combined_one_tbl.csv.gz')

# For IQA interaction energy I also need the two atom tables
two_tbl <- tibble(path = Sys.glob('aimall_tbl/*_twoatom.csv')) |>
    mutate(combination_id = str_match(path, 'aimall_tbl/(.*)_twoatom.csv')[,2]) |>
    group_by(combination_id) |>
    reframe(read_csv(path,
            col_names = c('section', 'atom_a', 'atom_b', 'property', 'value'),
            col_types = cols(.default = col_character(), value = col_double()))) |>
    mutate(atom_a = str_to_title(atom_a),
           atom_b = str_to_title(atom_b)) |>
    inner_join(accurate_integration, by = 'combination_id')

write_csv(two_tbl, 'combined_two_tbl.csv.gz')

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
    molecule_charge = col_integer()
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
# attractors, for filtering with a join. Should also exclude simulations
# without non-nuclear attractors, if the field is more extreme than a
# simulation with a non-nuclear attractor
# This filtering technique doesn't work if there's NNA's without a field
directional_field_nna <- nna_number |>
    # Calculate bounds for the field
    filter(field_value != 0) |>
    mutate(direction = sign(field_value),
           strength = abs(field_value))
nna_at_field_or_weaker <- directional_field_nna |>
    left_join(directional_field_nna,
              by = c('formula', 'direction'),
              suffix = c('', '_other'),
              relationship = 'many-to-many') |>
    # What's the strongest field under the constraint that it has no NNA's,
    # and no weaker field has any NNA?
    filter(strength >= strength_other) |>
    group_by(combination_id) |>
    summarize(nna_at_field_or_weaker = sum(nna_number_other != 0),
              .groups = 'drop')
write_csv(nna_at_field_or_weaker, 'nna_at_field_or_weaker.csv.gz')

too_strong <- nna_at_field_or_weaker |>
    filter(nna_at_field_or_weaker > 0) |>
    select(combination_id)

# Bader charge of each atom
bader_charge_atom <- one_tbl |>
    filter(section == 'Some Atomic Properties' & property == 'q(A)') |>
    group_by(combination_id, atom_id) |>
    transmute(charge = value) |>
    ungroup()

# Becke charge of each atom
becke_charge_atom <- read_csv('raw_becke_populations.csv.gz', col_types = cols(
    combination_id = col_character(),
    atom_number = col_double(),
    symbol = col_character(),
    excess_electrons = col_double(),
    population = col_double(),
    net_spin = col_double()
)) |>
    # Concatenate element symbol and number to obtain atom id
    # Numbers are one-indexed just like in AIMAll output so it should match
    # Element symbols are already in title case so there's no need to convert
    mutate(atom_id = glue('{symbol}{atom_number}')) |>
    # Convert excess electron populations to charge (negative excess population)
    # Actually that leads to a reversal relative to Bader charge
    mutate(charge = excess_electrons) |>
    # Select down to three columns like in the Bader charge table
    select(combination_id, atom_id, charge)

# Make a combined table with each kind of partial charge
charge_atom <- bind_rows(bader = bader_charge_atom, becke = becke_charge_atom,
                         .id = 'charge_type')
write_csv(charge_atom, 'charge_atom.csv.gz')

# Bader charge of each group
bader_charge_group <- charge_atom |>
    anti_join(too_strong, by = 'combination_id') |>
    left_join(simulation_table, by = 'combination_id') |>
    left_join(coordinates, by = c('formula', 'atom_id'),
              relationship = 'many-to-one') |>
    group_by(combination_id, formula, charge_type, donor_or_acceptor) |>
    summarize(
        heavy_atom_symbol = symbol[symbol != 'H'][1],
        total_charge = sum(charge),
        .groups = 'drop'
    ) |>
    rename(symbol = heavy_atom_symbol) |>
    # I previously used the column name "group_bader_charge"
    # To minimize changes, I'm going to continue using that column name
    mutate(column_name = glue('group_{charge_type}_charge')) |>
    select(-charge_type) |>
    pivot_wider(names_from = column_name, values_from = total_charge)

write_csv(bader_charge_group, 'bader_charge_group.csv.gz')

# The total energies of each molecule. These don't seem to be available
# directly, so I obtain them by summing the atomic energies
TO_EV <- 27.211386246
combined_atom_energies <- tibble(path = Sys.glob('aimall_tbl/*_oneatom.csv')) |>
    mutate(combination_id = str_match(path, 'aimall_tbl/(.*)_oneatom.csv')[,2]) |>
    anti_join(too_strong, by = 'combination_id') |>
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

# IQA energy of each atom
# Repetitive, same as previous bader charge code, should probably be combined
atom_iqa <- one_tbl |>
    filter(section == "IQA Intraatomic (\"Self\") Energy Components" &
                 property == "E_IQA_Intra(A)") |>
    group_by(combination_id, atom_id) |>
    transmute(iqa_energy = value * TO_EV) |>
    ungroup()
write_csv(atom_iqa, "atom_iqa.csv.gz")

group_iqa <- atom_iqa |>
    left_join(simulation_table, by = 'combination_id') |>
    left_join(coordinates, by = c('formula', 'atom_id'),
              relationship = 'many-to-one') |>
    group_by(combination_id, formula, donor_or_acceptor) |>
    summarize(
        heavy_atom_symbol = symbol[symbol != 'H'][1],
        total_iqa_energy = sum(iqa_energy),
        .groups = 'drop'
    ) |>
    rename(symbol = heavy_atom_symbol)
write_csv(group_iqa, "group_iqa.csv.gz")

# IQA interaction energy
atom_interaction <- two_tbl |>
    filter(section == 'IQA Diatomic "Interaction" Energy Components' &
           property == 'E_IQA_Inter(A,B)') |>
    group_by(combination_id, atom_a, atom_b) |>
    transmute(iqa_energy = value * TO_EV) |>
    ungroup()
write_csv(atom_interaction, 'atom_interaction.csv.gz')

group_interaction <- atom_interaction |>
    left_join(simulation_table, by = 'combination_id') |>
    left_join(coordinates, by = c('formula', 'atom_a' = 'atom_id'),
              relationship = 'many-to-one') |>
    left_join(coordinates, by = c('formula', 'atom_b' = 'atom_id'),
              relationship = 'many-to-one', suffix = c('_a', '_b')) |>
    filter(donor_or_acceptor_a != donor_or_acceptor_b) |>
    group_by(combination_id, formula, donor_or_acceptor_a, donor_or_acceptor_b) |>
    summarize(
        heavy_atom_symbol_a = symbol_a[symbol_a != 'H'][1],
        heavy_atom_symbol_b = symbol_b[symbol_b != 'H'][1],
        total_interaction_energy = sum(iqa_energy),
        .groups = 'drop'
    ) |>
    rename(symbol_a = heavy_atom_symbol_a,
           symbol_b = heavy_atom_symbol_b)
write_csv(group_interaction, 'group_interaction.csv.gz')

# Combine charge and energy into a single table mapping charge to energy
charge_energy <- coordinates |>
    # Heavy atoms only
    filter(symbol != 'H') |>
    select(formula, symbol, donor_or_acceptor) |>
    # Populate with the combination id of each simulation
    # Every simulation should match two rows in this table corresponding to the
    # two heavy atoms
    # And every atom should match multiple simulations corresponding to each
    # field strength
    left_join(simulation_table, by = 'formula',
              relationship = 'many-to-many') |>
    # Join with the charge of the donor and acceptor groups (as separate rows)
    left_join(bader_charge_group,
              by = c('formula', 'combination_id', 'symbol',
                     'donor_or_acceptor'),
              relationship = 'one-to-one') |>
    # Join with the total energy of the combination
    # This will produce some redundancy, since it will be the same for both
    # donor and acceptor
    left_join(total_atom_energies,
              by = 'combination_id',
              relationship = 'many-to-one') |>
    rename(total_energy = energy) |>
    # Join with the IQA energy of each individual group
    left_join(group_iqa,
              by = c('formula', 'combination_id', 'symbol',
                     'donor_or_acceptor')) |>
    rename(iqa_group_energy = total_iqa_energy) |>
    # Join with the IQA interaction energy between the groups
    # Again, this will be redundant. The interaction energy is the same on both
    # "sides"
    left_join(select(group_interaction, formula, combination_id,
                     total_interaction_energy),
              by = c('formula', 'combination_id'),
              relationship = 'many-to-one') |>
    rename(iqa_interaction_energy = total_interaction_energy)


# Probably I could have made a table labeling donor and acceptor by heavy atom
# symbol earlier and then used that information for the Bader charge table,
# but instead I did it as part of generating the Bader charge table, so
# I'll just re-use that information. Could refactor later
heavy_atom_element <- bader_charge_group |>
    select(formula, symbol, donor_or_acceptor) |>
    distinct()
# I want a table that you can join with, and based on the atom and the
# formula, it will fill in the symbol of the other heavy atom in the molecule
other_atom_symbol <- heavy_atom_element |>
    rename(other_symbol = symbol) |>
    mutate(donor_or_acceptor = ifelse(donor_or_acceptor == 'donor',
                                      'acceptor', 'donor'))

charge_energy_morelabels <- charge_energy |>
    # Add it into the charge and energy table, it's a useful thing to have
    left_join(other_atom_symbol, by = c('formula', 'donor_or_acceptor')) |>
    # Also a combined "group_id" column since I use it for grouping in some of
    # the plots
    mutate(group_id = glue('{formula}:{symbol}:{donor_or_acceptor}'))

write_csv(charge_energy_morelabels, 'charge_energy.csv.gz')
