# Make a table of parameters related to softness

library(tidyverse)

simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer(),
    gamess_input_file = col_character()
))

total_atom_energies <- read_csv('total_atom_energies.csv.gz', col_types = cols(
    combination_id = col_character(),
    energy = col_double()
))

# Copying and pasting which is no good. Should re-use the one from
# smooth_energy but I can't because it's filtered
bader_charge <- read_csv('bader_charge_group.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    total_bader_charge = col_double()
))

total_atom_energies <- read_csv('total_atom_energies.csv.gz', col_types = cols(
    combination_id = col_character(),
    energy = col_double()
))

charge_transfer <- bader_charge |>
    rename(charge = total_bader_charge) |>
    # Column 'donor_or_acceptor' has values 'donor' and 'acceptor'
    # Pivot the 'charge' and 'symbol' columns so we have columns
    # 'charge_donor', 'charge_acceptor', 'symbol_donor', 'symbol_acceptor'
    pivot_wider(names_from = donor_or_acceptor,
                values_from = c(charge, symbol))

# Make a table of energies and charges of the no-field simulations, including
# ions
nofield_charge_energy <- simulation_table |>
    filter(abs(field_value) < 1e-6) |>
    left_join(total_atom_energies, by = 'combination_id') |>
    left_join(charge_transfer, by = c('combination_id', 'formula')) |>
    select(formula, molecule_charge, energy,
           charge_acceptor, charge_donor) |>
    mutate(ion = ifelse(molecule_charge == -1, 'negative', ifelse(molecule_charge == 1, 'positive', ifelse(molecule_charge == 0, 'neutral', NA)))) |>
    arrange(formula, molecule_charge)

# Ionization potentials and electron affinities
ipea <- nofield_charge_energy |>
    pivot_wider(id_cols = formula, names_from = ion, values_from = energy) |>
    group_by(formula) |>
    transmute(ip = positive - neutral, ea = neutral - negative,
              electronegativity = (ip + ea) / 2,
              hardness = ip - ea) |>
    ungroup()
                              
# Upper and lower fukui function
fukui <- nofield_charge_energy |>
    pivot_wider(id_cols = formula, names_from = ion, values_from = c(charge_acceptor, charge_donor)) |>
    group_by(formula) |>
    transmute(
            lower_fukui_acceptor = charge_acceptor_positive - charge_acceptor_neutral,
            lower_fukui_donor = charge_donor_positive - charge_donor_neutral,
            upper_fukui_acceptor = charge_acceptor_neutral - charge_acceptor_negative,
            upper_fukui_donor = charge_donor_neutral - charge_donor_negative) |>
    ungroup()

ipea_fukui <- ipea |>
    left_join(fukui, by = 'formula')

write_csv(ipea_fukui, 'ipea_fukui.csv.gz')
