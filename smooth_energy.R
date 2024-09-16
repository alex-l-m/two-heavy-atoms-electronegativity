# Load energy as a function of charge, and then smooth it in order to take a
# derivative and obtain electronegativity
library(tidyverse)

charge_energy <- read_csv('charge_energy.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    charge = col_double(),
    structure_id = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    total_energy = col_double(),
    cdft = col_logical(),
    group_id = col_character(),
    category = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character(),
    crystal_structure = col_character(),
    formula = col_character()
))

# Make a table for smoothing energy
# It should have all the energies in the same column, so that I can smooth that
# column
# It should also be just the ones where I did constrained DFT
charge_energy_long <- charge_energy |>
    # I used to have IQA energies, those would be the other energy types
    pivot_longer(c(total_energy),
                 names_to = 'energy_type', values_to = 'energy')
energy_for_smoothing <- charge_energy_long |>
    add_count(structure_id) |>
    filter(n > 3) |>
    select(-n)
# Table of smoothed energies, sampled on a grid, and derivative
smoothed_energy_derivatives_long <- energy_for_smoothing |>
    group_by(
            # Identification of the system
            crystal_structure, category, structure_id, symbol_anion, symbol_cation, formula,
            # Identification of the atom
            symbol, other_symbol, donor_or_acceptor, group_id,
            # Which energy calculation is being smoothed
            energy_type) |>
    reframe(tibble(
            new_charge = seq(min(charge), max(charge),
                             length.out = 100),
            # Simple linear interpolation
            energy = approx(charge, energy, xout = new_charge)$y)
    ) |>
    rename(charge = new_charge) |>
    # Reframe loses the groups, but I need those exact same groups to calculate
    # the derivative
    group_by(crystal_structure, category, structure_id, symbol_anion, symbol_cation, formula,
             symbol, other_symbol, donor_or_acceptor, group_id,
             energy_type) |>
    # Calculate the derivative as discrete differences
    arrange(charge) |>
    mutate(derivative = c(NA, diff(energy) / diff(charge))) |>
    ungroup()

# I never have any need to plot both the energy and the derivative, I don't
# think. Therefore, it may make more sense to have separate tables of energy
# and derivative, with a column for each energy type. I don't need these in the
# same plot either, but then I can specify what energy type I'm plotting
# with the aesthetic. This requires separate plotting code for each energy
# type, but that's appropriate anyway because the data will be grouped
# differently
smoothed_energy <- smoothed_energy_derivatives_long |>
    select(-derivative) |>
    pivot_wider(names_from = energy_type, values_from = energy)
write_csv(smoothed_energy, 'smoothed_energy.csv.gz')
smoothed_energy_derivatives <- smoothed_energy_derivatives_long |>
    select(-energy) |>
    pivot_wider(names_from = energy_type, values_from = derivative)
write_csv(smoothed_energy_derivatives, 'smoothed_energy_derivatives.csv.gz')

# Assign derivatives to the unsmoothed energies
energy_derivatives <- smoothed_energy_derivatives_long |>
    # I'm going to join a table that also has a charge column, so I need a
    # separate name for the smoothed charges
    rename(smoothed_charge = charge) |>
    # Nesting is necessary because I want a single row, containing the single
    # Bader charge of this combination. It's the only way to handle the
    # different sizes of the smoothed data I'm interpolating, and the piont I'm
    # interpolating on
    group_by(crystal_structure, category, structure_id, symbol_anion, symbol_cation, formula,
             symbol, other_symbol, donor_or_acceptor, group_id,
             energy_type) |>
    nest() |>
    ungroup() |>
    left_join(charge_energy_long,
              by = c('crystal_structure', 'category', 'structure_id',
                     'symbol_anion', 'symbol_cation', 'formula',
                     'symbol', 'other_symbol', 
                     'donor_or_acceptor', 'group_id', 'energy_type')) |>
    group_by(
        # Have to map charge to derivative for each atom individually
        combination_id, symbol, other_symbol, donor_or_acceptor, group_id,
        charge, energy_type,
        # But I also want to keep all of the simulation metadata so I don't
        # have to join with it later
        crystal_structure, category, structure_id, symbol_anion, symbol_cation, formula,
        cdft, field_number, field_value
    ) |>
    # Using summarize rather than mutate just to get rid of the nested "data"
    # column
    summarize(derivative = sapply(data, function(df) {with(df, 
        approx(smoothed_charge, derivative, xout = charge)$y)
    }), .groups = 'drop') |>
    # Pivot wider like the energy table
    pivot_wider(
        names_from = energy_type,
        values_from = derivative
    )
write_csv(energy_derivatives, 'energy_derivatives.csv.gz')
