library(tidyverse)

simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer(),
    gamess_input_file = col_character()
))

neutral_combinations <- simulation_table |>
    filter(molecule_charge == 0) |>
    select(combination_id) |>
    distinct()

bader_charge <- read_csv("bader_charge_group.csv.gz", col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    total_bader_charge = col_double()
))

total_atom_energies <- read_csv("total_atom_energies.csv.gz", col_types = cols(
    combination_id = col_character(),
    energy = col_double()
))

charge_transfer <- bader_charge |>
    inner_join(neutral_combinations, by = "combination_id") |>
    rename(charge = total_bader_charge) |>
    # Column "donor_or_acceptor" has values "donor" and "acceptor"
    # Pivot the "charge" and "symbol" columns so we have columns
    # "charge_donor", "charge_acceptor", "symbol_donor", "symbol_acceptor"
    pivot_wider(names_from = donor_or_acceptor,
                values_from = c(charge, symbol))

write_csv(charge_transfer, "charge_transfer.csv.gz")

# Smooth energies
energy_charge <- total_atom_energies |>
    # Inner join so we filter according to how the charge transfer table has
    # been filtered
    inner_join(charge_transfer, by = "combination_id")

write_csv(energy_charge, "energy_charge.csv.gz")

smoothed_energy <- energy_charge |>
    group_by(formula, symbol_donor, symbol_acceptor) |>
    reframe(tibble(new_charge = seq(min(charge_acceptor), max(charge_acceptor),
                                    length.out = 100),
                   energy = predict(loess(energy ~ charge_acceptor),
                                    newdata = tibble(charge_acceptor = new_charge)))) |>
    ungroup() |>
    rename(charge_transfer = new_charge) |>
    group_by(formula) |>
    mutate(derivative = c(NA, diff(energy) / diff(charge_transfer)))

write_csv(smoothed_energy, "smoothed_energy.csv.gz")

# Assign derivatives to the unsmoothed energies
energy_derivatives <- smoothed_energy |>
    select(formula, charge_transfer, derivative) |>
    group_by(formula) |>
    nest() |>
    left_join(energy_charge, by = "formula") |>
    # Very hard to read. Relies on the fact that "charge_acceptor" has been
    # renamed "charge_transfer" during smoothing, and I shouldn't be grouping
    # by a value
    group_by(combination_id, charge_acceptor) |>
    summarize(derivative = sapply(data, function(df) {with(df, 
        approx(charge_transfer, derivative, xout = charge_acceptor)$y)})) |>
    rename(charge_transfer = charge_acceptor)
    
write_csv(energy_derivatives, "energy_derivatives.csv.gz")
    
    
