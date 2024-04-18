library(tidyverse)

simulation_table <- read_csv('simulations.csv.gz', col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    molecule_charge = col_integer()
))

neutral_combinations <- simulation_table |>
    filter(molecule_charge == 0) |>
    select(combination_id) |>
    distinct()

# Not really total atom energies but just copying and pasting for now
total_atom_energies <- read_csv("group_iqa.csv.gz", col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    total_iqa_energy = col_double()
)) |>
    rename(energy = total_iqa_energy) |>
    select(combination_id, donor_or_acceptor, symbol, energy)

bader_charge <- read_csv("bader_charge_group.csv.gz", col_types = cols(
    combination_id = col_character(),
    formula = col_character(),
    donor_or_acceptor = col_character(),
    symbol = col_character(),
    total_bader_charge = col_double()
))

# Smooth energies
energy_charge <- total_atom_energies |>
    # Inner join so we filter according to how the charge transfer table has
    # been filtered
    inner_join(bader_charge, by = c('combination_id', 'symbol', 'donor_or_acceptor')) |>
    filter(-1 < total_bader_charge & total_bader_charge < 1)

write_csv(energy_charge, "iqa_energy_charge.csv.gz")

smoothed_energy <- energy_charge |>
    group_by(formula, donor_or_acceptor, symbol) |>
    reframe(tibble(new_charge = seq(min(total_bader_charge), max(total_bader_charge),
                                    length.out = 100),
                   energy = predict(loess(energy ~ total_bader_charge),
                                    newdata = tibble(total_bader_charge = new_charge)))) |>
    group_by(formula, donor_or_acceptor, symbol) |>
    mutate(derivative = c(NA, diff(energy) / diff(new_charge))) |>
    ungroup() |>
    rename(charge = new_charge)

write_csv(smoothed_energy, "iqa_smoothed_energy.csv.gz")

# Assign derivatives to the unsmoothed energies
energy_derivatives <- smoothed_energy |>
    select(formula, donor_or_acceptor, charge, derivative) |>
    group_by(formula, donor_or_acceptor) |>
    nest() |>
    left_join(energy_charge, by = c('formula', 'donor_or_acceptor')) |>
    # Very hard to read. Relies on renaming charges
    group_by(combination_id, donor_or_acceptor, symbol) |>
    # Too repetitive
    summarize(
        newer_charge = sapply(data, function(df) {with(df, 
        approx(charge, derivative, xout = total_bader_charge)$x)}),
        derivative = sapply(data, function(df) {with(df, 
        approx(charge, derivative, xout = total_bader_charge)$y)})) |>
    rename(charge = newer_charge)

write_csv(energy_derivatives, "iqa_energy_derivatives.csv.gz")
