# Starting with a table of atomic numbers, make a table where each row
# represents simulation of a single atom
library(tidyverse)

atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_double()
))

simulation_metadata <- atomic_numbers |>
    # For now, all simulations are neutral
    mutate(charge = 0) |>
    # Multiplicity, based on the atomic number and the charge
    mutate(n_electrons = atomic_number - charge,
           multiplicity = 1 + n_electrons %% 2)

write_csv(simulation_metadata, 'single_atom_simulations.csv')
