# Starting with a table of atomic numbers, make a table where each row
# represents simulation of a single atom
library(tidyverse)
library(glue)

atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_double()
))

simulation_metadata <- atomic_numbers |>
    # For every element, charges -2 to 2
    expand_grid(charge = -2:2) |>
    # Multiplicity, based on the atomic number and the charge
    mutate(n_electrons = atomic_number - charge,
           multiplicity = 1 + n_electrons %% 2,
           job_id = glue('{symbol}_{charge}'))

write_csv(simulation_metadata, 'single_atom_simulations.csv')
