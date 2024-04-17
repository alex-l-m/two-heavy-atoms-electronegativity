library(readr)
library(dplyr)

# Atomic number of each element
number <- read_csv('wolframalpha_number.csv', col_names = c('1', '2', 'element', '3', 'atomic_number', '4'), col_types = cols(element = col_character(), atomic_number = col_integer(), .default = col_character())) |>
    select(element, atomic_number)

# One or two letter symbol of the element, corresponding to the full element string
symbol <- read_csv('wolframalpha_symbol.csv', col_names = c('element', '1', 'symbol', '2')) |>
    select(symbol, element)

# https://www.wolframalpha.com/input?i=coulomb%27s+constant+in+eV+*+angstrom+%2F+%28electron+charge%29%5E2
coulomb <- 14.39964548

# https://www.wolframalpha.com/input?i=convert+1+kJ%2Fmol+to+eV
kJ2eV <- 0.01036

# Table of electron affinity of each element
# https://www.wolframalpha.com/input?i=electron+affinity+of+%7Bhydrogen%2C+carbon%2C+nitrogen%2C+oxygen%2C+fluorine%2C+silicon%2C+phosphorus%2C+sulfur%2C+chlorine%7D
wa_ea <- read_csv('wolframalpha_ea.csv', col_names = c('rownum', '?', 'element', '?', 'electron_affinity', 'unit')) |>
    select(element, electron_affinity) |>
    mutate(electron_affinity = kJ2eV * electron_affinity)

# Table of ionization potential of each element
# https://www.wolframalpha.com/input?i=first+ionization+potential+of+%7Bhydrogen%2C+carbon%2C+nitrogen%2C+oxygen%2C+fluorine%2C+silicon%2C+phosphorus%2C+sulfur%2C+chlorine%7D
wa_ip <- read_csv('wolframalpha_ip.csv', col_names = c('rownum', '?', 'element', '?', 'ionization_potential', 'unit')) |>
    select(element, ionization_potential) |>
    mutate(ionization_potential = kJ2eV * ionization_potential)

# Join all tables of element properties into a single table
element_properties <- number |>
    left_join(symbol, by = 'element') |>
    left_join(wa_ea, by = 'element') |>
    left_join(wa_ip, by = 'element') |>
    mutate(
      mulliken_electronegativity = (electron_affinity + ionization_potential) / 2,
      mulliken_hardness = -(electron_affinity - ionization_potential)
    )

# Save as a csv file
write_csv(element_properties, 'element_properties.csv.gz')
