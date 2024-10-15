library(tidyverse)
library(cowplot)
library(ggdark)
library(robustbase)
theme_set(dark_mode(theme_cowplot()))
library(glue)
library(ggrepel)
library(broom)

# Raw energies and derivatives
charge_energy <- read_csv('charge_energy.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    charge = col_double(),
    bader_charge = col_double(),
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
    formula = col_character(),
    electronegativity_field_discrete = col_double()
))

# Compare Bader charges with my own charge (the "charge" column)
comparison_plot <- charge_energy |>
    filter(category == '3-5' & crystal_structure == 'zincblende') |>
    filter(donor_or_acceptor == 'acceptor') |>
    ggplot(aes(x = charge, y = bader_charge, color = formula)) +
    geom_line() +
    geom_point() +
    coord_fixed()

ggsave('charge_comparison.png', comparison_plot, width = unit(7, 'in'), height = unit(7, 'in'))
