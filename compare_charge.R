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
    cp2k_hirshfeld_charge = col_double(),
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

# Atomic numbers, for ordering
atomic_numbers <- read_csv('atomic_numbers.csv', col_types = cols(
    symbol = col_character(),
    atomic_number = col_double()
))
formula_order <- charge_energy |>
    filter(category == '3-5') |>
    select(formula, symbol_cation, symbol_anion) |>
    # Order based on atomic numbers of the elements
    left_join(atomic_numbers, by = c('symbol_cation' = 'symbol')) |>
    left_join(atomic_numbers, by = c('symbol_anion' = 'symbol'),
                              suffix = c('_cation', '_anion')) |>
    arrange(atomic_number_cation, atomic_number_anion) |>
    distinct(formula) |>
    pull(formula)
# Compare charges with my own charge (the "charge" column)
comparison_tbl <- charge_energy |>
    filter(category == '3-5') |>
    mutate(formula = factor(formula, levels = formula_order)) |>
    filter(donor_or_acceptor == 'acceptor') |>
    # Turn the Bader and Hirshfeld charges into a single column
    pivot_longer(cols = c(bader_charge, cp2k_hirshfeld_charge),
                 names_to = 'charge_type', values_to = 'other_charge') |>
    # Remove the suffix to make the names of the charge types shorter
    mutate(charge_type = str_remove(charge_type, '_charge'))
comparison_plot <- comparison_tbl |>
    ggplot(aes(x = charge, y = other_charge, color = crystal_structure)) +
    facet_grid(charge_type ~ formula) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    geom_line() +
    # Points for field value 0
    geom_point(data = filter(comparison_tbl, field_value == 0)) +
    coord_fixed() +
    theme(legend.position = 'bottom')

ggsave('charge_comparison.png', comparison_plot, 
       width = unit(11.5, 'in'), height = unit(4.76, 'in'))
