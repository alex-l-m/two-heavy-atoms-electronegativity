library(tidyverse)
library(cowplot)
library(ggrepel)
library(tidymodels)

theme_set(theme_cowplot(font_size = 24))

charge_energy <- read_csv('charge_energy.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    charge = col_double(),
    bader_charge = col_double(),
    cp2k_hirshfeld_charge = col_double(),
    weight_function_derivative = col_double(),
    structure_id = col_character(),
    field_number = col_integer(),
    field_value = col_double(),
    total_energy = col_double(),
    cdft = col_logical(),
    group_id = col_character(),
    category = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character(),
    crystal_structure = col_character(),
    scale_number = col_integer(),
    scale = col_double(),
    unscaled_structure_id = col_character(),
    formula = col_character(),
    first_iteration_charge = col_double(),
    last_iteration_charge = col_double(),
    lagged_first_iteration_charge = col_double(),
    lagged_last_iteration_charge = col_double(),
    electronegativity_field_discrete = col_double(),
    electronegativity_field_analytic = col_double()
))

# Charge of the acceptor atoms in ground state structures
acceptor_charge <- charge_energy |>
    filter(donor_or_acceptor == 'acceptor' & scale_number == 0 & !cdft) |>
    select(formula, charge) |>
    rename(acceptor_charge = charge)

loo <- read_csv('3-5:zincblende_loo.csv.gz', col_types = cols(
    formula = col_character(),
    electronegativity_donor = col_double(),
    electronegativity_acceptor = col_double(),
    hardness_donor = col_double(),
    hardness_acceptor = col_double(),
    eeq_acceptor_charge = col_double()
))

loo_pred <- loo |>
    select(formula, eeq_acceptor_charge)

comparison_table <- acceptor_charge |>
    left_join(loo_pred, by = 'formula')

comparison_plot <- comparison_table |>
    ggplot(aes(x = acceptor_charge, y = eeq_acceptor_charge, 
               label = formula)) +
    geom_point() +
    geom_label_repel() +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'red') +
    theme_cowplot(font_size = 24) +
    coord_obs_pred()

ggsave('3-5:zincblende_loo_comparison_plot.png', comparison_plot,
       height = unit(4.76, 'in'), width = unit(5.67, 'in'))
