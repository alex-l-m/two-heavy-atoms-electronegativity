# Plot interaction terms as a function of scale factor
library(tidyverse)
library(glue)
library(tidymodels)
library(ggrepel)
library(cowplot)
library(ggdark)

theme_set(dark_mode(theme_cowplot(font_size = 12)))

charge_energy <- read_csv('charge_energy.csv.gz', col_types = cols(
    combination_id = col_character(),
    symbol = col_character(),
    other_symbol = col_character(),
    donor_or_acceptor = col_character(),
    charge = col_double(),
    bader_charge = col_double(),
    cp2k_hirshfeld_charge = col_double(),
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
    electronegativity_field_discrete = col_double()
))

# Table of scale numbers and scale factors
# "scale_number" is an integer that indexes the scale. It is ordered, but
# doesn't otherwise have any guaranteed relationship to the lattice constant
# "scale" is the scaling factor: the ratio of the lattice constant to the
# equilibriums lattice constant. This is a physically meaningful quantity,
# measured in angstroms
scales <- charge_energy |>
    select(formula, scale_number, scale) |>
    distinct()

regression_estimate_coltypes <- cols(
    term = col_character(),
    estimate = col_double(),
    std.error = col_double(),
    statistic = col_double(),
    p.value = col_double()
)

regression_estimates <- read_csv('3-5:zincblende_electronegativity_regression_estimates.csv', col_types = regression_estimate_coltypes)

# Regular expression for recognizing interaction terms and extracting the formula and scale number
# Example: interaction_GaP_S3
term_name_regex = 'interaction_([A-Za-z]*)_S([0-9]+)'

# Table for plotting interaction parameter as a function of scale factor
# Needs to have columns for formula and scale factor, which can be obtained by
# extracting formula and scale number from the name of the term and joining to
# the table of scale factors
interaction_parameter_tbl <- regression_estimates |>
    filter(str_detect(term, term_name_regex)) |>
    mutate(formula = str_extract(term, term_name_regex, 1),
           scale_number = as.integer(str_extract(term, term_name_regex, 2))) |>
    select(formula, scale_number, estimate) |>
    left_join(scales, by = c('formula', 'scale_number'))

interaction_parameter_plot <- interaction_parameter_tbl |>
    ggplot(aes(x = scale, y = estimate)) +
    geom_point() +
    geom_line() +
    facet_wrap(~formula, scales = 'free_y') +
    labs(x = 'Scale factor', y = 'Interaction parameter')

ggsave('3-5:zincblende_interaction_parameter_plot.png',
       interaction_parameter_plot,
       height = unit(4.76, 'in'), width = unit(5.67, 'in'))
