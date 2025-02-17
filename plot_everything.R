library(tidyverse)
library(glue)
library(cowplot)
library(ggdark)
library(robustbase)
theme_set(dark_mode(theme_cowplot()))
library(glue)
library(ggrepel)
library(broom)

TO_EV <- 27.211386246

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
    field_number = col_double(),
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
    electronegativity_field_analytic = col_double(),
    unscaled_structure_id = col_character(),
    scale_number = col_integer(),
    scale = col_double()
))

everything_tbl <- charge_energy |>
    # Acceptor only
    filter(donor_or_acceptor == 'acceptor') |>
    select(
        # Identify a simulation
        combination_id, scale_number,
        # Cation and anion, for faceting
        symbol_cation, symbol_anion,
        # x axis
        field_number,
        # Everything I want to plot
        charge, field_value, total_energy,
        first_iteration_charge, last_iteration_charge,
        electronegativity_field_discrete,
        electronegativity_field_analytic) |>
    # Change total energy so that it's referenced to the minimum value
    group_by(symbol_cation, symbol_anion, scale_number) |>
    mutate(total_energy = total_energy - min(total_energy)) |>
    ungroup() |>
    # Since I'm plotting as a function of field number, I have to sort by field
    # number before taking differences
    arrange(symbol_cation, symbol_anion, scale_number, field_number) |>
    # Derived columns
    group_by(symbol_cation, symbol_anion, scale_number) |>
    mutate(field_value_ev = field_value * TO_EV,
           dE = c(NA, diff(total_energy)),
           dq = c(NA, diff(charge)),
           dq0 = first_iteration_charge - lag(last_iteration_charge),
           dq0pdq = dq0 / dq,
           dV = c(NA, diff(field_value_ev)),
           correction_factor = electronegativity_field_discrete / field_value_ev) |>
    ungroup() |>
    # Category column based on the rounded charge
    mutate(charge_category = as.factor(ceiling(charge))) |>
    # Pivot so everything I want the plot is one column
    pivot_longer(-c(combination_id, scale_number,
                    symbol_cation, symbol_anion,
                    field_number, charge_category),
                 names_to = 'y', values_to = 'y_value')

# Make an output directory for all the plots
output_dir <- 'plot_everything'
dir.create(output_dir, showWarnings = FALSE)
for (this_y in unique(everything_tbl$y))
{
    plt <- everything_tbl |>
        filter(y == this_y) |>
        ggplot(aes(x = field_number, y = y_value, color = charge_category)) +
        facet_grid(symbol_cation ~ symbol_anion) +
        geom_line() +
        geom_point(size = 0.5) + 
        labs(x = 'Field number',
             y = this_y) +
        ggtitle(this_y)
    ggsave(glue('{output_dir}/{this_y}.png'), plt,
           width = unit(11.5, 'in'), height = unit(4.76, 'in'))
}

