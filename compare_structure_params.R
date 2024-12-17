# Compare the electronegativity, hardness, and interaction parameters for
# different structures
library(tidyverse)
library(glue)
library(tidymodels)
library(ggrepel)
library(cowplot)
library(ggdark)

theme_set(dark_mode(theme_cowplot(font_size = 12)))

combination <- '3-5'
reference_element <- read_csv('3-5:zincblende_reference_element.csv', col_types = cols(symbol = col_character()))
reference_formulas <- read_csv('3-5:zincblende_reference_formulas.csv', col_types = cols(formula = col_character()))

infiles <- Sys.glob(glue('{combination}:*_electronegativity_regression_estimates.csv'))
regression_estimate_coltypes <- cols(
    term = col_character(),
    estimate = col_double(),
    std.error = col_double(),
    statistic = col_double(),
    p.value = col_double()
)
filename_regex <- '(.*):(.*)_electronegativity_regression_estimates.csv'
term_regex <- '([^_]*)_(.*)'
regression_estimates <- tibble(infile = infiles) |>
    group_by(infile) |>
    reframe(read_csv(infile, col_types = regression_estimate_coltypes)) |>
    # Extract information from the filename
    mutate(combination = str_extract(infile, filename_regex, 1),
           structure = str_extract(infile, filename_regex, 2)) |>
    # Extract information from the term name
    mutate(variable_type = str_extract(term, term_regex, 1),
           entity = str_extract(term, term_regex, 2)) |>
    # Hack: only plot scale zero
    filter(str_detect(entity, '_S0'))

# Make a plot comparing estimates the two structures, with separate scatter
# plots for each type of variable
comparison_table <- regression_estimates |>
    select(structure, variable_type, entity, estimate) |>
    pivot_wider(names_from = structure, values_from = estimate) |>
    bind_rows(mutate(rename(reference_element, entity = symbol), variable_type = 'electronegativity', zincblende = 0, rocksalt = 0)) |>
    bind_rows(mutate(rename(reference_element, entity = symbol), variable_type = 'hardness', zincblende = 0, rocksalt = 0)) |>
    bind_rows(mutate(rename(reference_formulas, entity = formula), variable_type = 'interaction', zincblende = 0, rocksalt = 0))
comparison_plot <- comparison_table |>
    ggplot(aes(x = zincblende, y = rocksalt)) +
    facet_wrap(~variable_type, scales = 'free') +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    geom_text_repel(aes(label = entity))

ggsave(glue('{combination}_electronegativity_comparison_plot.png'),
       comparison_plot,
       height = unit(4.76, 'in'), width = unit(5.67, 'in'))
