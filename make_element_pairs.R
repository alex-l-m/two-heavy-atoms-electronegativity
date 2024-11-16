library(readr)
library(dplyr)

element_roles <- read_csv('element_roles.csv', col_types = cols(
    symbol = col_character(),
    role = col_character(),
    category = col_character()))

element_roles |>
    left_join(element_roles,
              by = 'category',
              suffix = c('_cation', '_anion'),
              relationship = 'many-to-many') |>
    # This table has "switched" rows, where the cation symbol is in
    # "symbol_anion", etc. This shouldn't matter because in make_cp2k_jobs.py,
    # they will be skipped when there's no corresponding CONTCAR file. In fact,
    # grouping by category was not really necessary for this reason, but I want
    # to have category in a table so I can divide up structures for plots and
    # regressions later.
    select(category, symbol_cation, symbol_anion) |>
    write_csv('element_pairs.csv')
