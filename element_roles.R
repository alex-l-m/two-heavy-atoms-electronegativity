library(tidyverse)
library(glue)

group_numbers <- read_csv('group_numbers.csv', col_types = cols(
    symbol = col_character(),
    group = col_integer()
))

group_numbers |>
    group_by(symbol) |>
    transmute(
        role = ifelse(
            group < 4, 'cation', ifelse(
            group > 4, 'anion', ifelse(
            group == 4, 'neutral'))),
        category = glue('{min(group, 8-group)}-{max(group, 8-group)}')
        ) |>
    write_csv('element_roles.csv')
