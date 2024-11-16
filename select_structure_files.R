# Use the list of element pairs of the selected elements to select structure
# files containing those elements and make a table where each row corresponds
# to a system to simulate

library(readr)
library(dplyr)
library(glue)

element_pairs <- read_csv('element_pairs.csv', col_types = cols(
    category = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character()
))

structure_files <- read_csv('structure_files.csv', col_types = cols(
    crystal_structure = col_character(),
    structure_file_path = col_character(),
    symbol_cation = col_character(),
    symbol_anion = col_character()
))

selected_structures <- element_pairs |>
    inner_join(structure_files,
               by = c('symbol_cation', 'symbol_anion'),
               # Many categories, many files
               relationship = 'many-to-many') |>
    mutate(structure_id = glue('{crystal_structure}_{symbol_cation}_{symbol_anion}'))

# Check uniqueness of the structure id's
with(selected_structures,
     stopifnot(length(unique(structure_id)) == length(structure_id)))

write_csv(selected_structures, 'selected_structure_files.csv')
