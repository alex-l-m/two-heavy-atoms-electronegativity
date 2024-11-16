# Make table of files containing optimized structures of binary materials

library(readr)
library(dplyr)
library(stringr)
library(glue)

# Read the environment variable containing the path to the database of
# structures
environment_variable_name <- 'ELECTRONEGATIVITYSTRUCTUREPATH'
structure_path <- Sys.getenv(environment_variable_name)

structure_glob <- read_csv('structure_glob.csv', col_types = cols(
    crystal_structure = col_character(),
    infile_glob = col_character(),
    infile_regex = col_character()
))

# Find structure files for each crystal structure with the glob, and determine
# the anion and cation elements with the regex
structure_files <- structure_glob |>
    group_by(crystal_structure, infile_glob, infile_regex) |>
    reframe({
        glob_results <- Sys.glob(glue(infile_glob))
        match_results <- str_match(glob_results, glue(infile_regex))
        colnames(match_results) <- c('structure_file_path',
                                     'symbol_cation', 'symbol_anion')
        tbl <- as_tibble(match_results)
        tbl
    }) |>
    ungroup() |>
    select(-infile_glob, -infile_regex)

write_csv(structure_files, 'structure_files.csv')
