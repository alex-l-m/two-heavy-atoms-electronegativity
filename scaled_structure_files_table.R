library(readr)
library(dplyr)

selected_structures_unscaled <- read_csv('selected_structure_files_unscaled.csv',
        col_types = cols(.default = col_character())

scaled_structure_data_unannotated <- read_csv('lattice_constant_structures.csv',
        col_types = cols(structure_id = col_character(),
                         unscaled_structure_id = col_character(),
                         scale_number = col_integer(),
                         scale = col_double(),
                         structure_file_path = col_character())

annotations <- selected_structures_unscaled |>
    select(-structure_file_path) |>
    reneame(unscaled_structure_id = structure_id)

# Making the columns match before combining into a single table
scaled_structure_data <- scaled_structure_data_unannotated |>
    left_join(annotations, by = unscaled_structure_id)
equilibrium_structures <- selected_structures_unscaled |>
    mutate(scale_number = 0, scale = 1, unscaled_structure_id = structure_id)
# Combining into a single table
selected_structures <- bind_rows(unscaled = equilibrium_structures,
                                 scaled = scaled_structure_data,
                                 .id = 'scale_category')

write_csv(selected_structures, 'selected_structure_files.csv')
