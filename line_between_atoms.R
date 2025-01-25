# Sample density along a line that passes through both atoms

library(tidyverse)
library(glue)

simulations <- read_csv('simulations.csv', col_types = cols(
    simulation_id = col_character(),
    potential = col_character(),
    structure_id = col_character(),
    cation = col_character(),
    anion = col_character(),
    field_number = col_double(),
    field_value = col_double(),
    log_file_path = col_character(),
    cube_file_path = col_character(),
    pot_file_path = col_character(),
    hartree_pot_path = col_character(),
    acceptor_proatom_path = col_character(),
    donor_proatom_path = col_character(),
    promolecule_path = col_character()
))

all_cube_files <- simulations |>
    # Only include simulations with the field
    filter(potential == 'field') |>
    # Select only the columns I'm going to use later, as well as all the cube
    # file path columns
    select(simulation_id, field_number, structure_id,
           cube_file_path,
           acceptor_proatom_path, donor_proatom_path) |>
    # Pivot so there's one cube file path per row
    pivot_longer(cols = -c(simulation_id, field_number, structure_id),
                 names_to = 'cube_file_type', values_to = 'cube_file_path') |>
    # Filter for the cube files that exist
    filter(!is.na(cube_file_path)) |>
    filter(file.exists(cube_file_path)) |>
    # Make a "cube id" column combining simulation id and cube file type
    mutate(cube_id = glue('{simulation_id}_{cube_file_type}'))


# Convert the simulation table to a list of rows
inrows <- split(all_cube_files, seq(nrow(all_cube_files)))

# Empty list of tables containing slices
slice_list <- list()

# Loop over the rows of the simulation table
for (row in inrows)
{
    # Run the python script to convert the cube file into csv files
    # cube2csv.py
    system(glue('python cube2csv.py {row$cube_file_path}'))

    # Table of nuclear positions for each atom
    atoms <- read_csv('atoms.csv', col_types = cols(
        atom_id = col_double(),
        symbol = col_character(),
        x = col_double(),
        y = col_double(),
        z = col_double()
    ))
    
    # The atom not at the corner
    center_atom <- atoms |>
        filter(abs(x) > 1e-5, abs(y) > 1e-5, abs(z) > 1e-5)
    
    # Unit vector pointing along the line between the atoms
    # Since one of the atoms is at the origin, this is just the position of the
    # center atom, normalized
    u <- with(center_atom, c(x, y, z) / sqrt(x^2 + y^2 + z^2))
    
    # Table of electron density at each voxel
    density <- read_csv('density.csv', col_types = cols(
        i = col_integer(),
        j = col_integer(),
        k = col_integer(),
        x = col_double(),
        y = col_double(),
        z = col_double(),
        density = col_double()
    ))
    
    # Filter for points on the line by selecting, for each i, the point that's
    # closest to the line
    this_slice <- density |>
        group_by(i) |>
        mutate(
            squarednorm = x^2 + y^2 + z^2,
            projection_squarednorm = (x * u[1] + y * u[2] + z * u[3])^2,
            squared_distance = squarednorm - projection_squarednorm,
            min_squared_distance = min(squared_distance),
            closest = squared_distance == min_squared_distance
        ) |>
        filter(closest) |>
        select(-squared_distance, -min_squared_distance, -closest) |>
        ungroup()

    # Add to the list
    slice_list[[row$cube_id]] <- this_slice
}

# Combine the slice tables into a single table
slice_tbl <- bind_rows(slice_list, .id = 'cube_id') |>
    # Join with the original table of cube files so I have metadata
    left_join(all_cube_files, by = 'cube_id',
              relationship = 'many-to-one')

# Write the slice table to a CSV file
write_csv(slice_tbl, 'lines.csv.gz')
