# Integrate a density to get the Hirshfeld charges
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)

# Read arguments
# Input files:
donor_element <- commandArgs(trailingOnly = TRUE)[1]
acceptor_element <- commandArgs(trailingOnly = TRUE)[2]
density_path <- commandArgs(trailingOnly = TRUE)[3]
unit_cell_path <- commandArgs(trailingOnly = TRUE)[4]
coordinate_system_path <- commandArgs(trailingOnly = TRUE)[5]
atoms_path <- commandArgs(trailingOnly = TRUE)[6]
# Output files:
potential_path <- commandArgs(trailingOnly = TRUE)[7]
charges_path <- commandArgs(trailingOnly = TRUE)[8]
density_with_weights_path <- commandArgs(trailingOnly = TRUE)[9]
charge_ref_path <- commandArgs(trailingOnly = TRUE)[10]
# Charge to start with:
initial_acceptor_charge <- as.double(commandArgs(trailingOnly = TRUE)[11])

# Show more significant figures in tables, so I can see if the charges are
# converging in the logs
options(pillar.sigfig = 7)

# A table mapping element symbols to donor or acceptor status
elements <- tibble(symbol = c(donor_element, acceptor_element),
                   donor_or_acceptor = c('donor', 'acceptor'))


charge_ref <- read_csv(charge_ref_path, col_types = cols(
    symbol = col_character(),
    valence_electrons = col_integer()
))

# Table of nuclear positions for each atom
atoms <- read_csv(atoms_path, col_types = cols(
    atom_id = col_double(),
    symbol = col_character(),
    x = col_double(),
    y = col_double(),
    z = col_double()
))

# Table of electron density at each voxel
density <- read_csv(density_path, col_types = cols(
    i = col_integer(),
    j = col_integer(),
    k = col_integer(),
    x = col_double(),
    y = col_double(),
    z = col_double(),
    density = col_double()
))

# Vectors defining the coordinate system. This should be retrieved from the
# "spacing" key of ASE's dictionary format for cube files. So it's not a
# lattice vector. These vectors don't point between the cells in the crystal,
# they point between voxels
coordinate_system <- read_csv(coordinate_system_path, col_types = cols(
    vector = col_character(),
    x = col_double(),
    y = col_double(),
    z = col_double()
))

# Calculate the volume of a voxel
v1 <- with(filter(coordinate_system, vector == 'basis_1'), c(x, y, z))
v2 <- with(filter(coordinate_system, vector == 'basis_2'), c(x, y, z))
v3 <- with(filter(coordinate_system, vector == 'basis_3'), c(x, y, z))
basis_mat <- matrix(c(v1, v2, v3), ncol = 3)
dv <- det(basis_mat)

# Print the volume, so I can check if it matches the volume computed by other
# programs such as TopoMS
print(glue('Volume of a voxel: {dv}'))

# Assign the nuclei of the atoms to voxels, by projecting onto the basis
# vectors of the voxels and rounding the coefficients
atoms_with_voxel_index <- atoms |>
    group_by(atom_id, symbol) |>
    reframe({
        projection <- solve(basis_mat, c(x, y, z))
        i <- round(projection[1])
        j <- round(projection[2])
        k <- round(projection[3])
        tibble(nucleus_i = i, nucleus_j = j, nucleus_k = k)
    })

# Number of voxels along each basis vector direction
n_voxels <- density |>
    summarize(n_i = max(i) + 1,
              n_j = max(j) + 1,
              n_k = max(k) + 1)

# Vectors defining the unit cell
unit_cell <- read_csv(unit_cell_path, col_types = cols(
    vector = col_character(),
    x = col_double(),
    y = col_double(),
    z = col_double()
))
unit_cell_basis <- unit_cell |>
    pivot_wider(names_from = vector, values_from = c(x, y, z))

# Load unnormalized integer charge weight functions from files
donor_weight_functions <- tibble(path = Sys.glob(glue('{donor_element}_*_density_pbc.csv'))) |>
    mutate(reference_charge = as.integer(str_extract(path, glue('{donor_element}_(-?\\d+)_density_pbc.csv'), group = 1))) |>
    group_by(reference_charge) |>
    reframe(read_csv(path, col_types = cols(
    i = col_double(),
    j = col_double(),
    k = col_double(),
    density = col_double()
))) |>
    rename(unnormalized_weight = density)
acceptor_weight_functions <- tibble(path = Sys.glob(glue('{acceptor_element}_*_density_pbc.csv'))) |>
    mutate(reference_charge = as.integer(str_extract(path, glue('{acceptor_element}_(-?\\d+)_density_pbc.csv'), group = 1))) |>
    group_by(reference_charge) |>
    reframe(read_csv(path, col_types = cols(
    i = col_double(),
    j = col_double(),
    k = col_double(),
    density = col_double()
))) |>
    rename(unnormalized_weight = density)

# Make an initial charge table, initialized at neutral atoms
charges <- charge_ref |>
    filter(symbol %in% c(donor_element, acceptor_element)) |>
    mutate(donor_or_acceptor = ifelse(symbol == donor_element, 'donor', 'acceptor'),
           population = valence_electrons,
           charge = ifelse(donor_or_acceptor == 'acceptor',
                           initial_acceptor_charge,
                           -1 * initial_acceptor_charge))
print('Initial charges:')
print(charges)
finished <- FALSE

iteration <- 0
# Create a table that contains charges for all iterations
all_iteration_charges <- charges |>
    mutate(iteration = iteration)

while (!finished) {
    # Increment the iteration
    iteration <- iteration + 1

    # Decide the fraction contributions for the donor and acceptor weight
    # functions of each charge
    fraction_contributions <- charges |>
        group_by(donor_or_acceptor) |>
        # Decide the fraction contribution of the weight function with each
        # reference charge based on the integer and fractional parts of the
        # charge
        reframe({
            integer_part <- floor(charge)
            fractional_part <- charge - integer_part
            tibble(
                reference_charge = c(integer_part, integer_part + 1),
                reference_contribution = c(1 - fractional_part, fractional_part)
            )
        })
    # Average the integer charge weight functions to get matched charge weight
    # functions
    donor_weight_function <- fraction_contributions |>
        filter(donor_or_acceptor == 'donor') |>
        left_join(donor_weight_functions, by = 'reference_charge') |>
        mutate(contribution_at_point =
               reference_contribution * unnormalized_weight) |>
        group_by(i, j, k) |>
        summarize(unnormalized_weight = sum(contribution_at_point),
                  .groups = 'drop')

    acceptor_weight_function <- fraction_contributions |>
        filter(donor_or_acceptor == 'acceptor') |>
        left_join(acceptor_weight_functions, by = 'reference_charge') |>
        mutate(contribution_at_point =
               reference_contribution * unnormalized_weight) |>
        group_by(i, j, k) |>
        summarize(unnormalized_weight = sum(contribution_at_point),
                  .groups = 'drop')
    
    weight_functions <-
        bind_rows(acceptor = acceptor_weight_function,
                  donor = donor_weight_function,
                  .id = 'donor_or_acceptor') |>
        # Join the element symbols
        left_join(elements, by = 'donor_or_acceptor') |>
        # Join the information required for centering
        left_join(atoms_with_voxel_index, by = 'symbol',
                  relationship = 'many-to-one') |>
        cross_join(n_voxels) |>
        rename(untranslated_i = i,
               untranslated_j = j,
               untranslated_k = k) |>
        # Translate the weight function ends to be centered at the nucleus, by
        # calculating a new index
        mutate(i = (untranslated_i + nucleus_i) %% n_i,
               j = (untranslated_j + nucleus_j) %% n_j,
               k = (untranslated_k + nucleus_k) %% n_k) |>
        # Normalize the weight, so that the weight at each voxel sums to one
        group_by(i, j, k) |>
        mutate(weight = unnormalized_weight / sum(unnormalized_weight))
    
    density_with_weights <- density |>
        left_join(weight_functions, by = c('i', 'j', 'k'))
    
    # Integrate the total density and print it
    total_density <- density |>
        summarize(total_density = sum(density * dv)) |>
        pull(total_density)
    print(glue('Total density: {total_density}'))
    
    old_charges <- charges

    # Calculate the electron populations by integrating the weights against the
    # density
    charges <- density_with_weights |>
        # The information that determines a charge really this should be a atom id,
        # but since I'm not currently considering unit cells with multiple
        # atoms the same element, I can just use the element symbol
        group_by(symbol, donor_or_acceptor) |>
        summarize(population = sum(weight * density * dv),
                  .groups = 'drop') |>
        # Compute the charges by subtracting the electron population from the
        # charge after screening from the core (that is, the number of valence
        # electrons)
        left_join(charge_ref, by = 'symbol') |>
        mutate(charge = valence_electrons - population)

    # Add the charges to the table of all iterations
    charges_with_iteration <- charges |>
        mutate(iteration = iteration)
    all_iteration_charges <- bind_rows(all_iteration_charges,
                                       charges_with_iteration)
    
    print('Charges in this iteration:')
    print(charges)

    print('Total charge:')
    print(sum(charges$charge))

    # Decide if we're finished
    charge_comparison <- charges |>
        left_join(old_charges, by = c('symbol', 'donor_or_acceptor')) |>
        mutate(charge_difference = abs(charge.x - charge.y)) |>
        summarize(max_charge_difference = max(charge_difference)) |>
        pull(max_charge_difference)
    finished <- charge_comparison < 0.0001

    print('Change in charge:')
    print(charge_comparison)
}
    
# I want to write the charges for every iteration
write_csv(all_iteration_charges, charges_path)

# Write the weights as a csv so I can convert it to a cube file and use it as a
# potential
potential <- density_with_weights |>
    select(i, j, k, x, y, z, donor_or_acceptor, weight) |>
    pivot_wider(names_from = donor_or_acceptor,
                values_from = weight) |>
    # Negative electrical potential on acceptor (positive chemical potential)
    mutate(potential = donor - acceptor) |>
    select(i, j, k, x, y, z, potential)

write_csv(potential, potential_path)

# Also write the density with weights for plotting purposes
write_csv(density_with_weights, density_with_weights_path)

# Make a table of the acceptor pro-atom, that is, the unnormalized weight of
# the acceptor atom
acceptor_proatom <- density_with_weights |>
    filter(donor_or_acceptor == 'acceptor') |>
    select(i, j, k, x, y, z, unnormalized_weight)
write_csv(acceptor_proatom, 'acceptor_proatom.csv')
# Also donor pro-atom
donor_proatom <- density_with_weights |>
    filter(donor_or_acceptor == 'donor') |>
    select(i, j, k, x, y, z, unnormalized_weight)
write_csv(donor_proatom, 'donor_proatom.csv')
# Also the sum of both
promolecule <- density_with_weights |>
    group_by(i, j, k, x, y, z) |>
    summarize(unnormalized_weight = sum(unnormalized_weight),
              .groups = 'drop') |>
    select(i, j, k, x, y, z, unnormalized_weight)
write_csv(promolecule, 'promolecule.csv')
