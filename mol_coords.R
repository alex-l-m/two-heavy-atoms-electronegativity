library(tidyverse)
library(glue)

# Function to rotate to z axis
z_rot <- function(v, target)
{
    xz_vector <- c(v[1], v[3])
    xz_vec_len <- sqrt(sum(xz_vector^2))
    if (xz_vec_len > 0.000001)
    {
        xz_direction <- xz_vector / xz_vec_len
        # Angle with z axis
        xz_angle <- acos(xz_direction[2])
        # Rotation matrix around the y axis to eliminate x component
        xz_rotmat <- matrix(c(
                cos(xz_angle), 0, -sin(xz_angle),
                0, 1, 0,
                sin(xz_angle), 0, cos(xz_angle)),
                nrow=3, ncol=3, byrow=TRUE)
    } else {
        xz_rotmat <- diag(1, 3)
    }
    
    yz_vector <- c(v[2], v[3])
    yz_vec_len <- sqrt(sum(yz_vector^2))
    if (yz_vec_len > 0.000001)
    {
        yz_direction <- yz_vector / yz_vec_len
        # Angle with z axis
        yz_angle <- acos(yz_direction[2])
        # Rotation matrix around the x axis to eliminate y component
        yz_rotmat <- matrix(c(
                1, 0, 0,
                0, cos(yz_angle), -sin(yz_angle),
                0, sin(yz_angle), cos(yz_angle)),
                nrow=3, ncol=3, byrow=TRUE)
    } else {
        yz_rotmat <- diag(1, 3)
    }
    
    result <- xz_rotmat %*% yz_rotmat %*% matrix(target, ncol = 1)
}

# Load table of number of hydrogens required on each element, to saturate it
# and guarantee that the remaining bond is single
hs_required <- read_csv('hs_required.csv', col_types = cols(
    symbol = col_character(),
    hs_required = col_integer()
))

# Load table of element properties
element_properties <- read_csv('element_properties.csv.gz', col_types = cols(
    element = col_character(),
    atomic_number = col_integer(),
    symbol = col_character(),
    electron_affinity = col_double(),
    ionization_potential = col_double(),
    mulliken_electronegativity = col_double(),
    mulliken_hardness = col_double()
))

raw_coordinates <- read_csv('srd101_ExpCartesians.csv.gz', col_types = cols(
      atomnumber = col_double(),
      charge = col_double(),
      casno = col_character(),
      x = col_double(),
      y = col_double(),
      z = col_double(),
      state = col_double(),
      config = col_double(),
      atomtype = col_integer(),
      squib = col_character()
    )) |>
    rename(atomic_number = atomtype) |>
    # Use the atomic number to join with element properties
    left_join(element_properties, by = 'atomic_number') |>
    mutate(atom_id = glue('{symbol}{atomnumber}'))

# Make chemical formulas for each molecule
formulas <- raw_coordinates |>
    group_by(casno, state, config) |>
    count(symbol, name = 'num_atoms') |>
    group_by(casno, state, config) |>
    # Calculate formula
    summarise(formula = str_c(symbol, num_atoms, collapse = ''),
              .groups = 'drop')

# Count number of heavy atoms and number of hydrogens separately
coordinates <- raw_coordinates |>
    left_join(formulas, by = c('casno', 'state', 'config')) |>
    mutate(unique_id = str_c(casno, state, config, formula, sep = ':'))

# Filter for only those that have the required number of hydrogens, and two
# heavy atoms
selected_molecules <- coordinates |>
    left_join(hs_required, by = 'symbol') |>
    group_by(unique_id) |>
    summarize(num_hydrogens = sum(symbol == 'H'),
              hs_required = sum(hs_required, na.rm = FALSE),
              num_heavy_atoms = sum(symbol != 'H'),
              .groups = 'drop') |>
    filter(num_hydrogens == hs_required & num_heavy_atoms == 2) |>
    select(unique_id) |>
    # This one seems to be mislabeled? CAS number corresponds to HS, but it's
    # in here with formula H2S2
    filter(unique_id != '13940211:1:1:H2S2')

# Write table of unique id's with formulas, so I can verify that the formulas
# are unique
coordinates |>
    select(unique_id, formula, casno, state, config) |>
    distinct() |>
    inner_join(selected_molecules, by = 'unique_id') |>
    write_csv('selected_cccbdb_id.csv.gz')

# For each selected molecule, label each atom according to whether it is the
# electron donor or the electron acceptor
donor_or_acceptor <- selected_molecules |>
    left_join(coordinates, by = 'unique_id', relationship = 'one-to-many') |>
    # Filter for heavy atoms only
    filter(symbol != 'H') |>
    # Make columns listing the atom IDs of the electron donor (lower
    # electronegativity) and electron acceptor (higher electronegativity)
    select(unique_id, atom_id, mulliken_electronegativity) |>
    group_by(unique_id) |>
    # Detecting the minimum and maximum electronegativity atoms requires
    # special care for homonuclear molecules. We need to make sure that the
    # donor and acceptor atoms are distinct, even when their electronegativity
    # is the same. The which.min and which.max functions always detect the
    # location of the first minimum or maximum. What we want is for one to
    # detect the first, and the other the last, that meets the condition.
    # Therefore, for one of them, reverse the vectors
    summarize(donor = atom_id[which.min(mulliken_electronegativity)],
              acceptor = rev(atom_id)[which.max(rev(mulliken_electronegativity))],
              .groups = 'drop') |>
    # Convert this to a variable which labels each atom as a donor or acceptor
    pivot_longer(cols = c(donor, acceptor),
                 names_to = 'donor_or_acceptor',
                 values_to = 'atom_id')

# Coordinates of the two heavy atoms for purpose of centering
# It doesn't matter which atom is which or what order they're in, for centering
heavy_coordinates_uncentered <- selected_molecules |>
    left_join(coordinates, by = 'unique_id') |>
    filter(symbol != 'H') |>
    group_by(unique_id, formula) |>
    summarize(x1 = x[1], y1 = y[1], z1 = z[1],
              x2 = x[2], y2 = y[2], z2 = z[2],
              .groups = 'drop')

centered <- heavy_coordinates_uncentered |>
    left_join(coordinates, by = c('unique_id', 'formula'),
              relationship = 'one-to-many') |>
    # Centering by subtracting the average of the two heavy atoms
    mutate(x = x - (x1 + x2) / 2,
           y = y - (y1 + y2) / 2,
           z = z - (z1 + z2) / 2) |>
    select(unique_id, formula, atom_id, symbol, atomic_number, x, y, z)

# Centered coordinates of heavy atoms. This time it matters which one is the
# donor and which one is the acceptor, since we always want to rotate the
# acceptor to a particular side
# Actually, that doesn't seem to have worked, so don't rely on the z coordinate
# being meaningful
heavy_coordinates_centered <- centered |>
    # Filter for heavy atoms only
    filter(symbol != 'H') |>
    left_join(donor_or_acceptor, by = c('unique_id', 'atom_id'),
              relationship = 'one-to-one') |>
    group_by(unique_id, formula) |>
    # Atom 1 is donor, atom 2 is acceptor
    # If there's a way to do this with a pivot, I don't see it
    summarize(
        x1 = x[donor_or_acceptor == 'donor'],
        y1 = y[donor_or_acceptor == 'donor'],
        z1 = z[donor_or_acceptor == 'donor'],
        x2 = x[donor_or_acceptor == 'acceptor'],
        y2 = y[donor_or_acceptor == 'acceptor'],
        z2 = z[donor_or_acceptor == 'acceptor'],
        .groups = 'drop')

rotated_unlabeled <- heavy_coordinates_centered |>
    # Join the table containing all coordinates, creating a table of centered
    # coordinates that also has, for each atom, the coordinates of its molecule's
    # donor heavy atom and the coordinates of its molecule's acceptor heavy atom
    left_join(centered, by = c('unique_id', 'formula'),
              relationship = 'one-to-many') |>
    # Do the rotation
    # Have to give them new names or they'll replace each other
    group_by(unique_id, formula, atom_id) |>
    mutate(x_rot = z_rot(c(x2, y2, z2), c(x, y, z))[1],
           y_rot = z_rot(c(x2, y2, z2), c(x, y, z))[2],
           z_rot = z_rot(c(x2, y2, z2), c(x, y, z))[3]) |>
    ungroup() |>
    # Keep unique id for now, for the join with donor and acceptor table
    select(unique_id, formula, atom_id, symbol, atomic_number, x_rot, y_rot, z_rot) |>
    rename(x = x_rot, y = y_rot, z = z_rot) |>
    # Label top or bottom side, so heavy atom labels can be extended to hydrogen
    mutate(side = ifelse(z > 0, 'top', 'bottom'))

# Label sides as donor or acceptor
side_labels <- rotated_unlabeled |>
    inner_join(donor_or_acceptor, by = c('unique_id', 'atom_id'),
              relationship = 'one-to-one') |>
    select(unique_id, formula, side, donor_or_acceptor)

rotated_labeled <- rotated_unlabeled |>
    # Label by donor or acceptor
    left_join(side_labels, by = c('unique_id', 'formula', 'side'),
              relationship = 'many-to-one') |>
    # Select columns, now no longer requiring unique id
    select(formula, atom_id, symbol, atomic_number, x, y, z, donor_or_acceptor)

write_csv(rotated_labeled, 'coordinates.csv.gz')
