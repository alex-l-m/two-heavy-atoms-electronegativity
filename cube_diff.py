'''Differentiate densities with respect to field strength using the cube files'''

import sys
import csv
import os
import pandas as pd
from ase.io.cube import read_cube, write_cube

outdir = sys.argv[1]
# Create a directory for the density derivatives, if it doesn't already exist
try:
    os.mkdir(outdir)
except FileExistsError:
    pass

# Create a csv file to record the simulation id's and path to the density
# derivative
csv_file = f'{outdir}.csv'
with open(csv_file, 'w') as f:
    csv_writer = csv.writer(f)
    header = ['simulation_id', f'{outdir}_path']
    csv_writer.writerow(header)

# Open the table of simulations, which contains the field strength in the path
# to the cube
simulations = pd.read_csv('simulations.csv')
# Filter for only the simulations with the field applied
simulations_field = simulations[simulations['potential'] == 'field']
# Group by the structure id
grouped = simulations_field.groupby('structure_id')
# Loop over structures
for structure_id, group in grouped:
    # Sort by increasing field number
    group = group.sort_values('field_number')
    # Differentiation will proceed by subtracting a "right" density from a
    # "left" density and dividing by the difference in field values
    # Left density is undefined for the first simulation so most of the logic
    # has to be skipped
    right_density = None
    right_field = None
    previous_simulation_id = None
    for index, row in group.iterrows():
        left_density = right_density
        left_field = right_field

        # Read the simulation id
        simulation_id = row['simulation_id']
        # Read the field value
        right_field = row['field_value']
        # Make sure that it is the correct type
        assert isinstance(right_field, float)

        # Read the cube file
        # They should all be there, but some of them aren't, so for now stop
        # working on the structure if one is missing
        cube_path = row[sys.argv[2]]
        # Check for a missing value
        if pd.isnull(cube_path):
            # This is actually expected for the cube of a input potential,
            # since there is no input potential for the ground state
            # I should change this later though, since really that should be a
            # cube file of all zero
            # For now, continue, but in the future this should be a break and a
            # warning
            continue
        try:
            cube = read_cube(open(cube_path))
        except FileNotFoundError:
            print(f'Cube file not found for simulation {simulation_id}')
            break
        # Extract the density
        right_density = cube['data']
        if left_density is not None:
            # Path to the output cube file
            diff_path = f'{outdir}/{previous_simulation_id}.cube'
            # Calculate the derivative
            derivative = (right_density - left_density) / \
                         (right_field - left_field)
            # Write the derivative to a cube file
            # Docs say it can take a string but it actually can't
            # Maybe I should open an issue
            write_cube(open(diff_path, 'w'), cube['atoms'], derivative,
                    cube['origin'])
            # Write a row to the table of derivatives
            with open(csv_file, 'a') as f:
                csv_writer = csv.writer(f)
                csv_writer.writerow([previous_simulation_id, diff_path])
        previous_simulation_id = simulation_id
