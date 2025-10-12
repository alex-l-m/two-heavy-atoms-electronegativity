'''Read input files for the structures I'm going to be simulating, and create
versions with modified lattice constants'''

import csv
import os
import shutil
import numpy as np
import pandas as pd
import copy
from enutils import read_structure

def cell_size(crystal):
    '''From a crystal, return some measure of the size of the cell. Right now
    I'm using distance between the atoms, but I might change this later. Also
    currently assumes one atom is at the origin, so really it's returning the
    distance of the second atom from the origin.'''
    # Positions of the two atoms
    first_atom = crystal.get_positions()[0]
    second_atom = crystal.get_positions()[1]
    # Check the assumption that the first position vector is the zero vector
    assert np.allclose(first_atom, np.zeros(3))
    # Return the length of the second position vector
    return np.linalg.norm(second_atom)

outdir = 'scaled_structures'
# Create the output directory. If it already exists, erase it
try:
    os.mkdir(outdir)
except FileExistsError:
    shutil.rmtree(outdir)
    os.mkdir(outdir)

new_structures_table_path = 'lattice_constant_structures.csv'
# Write header row
with open(new_structures_table_path, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['structure_id', 'unscaled_structure_id',
        'scale_number', 'scale', 'structure_file_path'])

# A table of the sizes of the cells
# I'm making this a separate table because I also want to include the unscaled structures, which aren't a row in my previous table
cell_size_table_path = 'cell_sizes.csv'
with open(cell_size_table_path, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['structure_id', 'cell_size'])

intbl = pd.read_csv('selected_structure_files_unscaled.csv')
scales = []
for row in intbl.itertuples():
    unscaled_structure_id = row.structure_id
    # Load the original structure
    unscaled_crystal = read_structure(row.structure_file_path)

    # Save the size of the unscaled crystal
    unscaled_crystal_size = cell_size(unscaled_crystal)
    with open(cell_size_table_path, 'a') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([unscaled_structure_id, unscaled_crystal_size])
    # 1-index the scales because scale 0 will be equilibrium structure
    for scale_number, scale in zip(range(1, len(scales)+1), scales):
        # Create a new structure id by adding the scale number
        structure_id = f'{unscaled_structure_id}_scale_{scale_number}'
        # Copy the original structure
        crystal = copy.deepcopy(unscaled_crystal)
        # Modify the lattice constant
        crystal.set_cell(crystal.cell * scale, scale_atoms = True)
        # Save the modified structure
        outstructure_path = os.path.join(outdir, f'{structure_id}.xyz')
        crystal.write(outstructure_path, format = 'extxyz')
        with open(new_structures_table_path, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([structure_id, unscaled_structure_id,
                scale_number, scale, outstructure_path])
        # Save the size
        crystal_size = cell_size(crystal)
        with open(cell_size_table_path, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([structure_id, crystal_size])
