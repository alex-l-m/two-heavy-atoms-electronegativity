'''Read input files for the structures I'm going to be simulating, and create
versions with modified lattice constants'''

import csv
import os
import shutil
import pandas as pd
import copy
from enutils import read_structure

outdir = 'scaled_structures'
# Create the output directory. If it already exists, erase it
try:
    os.mkdir(outdir)
except FileExistsError:
    shutil.rmtree(outdir)
    os.mkdir(outdir)

outtbl_path = 'lattice_constant_structures.csv'
# Write header row
with open(outtbl_path, 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['structure_id', 'unscaled_structure_id',
        'scale_number', 'scale', 'structure_file_path'])

intbl = pd.read_csv('selected_structure_files_unscaled.csv')
scales = [.9, .95, 1.05, 1.1]
for row in intbl.itertuples():
    unscaled_structure_id = row.structure_id
    # Load the original structure
    unscaled_crystal = read_structure(row.structure_file_path)
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
        with open(outtbl_path, 'a') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow([structure_id, unscaled_structure_id,
                scale_number, scale, outstructure_path])
