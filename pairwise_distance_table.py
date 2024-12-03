'''Print pairwise distances of atoms in each structure'''

import pandas as pd
import numpy as np
import pymatgen.io.ase
from enutils import read_structure

intbl = pd.read_csv('selected_structure_files_unscaled.csv')

outpath = 'pairwise_distance_table.csv'
# Rows of output table, as dictionaries
outrows = []
for row in intbl.itertuples():
    unscaled_structure_id = row.structure_id
    atom_a_element = row.symbol_cation
    atom_b_element = row.symbol_anion
    formula = f'{atom_a_element}{atom_b_element}'
    crystal_structure = row.crystal_structure
    crystal_ase = read_structure(row.structure_file_path)
    # Convert from ASE to pymagen
    crystal_pymatgen = pymatgen.io.ase.AseAtomsAdaptor.get_structure(crystal_ase)
    # Get the minimum distance between any pair of atoms
    all_distances = crystal_pymatgen.distance_matrix
    distance = np.min(all_distances[all_distances > 0])
    outrows.append({'structure_id': unscaled_structure_id,
        'formula': formula, 'crystal_structure': crystal_structure,
        'distance': distance})
pd.DataFrame(outrows).to_csv(outpath, index = False)
