'''Module collecting functions that I might use in multiple scripts'''
import json
import re
from os.path import splitext
import numpy as np
from subprocess import run
from ase import Atoms
import ase.io

def cp2k2ase(inpath):
    '''Function to load a CP2K input file (given by inpath) and return a ASE
    Atoms object with periodic boundary conditions'''
    infile_json = json.loads(run(['fromcp2k', inpath], capture_output = True,
                                 text = True).stdout)

    # Read unit cell
    # Dictionary with keys 'A', 'B' and 'C' containing the lattice vectors as
    # lists
    cell_dict = infile_json['force_eval']['subsys']['cell']
    # Numpy array with lattice vectors as rows
    cell = np.array([cell_dict['A'], cell_dict['B'], cell_dict['C']],
                     dtype = float)

    # Read coordinates
    # Strings containing element and xyz coordinates for each atom
    # Are they supposed to be in a key called '*' and given as strings? This
    # may be a bug in fromcp2k
    atom_rows = infile_json['force_eval']['subsys']['coord']['*']
    # Element symbol for each atom
    symbols = [re.match(r'[A-Z][a-z]?', atom_row).group() for \
               atom_row in atom_rows]
    # Coordinates, as a list of strings
    coord_strings = [re.findall(r'-?\d+\.\d+', atom_row) for \
                     atom_row in atom_rows]
    # Coordinates, as a numpy array
    # Doesn't seem to be necessary to convert to float in Python, numpy will
    # convert from string
    coords = np.array(coord_strings, dtype = float)

    # Create the Atoms object
    crystal = Atoms(symbols = symbols, scaled_positions = coords, cell = cell,
                    pbc = True)

    return crystal

def read_structure(structure_path):
    # Read VASP CONTCAR with ASE
    if re.search('CONTCAR$', structure_path):
        structure = ase.io.read(structure_path)
    # Read CP2K input file (assumed to have extension .cp2k) using the Python
    # "cp2k-input-tools"
    elif splitext(structure_path)[1] == '.cp2k':
        structure = cp2k2ase(structure_path)
    # xyz also gets read with ASE
    elif splitext(structure_path)[1] == '.xyz':
        structure = ase.io.read(structure_path)
    # That should be everything, value error if I've added something else without
    # realizing it
    else:
        raise ValueError('Unrecognized file type for structure {structure_id} with path {structure_path}')
    return structure
