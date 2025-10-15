'''Module collecting functions that I might use in multiple scripts'''
import json
import re
from os.path import splitext
import numpy as np
from subprocess import run
from ase import Atoms
import ase.io

def read_structure(structure_path):
    # Read VASP CONTCAR with ASE
    if re.search('CONTCAR$', structure_path):
        structure = ase.io.read(structure_path)
    # xyz also gets read with ASE
    elif splitext(structure_path)[1] == '.xyz':
        structure = ase.io.read(structure_path)
    # That should be everything, value error if I've added something else without
    # realizing it
    else:
        raise ValueError('Unrecognized file type for structure {structure_id} with path {structure_path}')
    return structure
