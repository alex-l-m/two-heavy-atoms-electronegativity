'''Convert a cube file to csv files using ASE'''
import sys
import csv
from ase.io.cube import read_cube, write_cube
from ase.units import Bohr

# Path to input file
cube_path = sys.argv[1]

# Paths to output files
atoms_path = 'atoms.csv'
unit_cell_path = 'unit_cell.csv'
density_path = 'density.csv'
coordinate_system_path = 'coordinate_system.csv'

# Dictionary containing information from the cube file
cube = read_cube(open(cube_path))

# ASE atoms object for the structure
atoms = cube['atoms']

# Write position and element symbol of each atom
with open(atoms_path, 'w') as f:
    writer = csv.writer(f)
    header = ['atom_id', 'symbol', 'x', 'y', 'z']
    writer.writerow(header)
    for i, (symbol, (x, y, z)) in enumerate(zip(atoms.get_chemical_symbols(), atoms.get_positions())):
        writer.writerow([i, symbol, x, y, z])

# Write the vectors of the unit cell
# Based on resemblance to the spacing, which I know has basis vectors as
# rows, I'm assuming this does too
unit_cell = atoms.cell.tolist()
with open(unit_cell_path, 'w') as f:
    writer = csv.writer(f)
    header = ['vector', 'x', 'y', 'z']
    writer.writerow(header)
    writer.writerow(['basis_1'] + unit_cell[0])
    writer.writerow(['basis_2'] + unit_cell[1])
    writer.writerow(['basis_3'] + unit_cell[2])

# Read the origin and basis vectors of the coordinate system
origin = cube['origin']
v1 = cube['spacing'][0]
v2 = cube['spacing'][1]
v3 = cube['spacing'][2]

# Write the information about the coordinate system
with open(coordinate_system_path, 'w') as f:
    writer = csv.writer(f)
    header = ['vector', 'x', 'y', 'z']
    writer.writerow(header)
    writer.writerow(['origin'] + list(origin))
    writer.writerow(['basis_1'] + list(v1))
    writer.writerow(['basis_2'] + list(v2))
    writer.writerow(['basis_3'] + list(v3))

# Write the density
with open(density_path, 'w') as f:
    writer = csv.writer(f)
    header = ['i', 'j', 'k', 'x', 'y', 'z', 'density']
    writer.writerow(header)
    # I'm thinking of this as a matrix where each entry is a vector of values
    # Weird way to think of a 3D tensor but whatever
    for i, row in enumerate(cube['data']):
        for j, column in enumerate(row):
            for k, value in enumerate(column):
                x, y, z = i*v1 + j*v2 + k*v3
                # The value is density
                # Information on cube files:
                # https://gaussian.com/cubegen/
                # It says that all units are atomic units
                # My guess is that means that densities are in electrons per
                # cubic Bohr.  I want them in electrons per cubic Angstrom, so
                # I have to convert
                # "Bohr" I think means Angstrons per Bohr
                # (x electrons / y Bohr) * (1 Bohr / "Bohr" Angstroms)^3
                converted_density = value / Bohr**3
                writer.writerow([i, j, k, x, y, z, converted_density])
