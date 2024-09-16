'''Given a cube file and a wavefunction file, run Multiwfn to obtain density
data from the wavefunction file on the same grid as in the cube file'''

import sys
import csv
from subprocess import run
from ase.io.cube import read_cube
from ase.units import Bohr

cube_path = sys.argv[1]

wfn_path = sys.argv[2]

# Dictionary containing information from the cube file
cube = read_cube(open(cube_path))

# Read the origin and basis vectors of the coordinate system
# ASE converts these to Angstroms, but I need to convert them back to Bohr for
# Multiwfn
# Section 3.6 of the manual: "The coordinates must be given in Bohr."
origin = cube['origin'] / Bohr
v1 = cube['spacing'][0] / Bohr
v2 = cube['spacing'][1] / Bohr
v3 = cube['spacing'][2] / Bohr

# Decide the number of points to sample in each direction
# I'm assuming that the wavefunction is a atom centered on the origin. For now
# I want each unit cell that's adjacent to the origin, which is eight total. I
# can make this larger in the future if atoms have numerically important
# density outside of that
original_n = cube['data'].shape[0]
assert cube['data'].shape[1] == cube['data'].shape[2] == original_n
new_n = 2 * original_n

# Adapt the origin to fully sample
new_origin = origin - v1 * original_n - v2 * original_n - v3 * original_n

# Create a file of points
# From the manual:
# "In the file the first line should be the number of points recorded, and followed by X/Y/Z coordinates of all points."
# Also write as a csv file containing the indices of the points
with open('grid.txt', 'w') as f:
    with open('grid.csv', 'w') as g:
        csv_writer = csv.writer(g)
        header = ['i', 'j', 'k']
        csv_writer.writerow(header)
        f.write(f'{new_n ** 3}\n')
        for i in range(new_n):
            for j in range(new_n):
                for k in range(new_n):
                    point = new_origin + i * v1 + j * v2 + k * v3
                    f.write(f'{point[0]} {point[1]} {point[2]}\n')
                    csv_writer.writerow([i % original_n,
                                         j % original_n,
                                         k % original_n])

multiwfn_input = f'''5
1
100
grid.txt
density.txt
q
'''

# Run Multiwfn with the target file as an argument and Multiwfn input as stdin
multiwfn_command = 'Multiwfn_noGUI'
run([multiwfn_command, wfn_path],
    input=multiwfn_input.encode('utf-8'))
