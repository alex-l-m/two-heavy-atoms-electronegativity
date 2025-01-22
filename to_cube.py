'''Given a cube file and a csv, create a new cube file by replacing the data'''
import sys
import pandas as pd
from ase.io.cube import read_cube, write_cube

# Path to input cube file
cube_path = sys.argv[1]

# Path to input data as a csv file
data_path = sys.argv[2]

# Path to output cube to
output_cube_path = sys.argv[3]

# Name of the column containing the data
colname = sys.argv[4]

# Dictionary containing information from the cube file
cube = read_cube(open(cube_path))

# ASE atoms object for the structure
atoms = cube['atoms']

# Origin of coordinate system
origin = cube['origin']

# Create new data by overwriting old data, which already has the right shape
new_data = cube['data']

# Dataframe with the potential
potential_df = pd.read_csv(data_path)
for index, row in potential_df.iterrows():
    i = int(row['i'])
    j = int(row['j'])
    k = int(row['k'])
    new_data[i][j][k] = row[colname]

write_cube(open(output_cube_path, 'w'), atoms, new_data, origin)
