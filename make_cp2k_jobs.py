'''Output a shell script containing commands for all the CP2K jobs I want to run, by reading a file of anion and cation elements and outputting the jobs'''
import csv
import pandas as pd
from os.path import exists
from itertools import product
from os import mkdir

# Directory which will be used for log files
log_file_dir_path = 'cp2k_logs'
try:
    mkdir(log_file_dir_path)
except FileExistsError:
    pass

# Directory for the cube files of the densities
cube_file_dir_path = 'cp2k_cube'
try:
    mkdir(cube_file_dir_path)
except FileExistsError:
    pass

# Directory for the cube files of the potentials applied
pot_file_dir_path = 'cp2k_pot'
try:
    mkdir(pot_file_dir_path)
except FileExistsError:
    pass

# Path to simulation table to output
sim_tbl_path = 'simulations.csv'
# Column names for the simulation table
sim_tbl_header = ['simulation_id', 'potential',
                  'structure_id', 'cation', 'anion',
                  'field_number', 'field_value',
                  'log_file_path', 'cube_file_path', 'pot_file_path']
# Create (or overwrite) the simulation table and write the header
with open(sim_tbl_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(sim_tbl_header)

# Name of the csv file of charge potential pairs to write
charges_from_integration_path = f'charges_from_integration.csv'
# Create the file and write the header
with open(charges_from_integration_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['simulation_id', 'symbol', 'charge'])

with open('cp2k_jobs.sh', 'w') as f:
    for row in pd.read_csv('selected_structure_files.csv').itertuples():
        cation = row.symbol_cation
        anion = row.symbol_anion
        structure_path = row.structure_file_path
        assert exists(structure_path)
        structure_id = row.structure_id
        job_command = f'python apply_potential.py {structure_id} {structure_path} {cation} {anion} --kpoints 7\n'
        f.write(job_command)
