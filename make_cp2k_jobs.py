'''Output a shell script containing commands for all the CP2K jobs I want to run, by reading a file of anion and cation elements and outputting the jobs'''
import csv
import gzip
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

# Directory for the cube files
cube_file_dir_path = 'cp2k_cube'
try:
    mkdir(cube_file_dir_path)
except FileExistsError:
    pass

# Path to simulation table to output
sim_tbl_path = 'simulations.csv.gz'
# Column names for the simulation table
sim_tbl_header = ['simulation_id', 'potential',
                  'formula', 'cation', 'anion',
                  'field_number', 'field_value',
                  'log_file_path', 'cube_file_path']
# Create (or overwrite) the simulation table and write the header
with gzip.open(sim_tbl_path, 'wt') as f:
    writer = csv.writer(f)
    writer.writerow(sim_tbl_header)

# Name of the csv file of charge potential pairs to write
charges_from_integration_path = f'charges_from_integration.csv.gz'
# Create the file and write the header
with gzip.open(charges_from_integration_path, 'wt') as f:
    writer = csv.writer(f)
    writer.writerow(['simulation_id', 'symbol', 'charge'])

target_element_paths = 'target_elements.csv'
# Example
# symbol,role
# B,cation
# ...
anion_symbols = []
cation_symbols = []
for row in pd.read_csv(target_element_paths).itertuples():
    if row.role == 'anion':
        anion_symbols.append(row.symbol)
    elif row.role == 'cation':
        cation_symbols.append(row.symbol)

with open('cp2k_jobs.sh', 'w') as f:
    for anion, cation in product(anion_symbols, cation_symbols):
        structure_path = f'CP2Kset/{cation}{anion}.CONTCAR'
        if not exists(structure_path):
            continue
        else:
            job_command = f'python apply_potential.py {cation} {anion} True\n'
            f.write(job_command)
