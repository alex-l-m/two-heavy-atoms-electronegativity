'''Extract Becke population from Q-Chem log files'''
import os.path
from glob import glob
import re
import csv
import gzip

def idle(line, output):
    '''State that represents not being in the relevant section'''
    if re.search('CDFT Becke Populations', line) is not None:
        return find_header
    else:
        return idle

def find_header(line, output):
    '''State that finds the header of the Becke population section'''
    if re.search('Atom', line) is not None:
        return parse_data
    else:
        return find_header

def parse_data(line, output):
    '''State that parses the Becke population data'''
    if re.search('---', line) is not None:
        return idle
    else:
        output.append(line.split())
        return parse_data

full_output = []
for log_file_path in glob('qchem_logs/*.log'):
    combination_id, ext = os.path.splitext(os.path.basename(log_file_path))
    with open(log_file_path) as log_file:
        state = idle
        this_mol_output = []
        for line in log_file:
            state = state(line, this_mol_output)
        # Add the molecule id to each row and add to full output
        full_output += [[combination_id] + row for row in this_mol_output]

outfile = 'raw_becke_populations.csv.gz'
header = ['combination_id', 'atom_number', 'symbol', 'excess_electrons', 'population',
          'net_spin']
with gzip.open(outfile, 'wt') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    writer.writerows(full_output)
