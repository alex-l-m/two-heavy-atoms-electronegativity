'''Read charges output by the bader program and output a csv file'''
import re
import pandas as pd
from glob import glob
from os.path import basename

inpaths = glob('cp2k_bader/*.ACF.dat')

# Example input:
#    #         X           Y           Z       CHARGE      MIN DIST   ATOMIC VOL
# --------------------------------------------------------------------------------
#    1    0.000000    0.000000    0.000000   10.537902     1.138605    26.849137
#    2    1.287379    0.910315    2.229807    7.462098     1.309196    35.866635
# --------------------------------------------------------------------------------
#    VACUUM CHARGE:               0.0000
#    VACUUM VOLUME:               0.0000
#    NUMBER OF ELECTRONS:        18.0000

column_names = ['atom_number', 'x', 'y', 'z', 'population', 'min_dist', 'atomic_volume']

# Rows as dictionaries
rows = []
for inpath in inpaths:
    simulation_id = re.match(r'(.*)-ELECTRON_DENSITY-2.ACF.dat$', basename(inpath)).group(1)
    in_table = False
    for line in open(inpath):
        if re.search(r'-', line):
            in_table = not in_table
            continue
        if in_table:
            values = line.split()
            row = dict(zip(column_names, values))
            row['simulation_id'] = simulation_id
            rows.append(row)

df = pd.DataFrame(rows)
df.to_csv('bader_charges.csv.gz', index=False)
