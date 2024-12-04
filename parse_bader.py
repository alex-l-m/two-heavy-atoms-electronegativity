'''Read charges output by critic2 and output a csv file'''
import re
import pandas as pd

# Example input
# # Id   cp   ncp   Name  Z   mult     Volume            Pop             Lap       
#   1    1    1      Ga   31  --   5.52119041E+01  1.09868719E+01  1.27120802E-02 
#   2    2    2      N_   7   --   1.03713936E+02  7.01314410E+00 -1.27120802E-02 
# --------------------------------------------------------------------------------
#   Sum                            1.58925840E+02  1.80000160E+01 -6.98087421E-14 

column_names = ['atom_number', 'x', 'y', 'z', 'population', 'min_dist', 'atomic_volume']

# Table containing paths to the output files from critic2
critic2_tbl = pd.read_csv('critic2_bader_files.csv')

# Rows as dictionaries
rows = []
for row in critic2_tbl.itertuples():
    simulation_id = row.simulation_id
    inpath = row.bader_out_path
    # Variable tracking whether we have reached a table in the input
    in_table = False
    for line in open(inpath):
        if re.match(r'# Integrable properties', line):
            in_table = not in_table
            continue
        if in_table:
            if re.match(r'-+', line):
                # End of table
                # Setting in_table to false is redundant but whatever
                in_table = False
                # Break to avoid parsing the rest of the file, since each
                # contains only one table
                break
            elif re.match(r'#', line):
                # This is the header
                # Read everything but the initial "#" as column names. They
                # should be separated by an arbitrary number of spaces
                column_names = re.sub(r'#\s+', '', line).split()
            else:
                # The actual table
                # Read everything, again should be separated by spaces
                values = line.split()
                # Create a row as a dictionary, using column names as keys
                row = dict(zip(column_names, values))
                # Add the simulation id to the row
                row['simulation_id'] = simulation_id
                # Append the row
                rows.append(row)

df = pd.DataFrame(rows)
df.to_csv('bader_charges_raw.csv.gz', index=False)
