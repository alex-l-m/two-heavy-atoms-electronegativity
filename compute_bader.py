'''Use critic2 to compute Bader charge from cube files'''

from os import mkdir
import pandas as pd

# Table containing paths to the cube files
simulation_tbl = pd.read_csv('simulations.csv')

outdir = 'cp2k_bader'
# Create directory if it doesn't already exist
try:
    mkdir(outdir)
except FileExistsError:
    pass

# critic2 command
critic2_command = 'critic2'

# critic2 template, with a blank format field "inpath"
with open('bader_template.cri') as f:
    bader_template = f.read()

# Empty list of rows (as dictionaries) for a table recording the paths to the
# input and output files
bader_tbl_rows = []

# Make a file of commands that will be read by GNU parallel
with open('bader_commands.sh', 'w') as f:
    for row in simulation_tbl.itertuples():
        simulation_id = row.simulation_id
        cube_file_path = row.cube_file_path
        # critic2 input file, by filling in the template
        critic2_input_text = bader_template.format(inpath=cube_file_path)
        # File to write the critic2 input text to
        critic2_in_path = f'{outdir}/{simulation_id}_bader.cri'
        # File to write critic2 output to
        critic2_out_path = f'{outdir}/{simulation_id}.log'
        # Write the critic2 input file
        with open(critic2_in_path, 'w') as f2:
            f2.write(critic2_input_text)
        # Full command, that takes the input file as an argument and the output
        # file as stdout
        command = f'{critic2_command} {critic2_in_path} > {critic2_out_path}'
        # Write the command
        f.write(command + '\n')
        # Create a row recording the input and output files
        bader_tbl_rows.append({
            'simulation_id': simulation_id,
            'bader_in_path': critic2_in_path,
            'bader_out_path': critic2_out_path
        })

# Write the table of input and output files as a csv
bader_tbl_path = 'critic2_bader_files.csv'
bader_tbl = pd.DataFrame(bader_tbl_rows)
bader_tbl.to_csv(bader_tbl_path, index=False)
