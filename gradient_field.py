'''Compute the gradient field of a cube file with critic2'''
import subprocess
import os
import pandas as pd

# Create a directory for the log files, if it doesn't already exist
logfile_path = 'critic2_logs'
try:
    os.mkdir(logfile_path)
except FileExistsError:
    pass

# Create three directories, one for each of the components of the gradient
for component in ['x', 'y', 'z']:
    try:
        os.mkdir(f'dhartree_grad_{component}')
    except FileExistsError:
        pass

# Input table containing simulation id's and paths to cube files
intable = pd.read_csv('dhartree.csv')
# Rows of the output table containing paths to the cube files
outrows = []
for row in intable.itertuples():
    incube = row.dhartree_path
    simulation_id = row.simulation_id
    xpath = f'dhartree_grad_x/{simulation_id}.cube'
    ypath = f'dhartree_grad_y/{simulation_id}.cube'
    zpath = f'dhartree_grad_z/{simulation_id}.cube'
    
    # Text of the template of the critic2 input file
    template_path = 'gradient_template.cri'
    template_text = open(template_path, 'r').read()
    
    # Text of the critic2 input file
    cri_text = template_text.format(incube=incube,
                                    xpath=xpath,
                                    ypath=ypath,
                                    zpath=zpath)
    
    # Run critic2 with the input file text as stdin, and capture stdout
    run_results = subprocess.run(['critic2'],
            input=cri_text.encode('utf-8'),
            stdout=subprocess.PIPE)
    stdout = run_results.stdout.decode('utf-8')
    
    # Write the captured stdout as a log file
    logfile = f'{logfile_path}/{simulation_id}.log'
    with open(logfile, 'w') as f:
        f.write(stdout)

    # Add a row to the output table
    outrows.append({'simulation_id': simulation_id,
                    'dhartree_grad_x': xpath,
                    'dhartree_grad_y': ypath,
                    'dhartree_grad_z': zpath})

# Write the output table
outtable = pd.DataFrame(outrows)
