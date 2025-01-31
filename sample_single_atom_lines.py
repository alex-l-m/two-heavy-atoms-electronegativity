'''For each of the single atom reference wavefunctions, sample the density along the line using MultiWFN'''

from os.path import splitext, basename
from shutil import move
from glob import glob
from subprocess import run

# Name of MultiWFN executable
multiwfn_command = 'Multiwfn_noGUI'

# Commands to enter in MultiWFN, which are given as stdin
multiwfn_controls_path = 'line_commands.txt'

input_wfn_dir = 'single_atoms'
output_density_dir = 'single_atoms'

input_wfn_files = glob(f'{input_wfn_dir}/*.wfn')

# Execute a MultiWFN command
# They look like this, for example:
# Multiwfn_noGUI single_atoms/Al_0.wfn < line_commands.txt
# This will output a file "line.txt"
# Then, move the output file
for input_wfn_file in input_wfn_files:
    ion_id, ext = splitext(basename(input_wfn_file))
    run([multiwfn_command, input_wfn_file], stdin=open(multiwfn_controls_path))
    move('line.txt', f'{output_density_dir}/{ion_id}_line.txt')
