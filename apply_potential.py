'''Run CP2K without a field, and then with increasing field strengths'''
import argparse
import sys
import re
import gzip
from os.path import join
from subprocess import run
import os
import shutil
import csv
import pandas as pd
import ase.io
from ase.calculators.cp2k import CP2K

# The directories where the files should be moved to,  assumed to already exist
log_file_dir_path = 'cp2k_logs'
cube_file_dir_path = 'cp2k_cube'

# The pseudopotential that I'm using, since Zhibo used it
pseudopotential = 'GTH-PBE'

# Write table of number of valence electrons for each element, which will be
# used to convert electron populations into charges
# Should probably read this directory from an environment variable
cp2k_dir = '/usr/share/cp2k'
infile = join(cp2k_dir, 'GTH_POTENTIALS')

# Regex for extracting number of valence electrons
# Examples for pseudopotential GTH-PBE
# Has to include lines like this:
# Cu GTH-PBE-q11 GTH-PBE
# But exclude lines like this that aren't default for GTH-PBE:
# Cu GTH-PBE-q19
# And exclude lines like this that are from a different pseudopotential:
# Cu GTH-BLYP-q11 GTH-BLYP
# The regex must have groups for both the element and the number of valence
# electrons (after q)
ve_re = re.compile(rf'([A-Z][a-z]?) +{pseudopotential}-q(\d+) +{pseudopotential}')

# Dictionaries to be converted to rows of a dataframe, containing pairs of
# element symbol and number of valence electrons
rows = []
with open(infile) as f:
    for line in f:
        m = ve_re.match(line)
        if m is not None:
            symbol, ve = m.groups()
            rows.append({'symbol': symbol, 'valence_electrons': int(ve)})

df = pd.DataFrame(rows)
outfile = 'n_valence_electrons.csv'
df.to_csv(outfile, index=False)

# Path to simulation table to append to
sim_tbl_path = 'simulations.csv.gz'

# cation = sys.argv[1]
# anion = sys.argv[2]
# apply_potential = bool(int(sys.argv[3]))
# 
# # Size of the k-point grid
# # "1" means just the gamma point, n>1 means nxnxn grid
# kpoint_grid_size = sys.argv[4]

# Read arguments with argparse
parser = argparse.ArgumentParser()
parser.add_argument('cation', help='Cation element symbol')
parser.add_argument('anion', help='Anion element symbol')
parser.add_argument('--potential', help='Whether to apply potential (0 for no, 1 for yes)', type=int, choices=[0,1], default=1)
parser.add_argument('--kpoints', help='Size of the k-point grid', type=int, default=1)
parser.add_argument('--maxiter', help='Maximum number of iterations to run', type=int, default=None)
args = parser.parse_args()
cation = args.cation
anion = args.anion
apply_potential = args.potential
kpoint_grid_size = args.kpoint
max_iter = args.maxiter

print(f'Simulating {cation}{anion}')

field_strength = 0.0

# Read the crystal structure
structure_path = f'CP2Kset/{cation}{anion}.CONTCAR'
print(f'Reading structure from {structure_path}')
structure = ase.io.read(structure_path)
# Assert that there's only two atoms
assert len(structure) == 2

# Read the template with the CP2K settings
template_path = 'ase_template_preformat.cp2k'
template_text = open(template_path).read()

potential_section_template = '''&EXTERNAL_POTENTIAL
      READ_FROM_CUBE .TRUE.
      SCALING_FACTOR {field_strength}
    &END EXTERNAL_POTENTIAL'''

guess_options = {'fresh_start': 'ATOMIC', 'restart': 'RESTART'}

# Name of the csv file of charges to write
# Should already exist
charges_from_integration_path = f'charges_from_integration.csv.gz'

def simulate(structure : ase.Atoms,
        current_field_number : int,
        field_strength : float,
        restart : bool,
        donor_element : str, acceptor_element : str) -> float:
    '''Run DFT and iterate SCF to convergence, with a certain field strength,
    recording whatever information I'm going to need later'''

    formula = f'{donor_element}{acceptor_element}'

    simulation_id = f'{formula}_F{current_field_number}_Vfield'

    # File names, needed for moving them later
    # Name of the cube file that gets written
    # Changing this variable doesn't set the name of the cube file, it just has
    # to match whatever gets used
    cube_file_name = f'{simulation_id}-ELECTRON_DENSITY-2.cube'
    # Path to the cube file, after it's moved
    cube_file_path = join(cube_file_dir_path, cube_file_name)
    # Name of the log file
    log_file_name = f'{simulation_id}.out'
    # Path to the log file, after it's moved
    log_file_path = join(log_file_dir_path, log_file_name)


    # Write a row to the simulation table
    with gzip.open(sim_tbl_path, 'at') as f:
        writer = csv.writer(f)
        writer.writerow([simulation_id, 'field',
                         formula, donor_element, acceptor_element,
                         current_field_number, field_strength,
                         log_file_path, cube_file_path])

    guess = guess_options['restart'] if restart else guess_options['fresh_start']
    if current_field_number > 0:
        # Construct the previous combination id in order to find the cube file
        previous_field_number = current_field_number - 1
        previous_simulation_id = f'{formula}_F{previous_field_number}_Vfield'
        previous_cube_name = f'{previous_simulation_id}-ELECTRON_DENSITY-2.cube'
        previous_cube_path = join(cube_file_dir_path, previous_cube_name)
        run(['python', 'to_cube.py', previous_cube_path, 'potential.csv', 'pot.cube'])
        potential_section = \
                potential_section_template.format(field_strength = field_strength)
    else:
        potential_section = ''
    print(f'Running simulation for {formula} with field strength {field_strength}')
    print(f'using guess {guess} and potential section {potential_section}')
    text = template_text.format(potential_section = potential_section,
            guess = guess,
            restart_file = 'potrestart',
            max_outer_scf = 20)
    calc = CP2K(label=simulation_id,
            inp=text,
            # Copying Zhibo's settings
            basis_set = 'DZVP-MOLOPT-SR-GTH', pseudo_potential = 'GTH-PBE')
    structure.calc = calc
    # I don't use this energy variable, but the energy calculation produces the
    # log file and cube file that I use
    energy = structure.get_potential_energy()


    # Make the single atom density tables to use
    for charge in range(-2,3):
        run(['python', 'cube2multiwfn.py', cube_file_name,
            f'single_atoms/{donor_element}_{charge}.wfn'])
        run(['Rscript', 'overlay_density.R', f'{donor_element}_{charge}_density_pbc.csv'])
        run(['python', 'cube2multiwfn.py', cube_file_name,
            f'single_atoms/{acceptor_element}_{charge}.wfn'])
        run(['Rscript', 'overlay_density.R', f'{acceptor_element}_{charge}_density_pbc.csv'])

    # Compute charges by integration, and generate a potential to apply
    # First, convert electron density of crystal from cube to csv
    # Writes files:
    # atoms.csv
    # unit_cell.csv
    # density.csv
    # coordinate_system.csv
    run(['python', 'cube2csv.py', cube_file_name])
    # Then, integrate to obtain charges, and generate potential as a csv
    # Writes files:
    # charges.csv
    # potential.csv
    # density_with_weights.csv
    run(['Rscript', 'calculate_charge.R',
         donor_element, acceptor_element,
         'density.csv', 'unit_cell.csv', 'coordinate_system.csv', 'atoms.csv',
         'potential.csv', 'charges.csv', 'density_with_weights.csv'])
    # Read the charges file
    # It has columns "symbol", "population", and "charge"
    charge_tbl = pd.read_csv('charges.csv')
    print(charge_tbl)
    # Set the current charge based on the charge for the electron donor, in
    # this case boron
    # We can assume there's only one matching row
    current_charge = charge_tbl[charge_tbl['symbol'] == acceptor_element]['charge'].values[0]
    # Assert not nan
    assert not pd.isna(current_charge)
    print(f'Current charge: {current_charge}')

    # Write last charge to the csv file
    with gzip.open(charges_from_integration_path, 'at') as f:
        writer = csv.writer(f)
        for row in charge_tbl.itertuples():
            writer.writerow([simulation_id, row.symbol, row.charge])

    # Move the output files
    shutil.move(log_file_name, log_file_path)
    shutil.move(cube_file_name, cube_file_path)

    # Do a run with a single iteration just to get an energy in the absence of
    # a field
    simulation_id = f'{formula}_F{current_field_number}_Vnuclei'

    # Recompute filenames with the new combination id
    cube_file_name = f'{simulation_id}-ELECTRON_DENSITY-2.cube'
    cube_file_path = join(cube_file_dir_path, cube_file_name)
    log_file_name = f'{simulation_id}.out'
    log_file_path = join(log_file_dir_path, log_file_name)

    # Write a row on the no field simulation to the simulation table
    with gzip.open(sim_tbl_path, 'at') as f:
        writer = csv.writer(f)
        writer.writerow([simulation_id, 'nuclei',
                         formula, donor_element, acceptor_element,
                         current_field_number, field_strength,
                         log_file_path, cube_file_path])
    print(f'Running no field simulation with combination id {simulation_id}')
    # Copy the restart file to a temporary one so we're not messing up the
    # central one
    shutil.copy('potrestart', 'energyrestart')
    nofield_text = template_text.format(potential_section = '',
                                        guess = guess_options['restart'],
                                        restart_file = 'energyrestart',
                                        max_outer_scf = 1)
    calc = CP2K(label=simulation_id,
                inp=nofield_text, max_scf = 1)
    structure.calc = calc
    energy = structure.get_potential_energy()

    # Move the new output files
    shutil.move(log_file_name, log_file_path)
    shutil.move(cube_file_name, cube_file_path)

    return current_charge


field_strength_increment = 0.01

# Run simulations at higher field strengths until the charge becomes positive
# Later I want to put initial simulations to decide which is donor and which is acceptor
current_charge = -1.0
# The field number, which is just some integer that indexes fields
# I'm keeping it distinct from number of iterations completed in case I want to
# use negative ones for reversed charges
current_field_number = 0
first = True
n_iterations_completed = 0
while current_charge < 0:
    # Do the simulation
    current_charge = simulate(structure, current_field_number, field_strength,
            restart = not first, donor_element = cation, acceptor_element =
            anion)
    print(f'updated charge to {current_charge}')

    # Increment the field strength
    field_strength += field_strength_increment

    # Increment the field number
    current_field_number += 1

    first = False
    n_iterations_completed += 1
    if max_iter is not None and n_iterations_completed >= max_iter:
        break

    # Exit the loop if I don't want to apply potentials to this material
    if not apply_potential:
        break

print(f'Finished formula {cation}{anion} after {n_iterations_completed} iterations')
