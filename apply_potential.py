'''Run CP2K without a field, and then with increasing field strengths'''
import argparse
import re
from glob import glob
from os.path import join
from subprocess import run
import os
import shutil
import csv
import pandas as pd
import ase.io
from ase.calculators.cp2k import CP2K
from ase.io.cube import read_cube, write_cube
from enutils import read_structure

# The directories where the files should be moved to, assumed to already exist
log_file_dir_path = 'cp2k_logs'
cube_file_dir_path = 'cp2k_cube'
pot_file_dir_path = 'cp2k_pot'
hartree_pot_dir_path = 'v_hartree_cube'

# The pseudopotential that I'm using, since Zhibo used it
pseudopotential = 'GTH-PBE'

# Read path of CP2K pseudopotentials from the environment variable
cp2k_psuedopotential_varname = 'CP2KPSEUDOPOTENTIALPATH'
cp2k_pseudopotential_dir = os.environ[cp2k_psuedopotential_varname]
infile = join(cp2k_pseudopotential_dir, 'GTH_POTENTIALS')

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
sim_tbl_path = 'simulations.csv'

# Read arguments with argparse
parser = argparse.ArgumentParser()
parser.add_argument('structure_id', help='Label for the structure')
parser.add_argument('structure_path', help='Path to CP2K input file or VASP CONTCAR containing the structure')
parser.add_argument('cation', help='Cation element symbol')
parser.add_argument('anion', help='Anion element symbol')
parser.add_argument('--potential', help='Whether to apply potential (0 for no, 1 for yes)', type=int, choices=[0,1], default=1)
parser.add_argument('--kpoints', help='Size of the k-point grid', type=int, default=1)
parser.add_argument('--maxiter', help='Maximum number of iterations to run', type=int, default=None)
parser.add_argument('--initfield', help='Initial value of the field increment', type=float, default=0.001)
parser.add_argument('--dcharge', help='Target charge increment', type=float, default=0.01)
args = parser.parse_args()
structure_id = args.structure_id
structure_path = args.structure_path
cation = args.cation
anion = args.anion
apply_potential = args.potential
kpoint_grid_size = args.kpoints
max_iter = args.maxiter
# Initial field strength increment
field_strength_increment = args.initfield
# Target charge increment
# Field strength increment will be adjusted to achieve this
target_charge_increment = args.dcharge

print(f'Simulating {structure_id}')

field_strength = 0.0

# Read the crystal structure
print(f'Reading structure from {structure_path}')
structure = read_structure(structure_path)
    
# Assert that there's only two atoms
assert len(structure) == 2

# Read the template with the CP2K settings
template_path = 'ase_template_preformat.cp2k'
template_text = open(template_path).read()

potential_section_template = '''&EXTERNAL_POTENTIAL
      READ_FROM_CUBE .TRUE.
      SCALING_FACTOR {field_strength}
    &END EXTERNAL_POTENTIAL'''

if kpoint_grid_size == 1:
    kpoint_scheme = 'GAMMA'
elif kpoint_grid_size > 1:
    kpoint_scheme = f'MONKHORST-PACK {kpoint_grid_size} {kpoint_grid_size} {kpoint_grid_size}'
else:
    raise ValueError('kpoint grid size must be at least 1')

guess_options = {'fresh_start': 'ATOMIC', 'restart': 'RESTART'}

# Name of the csv file of charges to write
# Should already exist
charges_from_integration_path = f'charges_from_integration.csv'

def simulate(structure : ase.Atoms,
        structure_id : str,
        current_field_number : int,
        field_strength : float,
        first : bool,
        donor_element : str, acceptor_element : str,
        sim_working_dir : str,
        initial_charge : float) -> float:
    '''Run DFT and iterate SCF to convergence, with a certain field strength,
    recording whatever information I'm going to need later'''

    # Whether to run from a restart file
    restart = not first

    # Read the CP2K command from an environment variable
    cp2k_command_varname = 'CP2KCOMMAND'
    cp2k_command = os.environ[cp2k_command_varname]

    simulation_id = f'{structure_id}_F{current_field_number}_Vfield'

    # File names, needed for moving them later
    # Name of the cube file that gets written
    # Changing this variable doesn't set the name of the cube file, it just has
    # to match whatever gets used
    cube_file_name = f'{simulation_id}-ELECTRON_DENSITY-2.cube'
    # Path to copy the density to
    cube_file_path = join(cube_file_dir_path, cube_file_name)
    # Path to copy the potential to
    potential_file_name = f'{simulation_id}-potential.cube'
    potential_file_path = join(pot_file_dir_path, potential_file_name)
    # Name of the cube file output from the V_HARTREE_CUBE print setting
    hartree_pot_name = f'{simulation_id}-v_hartree-2.cube'
    # Path it will be moved to
    hartree_pot_path = join(hartree_pot_dir_path, hartree_pot_name)
    # Name of the log file
    log_file_name = f'{simulation_id}.out'
    # Path to the log file, after it's moved
    log_file_path = join(log_file_dir_path, log_file_name)

    # Write a row to the simulation table
    with open(sim_tbl_path, 'a') as f:
        writer = csv.writer(f)
        writer.writerow([simulation_id, 'field',
                         structure_id, donor_element, acceptor_element,
                         current_field_number, field_strength,
                         log_file_path, cube_file_path,
                         potential_file_path, hartree_pot_path])

    # Save current working directory, which is where all the scripts are located
    # This is so I can include it in the subprocess commands
    project_dir = os.getcwd()

    # Move to the working directory of the simulation
    os.chdir(sim_working_dir)

    guess = guess_options['restart'] if restart else guess_options['fresh_start']
    if not first:
        # Construct the previous combination id in order to find the cube file
        previous_field_number = current_field_number - 1
        previous_simulation_id = f'{structure_id}_F{previous_field_number}_Vfield'
        previous_cube_name = f'{previous_simulation_id}-ELECTRON_DENSITY-2.cube'
        run(['python', join(project_dir, 'to_cube.py'), previous_cube_name,
             'potential.csv', 'pot.cube'])
        potential_section = \
                potential_section_template.format(field_strength = field_strength)
    else:
        potential_section = ''
    print(f'Running simulation for {structure_id} with field strength {field_strength}')
    print(f'using guess {guess} and potential section {potential_section}')
    text = template_text.format(potential_section = potential_section,
            guess = guess,
            restart_file = 'potrestart',
            max_outer_scf = 20,
            kpoint_scheme = kpoint_scheme)
    with CP2K(label=simulation_id,
              inp=text,
              command = cp2k_command,
              # Copying Zhibo's settings
              basis_set = 'DZVP-MOLOPT-SR-GTH', pseudo_potential = 'GTH-PBE') \
            as calc:
        # I don't use this energy variable, but the energy calculation produces
        # the log file and cube file that I use
        energy = calc.get_potential_energy(structure)

    # Make the single atom density tables to use
    if first:
        donor_ions = glob(join(project_dir,
                               f'single_atoms/{donor_element}_*.wfn'))
        donor_charges = [int(re.search(r'_(-?\d+)\.wfn$', donor_ion).group(1)) \
                         for donor_ion in donor_ions]
        for path, charge in zip(donor_ions, donor_charges):
            run(['python', join(project_dir, 'cube2multiwfn.py'),
                 cube_file_name, path])
            run(['Rscript', join(project_dir, 'overlay_density.R'),
                f'{donor_element}_{charge}_density_pbc.csv'])

        acceptor_ions = glob(join(project_dir,
                                  f'single_atoms/{acceptor_element}_*.wfn'))
        acceptor_charges = [int(re.search(r'_(-?\d+)\.wfn$', acceptor_ion).group(1)) \
                            for acceptor_ion in acceptor_ions]
        for path, charge in zip(acceptor_ions, acceptor_charges):
            run(['python', join(project_dir, 'cube2multiwfn.py'),
                 cube_file_name, path])
            run(['Rscript', join(project_dir, 'overlay_density.R'),
                 f'{acceptor_element}_{charge}_density_pbc.csv'])

    # Compute charges by integration, and generate a potential to apply
    # First, convert electron density of crystal from cube to csv
    # Writes files:
    # atoms.csv
    # unit_cell.csv
    # density.csv
    # coordinate_system.csv
    run(['python', join(project_dir, 'cube2csv.py'), cube_file_name])
    # Then, integrate to obtain charges, and generate potential as a csv
    # Writes files:
    # charges.csv
    # potential.csv
    # density_with_weights.csv
    run(['Rscript', join(project_dir, 'calculate_charge.R'),
         donor_element, acceptor_element,
         'density.csv', 'unit_cell.csv', 'coordinate_system.csv', 'atoms.csv',
         'potential.csv', 'charges.csv', 'density_with_weights.csv',
         join(project_dir, 'n_valence_electrons.csv'),
         str(initial_charge)])
    # Read the charges file
    # It has columns "symbol", "valence_electrons", "donor_or_acceptor",
    # "population", "charge", and "iteration"
    charge_iterations_tbl = pd.read_csv('charges.csv')

    # Filter the charge table for the maximum iteration
    max_iteration = charge_iterations_tbl['iteration'].max()
    charge_tbl = charge_iterations_tbl[charge_iterations_tbl['iteration'] == max_iteration]

    print(charge_tbl)
    # Set the current charge based on the charge for the electron donor
    # We can assume there's only one matching row
    current_charge = charge_tbl[charge_tbl['symbol'] == acceptor_element]['charge'].values[0]
    # Assert not nan
    assert not pd.isna(current_charge)
    print(f'Current charge: {current_charge}')

    # If no potential was applied, create a cube file of all zeroes
    # representing the potential applied
    if first:
        # First, load the cube file of the density, since the new cube file
        # will re-use the unit cell, origin and Atoms object
        template_cube = read_cube(open(cube_file_name))
        # Create cube data for the zero potential, all zeroes but in the same
        # shape as the data from the density
        zero_potential = template_cube['data'] * 0
        # Write a new cube file to 'pot.cube', which is ordinarily the name of
        # the potential file read by CP2K
        write_cube(open('pot.cube', 'w'), template_cube['atoms'],
                zero_potential, template_cube['origin'])

    # Copy the output files to the project directory
    shutil.copy(log_file_name, project_dir)
    shutil.copy(cube_file_name, project_dir)
    shutil.copy(hartree_pot_name, project_dir)
    # Copying the potential that was used as input
    # Although the charge calculation output a potential to be used in the next
    # iteration of the loop, it is actually still a csv file, and 'pot.cube'
    # has not been overwritten yet
    shutil.copy('pot.cube', os.path.join(project_dir, potential_file_name))

    # Go back to the project directory
    os.chdir(project_dir)

    # Move the output files
    shutil.move(log_file_name, log_file_path)
    shutil.move(cube_file_name, cube_file_path)
    shutil.move(hartree_pot_name, hartree_pot_path)
    shutil.move(potential_file_name, potential_file_path)

    # Write charges for all iterations
    with open(charges_from_integration_path, 'a') as f:
        writer = csv.writer(f)
        for row in charge_iterations_tbl.itertuples():
            writer.writerow([simulation_id, row.symbol, row.charge, row.iteration])

    # Do a run with a single iteration just to get an energy in the absence of
    # a field
    simulation_id = f'{structure_id}_F{current_field_number}_Vnuclei'

    # Recompute filenames with the new combination id
    cube_file_name = f'{simulation_id}-ELECTRON_DENSITY-2.cube'
    cube_file_path = join(cube_file_dir_path, cube_file_name)
    log_file_name = f'{simulation_id}.out'
    log_file_path = join(log_file_dir_path, log_file_name)

    # Write a row on the no field simulation to the simulation table
    with open(sim_tbl_path, 'a') as f:
        writer = csv.writer(f)
        writer.writerow([simulation_id, 'nuclei',
                         structure_id, donor_element, acceptor_element,
                         current_field_number, field_strength,
                         log_file_path, cube_file_path,
                         # Don't bother writing applied potential and
                         # electrostatic potential, I'm not using them
                         '', ''])
    
    # Go to the simulation directory
    os.chdir(sim_working_dir)

    print(f'Running no field simulation with combination id {simulation_id}')
    # Copy the restart file to a temporary one so we're not messing up the
    # central one
    shutil.copy('potrestart', 'energyrestart')
    nofield_text = template_text.format(potential_section = '',
                                        guess = guess_options['restart'],
                                        restart_file = 'energyrestart',
                                        max_outer_scf = 1,
                                        kpoint_scheme = kpoint_scheme)
    with CP2K(label=simulation_id,
              inp=nofield_text,
              command = cp2k_command,
              max_scf = 1,
              # Copying Zhibo's settings
              basis_set = 'DZVP-MOLOPT-SR-GTH', pseudo_potential = 'GTH-PBE') \
            as calc:
        energy = calc.get_potential_energy(structure)

    # Copy the new output files to the project directory
    shutil.copy(log_file_name, project_dir)
    shutil.copy(cube_file_name, project_dir)
    
    # Go back to the project directory
    os.chdir(project_dir)

    # Move the new output files
    shutil.move(log_file_name, log_file_path)
    shutil.move(cube_file_name, cube_file_path)

    return current_charge

# Run simulations at higher field strengths until the charge becomes positive
# Later I want to put initial simulations to decide which is donor and which is acceptor
current_charge = -1.0
# The field number, which is just some integer that indexes fields
# I'm keeping it distinct from number of iterations completed in case I want to
# use negative ones for reversed charges
current_field_number = 0
first = True
n_iterations_completed = 0
# Create a directory to do the simulation in
sim_working_dir = f'{structure_id}_simulation'
os.makedirs(sim_working_dir, exist_ok = True)
# Sometimes the "cation" by position in the periodic table is actually the
# electron acceptor. I call these "reversed" materials. If it's reversed, we
# want to actually do everything in the opposite direction.
reverse = False
while (not reverse and current_charge < 0) or (reverse and current_charge > 0):
    # Store the previous charge, so I can calculate a charge difference, which
    # will be used to update the field strength increment
    previous_charge = current_charge
    # Do the simulation
    current_charge = simulate(structure, structure_id,
            current_field_number, field_strength,
            first = first,
            donor_element = cation,
            acceptor_element = anion,
            sim_working_dir = sim_working_dir,
            initial_charge = 0 if first else current_charge)
    print(f'updated charge to {current_charge}')

    # If this is the ground state (first simulation), detect reversal
    # "current_charge" is really the charge of the anion
    # So if it's positive in the ground state that's a reversal
    if first and current_charge > 0:
        print('"Anion" is positive in the ground state, applying the field in the opposite direction')
        reverse = True
        # The field also has to be applied in the opposite direction
        field_strength_increment = -field_strength_increment
        # Also, this target charge increment should be negative instead of
        # positive, since we now have a positive charge that needs to be
        # reduced to zero
        target_charge_increment = -target_charge_increment

    # If this wasn't the first stimulation, adjust the field strength increment
    # First, calculate the derivative of charge with respect to field strength
    if not first:
        charge_increment = current_charge - previous_charge
        # The charge of the anion is negative in the ground state, and should
        # be getting more positive, unless it's reversed
        assert (not reverse and charge_increment > 0) or (reverse and charge_increment < 0)
        charge_derivative = charge_increment / field_strength_increment
        # Then, calculate the required field strength increment
        field_strength_increment = target_charge_increment / charge_derivative

    # Increment the field strength
    field_strength += field_strength_increment

    # Increment the field number
    current_field_number += 1

    n_iterations_completed += 1
    if max_iter is not None and n_iterations_completed >= max_iter:
        break

    # Exit the loop if I don't want to apply potentials to this material
    if not apply_potential:
        break

    first = False

print(f'Finished {structure_id} after {n_iterations_completed} iterations')
