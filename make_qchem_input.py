import pandas as pd
from ase import Atoms
import ase.io
from ase.calculators.qchem import QChem

# Q-Chem parameters
# The default parameters use lowercase parameter names
# However, parameter names are converted to uppercase when writing the input file
# Using lowercase parameter names to overwrite defaults and stay consistent
parameters = {'method': 'B3LYP',
              'basis': 'TZV',
              'unrestricted': 'TRUE',
              'scf_guess': 'GWH',
              'scf_guess_mix': 'TRUE'}

# Format template for the CDFT section of the input file, which doesn't seem
# writable with ASE's calculator interface and has to be added separately
cdft_section_template = '$cdft\n{target_becke_population}\n  1 1 {n_donor_atoms}\n -1 {n_donor_atoms_plus_one} {n_atoms}\n$end\n'

def write_qchem_input(combination_id, ase_atoms,
                      field_value = None, n_donor_atoms = None, n_atoms = None,
                      charge = None, multiplicity = None):
    # Write Q-Chem input file
    # Path to Q-Chem input file. ASE will add ".inp" as an extension
    outpath_noext = f'qchem_input/{combination_id}'
    outpath = f'{outpath_noext}.inp'
    # Write the input file through ASE, without CDFT section
    calc = QChem(label = outpath_noext)
    calc.parameters.update(parameters)
    # The "write_wfn" parameter contains the prefix of the wavefunction. This
    # will end up converted to all caps, and will have to be converted back to
    # title case later. I don't include a directory, because I'm not sure if it
    # can handle a directory, because I don't want the directory to be all
    # caps, and because I can easily just move the resulting wavefunction files
    # to a new directory
    calc.parameters['write_wfn'] = combination_id
    if field_value is not None:
        calc.parameters['cdft'] = 'TRUE'
    # Set charge
    if charge is not None:
        assert multiplicity is not None
        calc.parameters['charge'] = charge
    # I could remove the multiplicity argument and calculate it from charge,
    # but whatever
    if multiplicity is not None:
        calc.parameters['multiplicity'] = multiplicity

    calc.write_input(ase_atoms, properties = [])

    # Add CDFT section
    if field_value is not None:
        assert n_donor_atoms is not None
        cdft_section = cdft_section_template.format(
                target_becke_population = field_value,
                n_donor_atoms = n_donor_atoms,
                n_donor_atoms_plus_one = n_donor_atoms + 1,
                n_atoms = n_atoms)
        open(outpath, 'a').write(cdft_section)

field_numbers = list(range(61))
# Normalize to the range [-1, 1]
field_values = [i/30 - 1 for i in field_numbers]

# Coordinates of each molecule
# Column names: unique_id, formula, atomnumber, symbol, atomic_number, x, y, z
coordinates = pd.read_csv('coordinates.csv.gz')

outrows = []

for formula, table in coordinates.groupby('formula'):
    this_coordinates = []
    atomic_numbers = []
    symbols = []
    n_donor_atoms = 0
    n_atoms = 0
    for row in table.itertuples():
        this_coordinates.append([row.x, row.y, row.z])
        atomic_numbers.append(row.atomic_number)
        symbols.append(row.symbol)
        if row.donor_or_acceptor == 'donor':
            n_donor_atoms += 1
        n_atoms += 1
    ase_atoms = Atoms(symbols, this_coordinates)

    # Write XYZ file containing the structure
    # Currently this is redundant with the same operation in the script to make
    # GAMESS input. I'm not sure if I need both scripts, or if they should be
    # combined, or what, but it's not good to have both, because then it's not
    # clear reading the code where the file came from
    ase.io.write(f'xyz/{formula}.xyz', ase_atoms)

    # Generate input files for constrained DFT for various charges
    for field_number, field_value in zip(field_numbers, field_values):
        charge = 0
        combination_id = f'{formula}_F{field_number:02}_{charge}'

        write_qchem_input(combination_id, ase_atoms,
                field_value = field_value,
                n_donor_atoms = n_donor_atoms, n_atoms = n_atoms)

        # Add record of simulation attempt to table
        outrows.append({
            'combination_id': combination_id,
            'formula': formula,
            'cdft': True,
            'field_number': field_number,
            'field_value': field_value,
            'molecule_charge': charge
        })

    # Generate input files for simulations of the ground state, without CDFT
    charges = [-1, 0, 1]
    multiplicities = [2, 1, 2]
    for charge, multiplicity in zip(charges, multiplicities):
        combination_id = f'{formula}_NF_{charge}'
        write_qchem_input(combination_id, ase_atoms,
                charge = charge, multiplicity = multiplicity)

        # Add record of simulation attempt to table
        # Doesn't seem worth it to reduce duplication, but if I end up having
        # to modify both I'll find a way to eliminate duplication rather than
        # risk inconsistency
        outrows.append({
            'combination_id': combination_id,
            'formula': formula,
            'cdft': False,
            'field_number': None,
            'field_value': None,
            'molecule_charge': charge
        })

pd.DataFrame(outrows).to_csv('simulations.csv.gz', index=False)
