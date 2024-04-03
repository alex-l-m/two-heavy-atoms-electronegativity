import pandas as pd
from ase import Atoms
import ase.io
from create_gamess_inp import parse_inp, dict2text, atoms2data

field_numbers = list(range(61))
field_values = [i / 300 - 0.1 for i in field_numbers]

# Coordinates of each molecule
# Column names: unique_id, formula, atomnumber, symbol, atomic_number, x, y, z
coordinates = pd.read_csv('coordinates.csv.gz')

outrows = []

for formula, table in coordinates.groupby('formula'):
    this_coordinates = []
    atomic_numbers = []
    symbols = []
    for row in table.itertuples():
        this_coordinates.append([row.x, row.y, row.z])
        atomic_numbers.append(row.atomic_number)
        symbols.append(row.symbol)
    ase_atoms = Atoms(symbols, this_coordinates)
    ase.io.write(f'xyz/{formula}.xyz', ase_atoms)
    this_data = atoms2data(formula, this_coordinates, symbols, atomic_numbers)
    for field_number, field_value in zip(field_numbers, field_values):
        # Check if there's a nonzero electric field
        has_field = abs(field_value) > 1e-6
        # Simulate alternate charges for neutral molecules only
        if has_field:
            charges = [0]
            multiplicities = [1]
        else:
            charges = [-1, 0, 1]
            multiplicities = [2, 1, 2]
        for charge, multiplicity in zip(charges, multiplicities):
            combination_id = f'{formula}_{field_number:02}_{charge}'
            outpath = f'gamess_input/{combination_id}.inp'
            # Template for GAMESS input
            template = parse_inp('template.inp')
            if has_field:
                template['EFIELD'] = {'EVEC(3)': str(field_value)}
            template['CONTRL']['ICHARG'] = charge
            template['CONTRL']['MULT'] = multiplicity
            template_text = dict2text(template)
            outtext = template_text + '\n' + this_data
            with open(outpath, 'w') as outfile:
                outfile.write(outtext)
            outrows.append({
                'combination_id': combination_id,
                'formula': formula,
                'field_number': field_number,
                'field_value': field_value,
                'molecule_charge': charge,
                'gamess_input_file': outpath
            })

pd.DataFrame(outrows).to_csv('gamess_simulations.csv.gz', index=False)
