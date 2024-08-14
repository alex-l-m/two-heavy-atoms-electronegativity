'''Make input files for GAMESS single atom simulations by filling in a template
with the symbol, atomic number, charge and multiplicity'''

import pandas as pd

simulation_table = pd.read_csv('single_atom_simulations.csv')

template_text = open('template.inp', 'r').read()

for row in simulation_table.itertuples():
    symbol = row.symbol
    atomic_number = row.atomic_number
    charge = row.charge
    multiplicity = row.multiplicity

    input_text = template_text.format(symbol=symbol,
                                      atomic_number=atomic_number,
                                      charge=charge, multiplicity=multiplicity)

    with open(f'{symbol}.inp', 'w') as f:
        f.write(input_text)
