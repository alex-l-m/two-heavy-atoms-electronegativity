'''Read the number of valence electrons for each element from the pseudopotentials file'''

from os.path import join
import re
import pandas as pd

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
pseudopotential = 'GTH-PBE'
ve_re = re.compile(rf'([A-Z][a-z]) +{pseudopotential}-q(\d+) +{pseudopotential}')

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
