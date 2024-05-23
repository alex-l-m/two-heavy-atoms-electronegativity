'''Remove all qchem jobs that correspond to already existing wavefunction files'''

from os.path import basename, splitext
import re
from glob import glob

wfn_files = glob('aimall/*.wfn')
completed_job_mol_ids = []
for wfn_file in wfn_files:
    mol_id, _ = splitext(basename(wfn_file))
    completed_job_mol_ids.append(mol_id)

with open('qchem_jobs_filtered.sh', 'w') as outfile:
    with open('qchem_jobs.sh') as infile:
        for line in infile:
            mol_id = re.search('qchem_input/(.*).inp', line).group(1)
            if mol_id not in completed_job_mol_ids:
                outfile.write(line)
