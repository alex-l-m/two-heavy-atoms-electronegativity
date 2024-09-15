'''Given stdin where each line is a shell command, add a header to each command with slurm scheduler settings and output each one as a separate shell script'''
from sys import stdin

slurm_header_path = 'slurm_header.sh'
slurm_header = open(slurm_header_path).read()

for i, line in enumerate(stdin):
    script_path = 'slurm_job_{}.sh'.format(i)
    with open(script_path, 'w') as script:
        script.write(slurm_header)
        script.write(line)
