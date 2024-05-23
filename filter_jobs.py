'''Filter the list of AIMAll jobs for only those that didn't already finish'''

from glob import glob
from os.path import basename, splitext
import re

# Make list of finished jobs
sum_file_paths = glob('aimall/*.sum')
finished_combination_id = set()
for sum_file_path in sum_file_paths:
    combination_id, ext = splitext(basename(sum_file_path))
    finished_combination_id.add(combination_id)

# Print filtered list of commands
for line in open('aimall_jobs.sh'):
    found = False
    for this_combination_id in finished_combination_id:
        if re.search(this_combination_id, line):
            found = True
            break
    if not found:
        print(line, end='')
