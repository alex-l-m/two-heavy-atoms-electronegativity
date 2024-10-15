'''Use 'bader' to compute Bader charge from cube files'''

from glob import glob
from subprocess import run
import shutil
from os import mkdir
from os.path import basename, splitext

indir = 'cp2k_cube'
infiles = glob(f'{indir}/*.cube')

outdir = 'cp2k_bader'
# Create directory if it doesn't already exist
try:
    mkdir(outdir)
except FileExistsError:
    pass

suffix = '.ACF.dat'

for infile in infiles:
    prefix, ext = splitext(basename(infile))
    outfile = f'{outdir}/{prefix}{suffix}'
    run(['bader', infile])
    shutil.move('ACF.dat', outfile)
