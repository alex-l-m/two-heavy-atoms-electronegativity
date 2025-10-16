#!/bin/sh
#SBATCH --account=rrg-ovoznyy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --partition=debug
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=01:00:00
module --force purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cp2k/2023.1
module load udunits
module load r
module load xml-libxml
source ~/.bashrc
conda activate chem
CHEM="$HOME/miniconda3/envs/chem"
export PATH="$CHEM/bin:$PATH"
export OMP_DISPLAY_ENV=VERBOSE
export OMP_NUM_THREADS=6
export CP2KCOMMAND='srun -n 1 -c ${OMP_NUM_THREADS} --cpu-bind=cores cp2k_shell.psmp'
