#!/bin/bash -l
#SBATCH --account=rrg-ovoznyy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --export=ALL
#SBATCH --time=03:00:00
module --force purge
module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cp2k/2025.2
module load udunits
module load r
# Not sure if necessary but I did it before installing R packages so just in case
module load xml-libxml
source ~/.bashrc
conda activate chem
CHEM="$HOME/miniconda3/envs/chem"
export CP2KCOMMAND="srun -n 1 -c 7 --mem=25GB --exclusive --cpu-bind=cores --hint=nomultithread cp2k_shell.psmp"
export PATH="$CHEM/bin:$PATH"
export OMP_DISPLAY_ENV=VERBOSE
export OMP_NUM_THREADS=7
export SLURM_CPU_BIND_VERBOSE=1
sh run_cp2k.sh
while read COMMAND
do
    $COMMAND &
done < cp2k_jobs.sh
wait $(jobs -p)
sh extract_charge_energy.sh
