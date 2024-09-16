#!/bin/sh
#SBATCH --account=rrg-ovoznyy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --exclusive
#SBATCH --export=ALL
#SBATCH --time=23:59:59
module load NiaEnv/2022a
module load gcc/11.3.0 openmpi/4.1.4+ucx-1.11.2
module load cp2k/2024.2
module load r
