module load NiaEnv/2022a
module load gcc/11.3.0 openmpi/4.1.4+ucx-1.11.2
module load cp2k/2024.2
module load r
# Partial charges
# CP2K Hirshfeld:
python parse_cp2k_hirshfeld.py
# Energies
Rscript parse_energies.R

module load intel
module load gnu-parallel
# Bader:
python compute_bader.py
cd cp2k_cube
parallel -j $(nproc) < bader_commands.sh
cd ..
python parse_bader.py

# Package into a tarball
# This step isn't necessary for further analysis, I'm just including it because
# I usually run the script on a remote server for the subsequent steps locally
# and this gives me a single file to transfer
tar -czf results.tar.gz *.csv *.csv.gz
