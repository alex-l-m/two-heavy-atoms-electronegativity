# Partial charges
# CP2K Hirshfeld:
python parse_cp2k_hirshfeld.py
# Energies
Rscript parse_energies.R

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
