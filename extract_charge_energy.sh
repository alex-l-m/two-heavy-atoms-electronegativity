# Partial charges
# CP2K Hirshfeld:
python parse_cp2k_hirshfeld.py
# Bader:
python compute_bader.py
cd cp2k_cube
parallel -j $(nproc) < bader_commands.sh
cd ..
python parse_bader.py

# Energies
Rscript parse_energies.R
