# Extraction and analysis of CP2K results
python compute_bader.py
python parse_bader.py
Rscript charge_energy_cp2k.R
Rscript smooth_energy.R
Rscript smoothed_energy_plots.R
Rscript electronegativity_regression.R
