# Extraction and analysis of CP2K results

# Partial charges
python parse_cp2k_hirshfeld.py

# Analysis of charges
Rscript charge_energy_cp2k.R
Rscript compare_charge.R
Rscript smooth_energy.R
Rscript smoothed_energy_plots.R

# Model fitting and analysis of parameters
Rscript electronegativity_regression.R
Rscript compare_structure_params.R
