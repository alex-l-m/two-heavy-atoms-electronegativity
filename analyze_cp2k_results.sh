# Extraction and analysis of CP2K results

# Derivative cube files
python cube_diff.py density_derivatives cube_file_path
python cube_diff.py vdiff pot_file_path
python cube_diff.py dhartree hartree_pot_path
python gradient_field.py

# Partial charges
python parse_cp2k_hirshfeld.py
python compute_bader.py
python parse_bader.py

# Analysis of charges
Rscript charge_energy_cp2k.R
Rscript compare_charge.R
Rscript smooth_energy.R
Rscript smoothed_energy_plots.R

# Model fitting and analysis of parameters
Rscript electronegativity_regression.R
Rscript compare_structure_params.R

# Making images
Rscript slice.R
Rscript slice_animations.R
