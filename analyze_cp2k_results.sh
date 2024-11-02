# Extraction and analysis of CP2K results
python cube_diff.py density_derivatives cube_file_path
python cube_diff.py vdiff pot_file_path
python cube_diff.py dhartree hartree_pot_path
python parse_cp2k_hirshfeld.py
python compute_bader.py
python parse_bader.py
Rscript charge_energy_cp2k.R
Rscript compare_charge.R
Rscript smooth_energy.R
Rscript smoothed_energy_plots.R
Rscript electronegativity_regression.R
Rscript compare_structure_params.R
