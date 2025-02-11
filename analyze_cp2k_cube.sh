# Derivative cube files
python cube_diff.py density_derivatives cube_file_path
python cube_diff.py vdiff pot_file_path
python cube_diff.py dhartree hartree_pot_path
python gradient_field.py

# Making images
Rscript slice.R
Rscript slice_animations.R
Rscript plot_vector_field.R
Rscript line_between_atoms.R
