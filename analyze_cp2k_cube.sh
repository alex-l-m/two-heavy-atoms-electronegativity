module load NiaEnv/2022a
module load gcc/11.3.0 openmpi/4.1.4+ucx-1.11.2
module load cp2k/2024.2
module load r

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
