module load NiaEnv/2022a
module load gcc/11.3.0 openmpi/4.1.4+ucx-1.11.2
module load cp2k/2024.2
module load r
# Run CP2K
Rscript make_structure_file_table.R
Rscript make_element_pairs.R
Rscript select_structure_files.R
python adjust_lattice_constant.py
Rscript scaled_structure_files_table.R
python make_cp2k_jobs.py
