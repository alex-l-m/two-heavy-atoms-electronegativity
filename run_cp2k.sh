# Run CP2K
Rscript make_structure_file_table.R
Rscript make_element_pairs.R
Rscript select_structure_files.R
python adjust_lattice_constant.py
Rscript scaled_structure_files_table.R
python make_cp2k_jobs.py
