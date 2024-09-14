# Run CP2K
Rscript make_structure_file_table.R
Rscript make_element_pairs.R
Rscript select_structure_files.R
python make_cp2k_jobs.py
sh cp2k_jobs.sh
