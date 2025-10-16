module load StdEnv/2023 gcc/12.3 openmpi/4.1.5 cp2k/2023.1
module load udunits
module load r
# Not sure if necessary but I did it before installing R packages so just in case
module load xml-libxml
source ~/.bashrc
conda activate chem
CHEM="$HOME/miniconda3/envs/chem"
export PATH="$CHEM/bin:$PATH"
# Run CP2K
Rscript make_structure_file_table.R
Rscript make_element_pairs.R
Rscript select_structure_files.R
python adjust_lattice_constant.py
Rscript scaled_structure_files_table.R
python make_cp2k_jobs.py
