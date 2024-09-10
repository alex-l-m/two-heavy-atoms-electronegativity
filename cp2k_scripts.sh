# Generate reference data
# Parse CP2K pseudopotentials, required for assigning charges from CP2K
python parse_pseudopotentials.py
# Element data required for generating CP2K jobs
wolframscript group_numbers.wls
Rscript element_roles.R
# Electronegavities from standard tables, used as a sanity check on regression
# results
wolframscript retrieve_electronegativity.wls
# Atomic numbers, used to set up atomic reference simulations
wolframscript atomic_numbers.wls

# Create atomic references for Hirshfeld charge calculation
Rscript single_atom_simulation_table.R
python make_single_atom_input_files.py
# Simulate a single atom, for each selected element
> gamess_jobs.sh
> wavefunction_jobs.sh
csvcut single_atom_simulations.csv -c job_id | tail +2 > job_ids.txt
while read JOBID
do
    echo "rungms single_atoms/${JOBID}.inp > single_atoms/${JOBID}.log" >> gamess_jobs.sh
    GAMESS_OUT=~/gamess/restart/${JOBID}.dat
    # Command to convert to a wavefunction file
    echo "python ~/repos/qtaim-utilities/extract_wfn.py ${GAMESS_OUT} > single_atoms/${JOBID}.wfn" >> wavefunction_jobs.sh
    # Command to delete that wavefunction file if it's empty
    echo "if [ ! -s single_atoms/${JOBID}.wfn ]; then rm single_atoms/${JOBID}.wfn; fi" >> wavefunction_jobs.sh
done < job_ids.txt
NPROC=$(nproc)
parallel --jobs $NPROC < gamess_jobs.sh
# Wavefunction jobs are fast and don't need to be parallelized. Also now that I
# include two commands per wavefunction, they can't be parallelized, I would
# have to separate into two shell scripts
sh wavefunction_jobs.sh

# Run CP2K
Rscript make_structure_file_table.R
Rscript make_element_pairs.R
Rscript select_structure_files.R
python make_cp2k_jobs.py
sh cp2k_jobs.sh

# Extraction and analysis of CP2K results
Rscript charge_energy_cp2k.R
Rscript smooth_energy.R
Rscript smoothed_energy_plots.R
Rscript electronegativity_regression.R
