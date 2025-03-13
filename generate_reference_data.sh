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
mkdir -p single_atoms
Rscript single_atom_simulation_table.R
rm -f single_atoms/*.inp
python make_single_atom_input_files.py
# Simulate a single atom, for each selected element
rm -f single_atoms/*.log
rm -f single_atoms/*.wfn
rm -f single_atoms/*_line.txt
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
# Make plots of the reference atoms
python sample_single_atom_lines.py
Rscript plot_single_atom_lines.R
