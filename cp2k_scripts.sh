cp selected_elements.txt single_atoms
cd single_atoms
wolframscript atomic_numbers.wls
Rscript single_atom_simulation_table.R
python make_input_files.py

# Simulate a single atom, for each selected element
> gamess_jobs.sh
> wavefunction_jobs.sh
csvcut single_atom_simulations.csv -c job_id | tail +2 > job_ids.txt
while read JOBID
do
    echo "rungms ${JOBID}.inp > ${JOBID}.log" >> gamess_jobs.sh
    GAMESS_OUT=~/gamess/restart/${JOBID}.dat
    echo "python ~/repos/qtaim-utilities/extract_wfn.py ${GAMESS_OUT} > ${JOBID}.wfn" >> wavefunction_jobs.sh
done < job_ids.txt
NPROC=$(nproc)
parallel --jobs $NPROC < gamess_jobs.sh
sh wavefunction_jobs.sh
cd ..

# Parse CP2K pseudopotentials, required for assigning charges from CP2K
python parse_pseudopotentials.py

# Element data required for generating CP2K jobs
wolframscript group_numbers.wls
Rscript element_roles.R

# Run CP2K and extract charges and energies
Rscript make_element_pairs.R
python make_cp2k_jobs.py
sh cp2k_jobs.sh
Rscript charge_energy_cp2k.R

Rscript smooth_energy.R
Rscript smoothed_energy_plots.R

wolframscript retrieve_electronegativity.wls
Rscript electronegativity_regression.R
