cd single_atoms
wolframscript atomic_numbers.wls
Rscript single_atom_simulation_table.R
python make_input_files.py

# Simulate a single atom, for each selected element
> gamess_jobs.sh
> wavefunction_jobs.sh
INPATH=selected_elements.txt
while read ELEMENT
do
    echo "rungms ${ELEMENT}.inp > ${ELEMENT}.log" >> gamess_jobs.sh
    GAMESS_OUT=~/gamess/restart/${ELEMENT}.dat
    echo "python ~/repos/qtaim-utilities/extract_wfn.py ${GAMESS_OUT} > ${ELEMENT}.wfn" >> wavefunction_jobs.sh
done < $INPATH
NPROC=$(nproc)
parallel --jobs $NPROC < gamess_jobs.sh
sh wavefunction_jobs.sh
cd ..

# Parse CP2K pseudopotentials, required for assigning charges from CP2K
python parse_pseudopotentials.py

# Run CP2K and extract charges and energies
python make_cp2k_jobs.py
sh cp2k_jobs.sh
Rscript charge_energy_cp2k.R

Rscript smooth_energy.R
Rscript smoothed_energy_plots.R

wolframscript retrieve_electronegativity.wls
Rscript electronegativity_regression.R
