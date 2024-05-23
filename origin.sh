Rscript element_properties.R
python cccbdb_coords_to_csv.py
Rscript mol_coords.R
mkdir -p xyz
mkdir -p qchem_input
mkdir -p qchem_logs
mkdir -p qchem_wfn
python make_qchem_input.py

# AIMAll directory has to be made now, not later, because the Q-Chem loop moves
# wavefunction files there immediately after output
mkdir -p aimall

NPROC=$(nproc)
# Leave a few processors open for AIMAll
# Using all processors for AIMAll seems to freeze my computer
# Subtract 4 from the total
AIMALL_NPROC=$(($NPROC - 4))

# Run Q-Chem on all inputs and save logs
# Clear jobs written on previous runs
> qchem_jobs.sh
> copy_qchem_wfn_jobs.sh
for INPATH in qchem_input/*.inp
do
    BASENAME=$(basename $INPATH)
    COMBINATION_ID=${BASENAME%.inp}
    OUTPATH=~/qchem/restart/$COMBINATION_ID.dat
    LOGPATH=qchem_logs/$COMBINATION_ID.log
    # Since I'm parallelizing across molecules, trying one core per job
    echo "qchem -nt 1 $INPATH > $LOGPATH" >> qchem_jobs.sh
    # Combination IDs are not all uppercase because of two letter element
    # symbols like "Cl". However, wavefunctions will was be saved uppercase.
    COMBINATION_ID_UPPERCASE=$(echo $COMBINATION_ID | tr '[:lower:]' '[:upper:]')
    INWFN=$COMBINATION_ID_UPPERCASE.wfn
    OUTWFN=aimall/$COMBINATION_ID.wfn
    echo "python modify_wfn.py $INWFN $OUTWFN" >> copy_qchem_wfn_jobs.sh
    # Moving, not deleting, because once the modification failed and I ended up
    # with nothing, and also it would be good to have the raw ones
    echo "mv $INWFN qchem_wfn" >> copy_qchem_wfn_jobs.sh
done
# Filter jobs based on whether there's already a wfn file saved
python filter_qchem_jobs.py
parallel --jobs $NPROC < qchem_jobs_filtered.sh
sh copy_qchem_wfn_jobs.sh
# Parsing information from the Q-Chem logs
grep "Convergence criterion met" qchem_logs/*.log > qchem_converged.txt
grep "SCF failed to converge" qchem_logs/*.log > qchem_unconverged.txt
python extract_becke_population.py
Rscript parse_lam.R
# Extract "Lam" values from Q-Chem simulation, which is presumably the
# potential
grep -E "Lam *-?[0-9.]+" qchem_logs/*.log > lamvals.txt

> aimall_jobs.sh
for INPATH in aimall/*.wfn
do
     BASENAME=$(basename $INPATH)
     COMBINATION_ID=${BASENAME%.wfn}
     # Directory of the wavefunction file from within the AIMAll output folder
     # Currently I'm putting them directly into the AIMAll folder
     # This is because AIMAll outputs to the input dir, not the working dir
     RELATIVE_INPATH=$BASENAME
     echo "aimqb.ish -nogui -encomp=4 $RELATIVE_INPATH" >> aimall_jobs.sh
done
cd aimall
parallel --jobs $AIMALL_NPROC < ../aimall_jobs.sh
cd ..
# -a Argument is needed because grep incorrectly infers some files are binary
grep -a -E "Warning! *Significant cumulative integration error." aimall/*.sum | grep -E -o "[A-Za-z0-9]*_[0-9]*_-?[01]*" > bad_integration_combination_id.txt

mkdir -p aimall_tbl
for INPATH in aimall/*.sum
do
    MOL=$(basename $INPATH .sum)
    python ~/repos/qtaim-utilities/parse_sum.py $INPATH aimall_tbl/$MOL
done

Rscript simulation_status.R

Rscript charge_energy_tables.R

Rscript smooth_energy.R
Rscript softness_table.R
Rscript smoothed_energy_plots.R
# Evaluate electron density on a grid in a plane
# Create directory for density results
mkdir -p densities_plane
# Run on every wavefunction file that was used for AIMAll
for INPATH in aimall/*.wfn
do
    MOLID=$(basename $INPATH .wfn)
    Multiwfn_noGUI $INPATH < density_multiwfn_input.txt
    # Rename and move output file
    mv plane.txt densities_plane/$MOLID.txt
done
# Combine densities and calculate differences
mkdir -p density_difference_tables
Rscript combine_compress_densities_plane.R
# Make gif animations
mkdir -p density_animation
for INPATH in density_difference_tables/*.csv.gz
do
    Rscript density_images.R $INPATH
done
