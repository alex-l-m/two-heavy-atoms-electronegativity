Rscript element_properties.R
python cccbdb_coords_to_csv.py
Rscript mol_coords.R
mkdir xyz
mkdir gamess_input
mkdir gamess_logs
mkdir gamess_output
python make_gamess_input.py

mkdir qchem_input
mkdir qchem_logs
python make_qchem_input.py

# Divide number of processors by two
# Using all processors for AIMAll seems to freeze my computer
NPROC=$(nproc)
HALF_NPROC=$(($NPROC / 2))

# Run GAMESS on all input files and save logs
# Clear jobs written on previous runs
> gamess_jobs.sh
> copy_jobs.sh
for INPATH in gamess_input/*.inp
do
    BASENAME=$(basename $INPATH)
    COMBINATION_ID=${BASENAME%.inp}
    OUTPATH=~/gamess/restart/$COMBINATION_ID.dat
    OUTPUT_DIR=gamess_output
    LOGPATH=gamess_logs/$COMBINATION_ID.log
    echo "rungms $INPATH > $LOGPATH" >> gamess_jobs.sh
    echo "cp $OUTPATH $OUTPUT_DIR" >> copy_jobs.sh
done
parallel --jobs $NPROC < gamess_jobs.sh
sh copy_jobs.sh
grep UNCONVERGED gamess_logs/*.log | grep -E -o "[A-Za-z0-9]*_[0-9]*_-?[01]*" > gamess_unconverged_combination_id.txt

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
    echo "rm $INWFN" >> copy_qchem_wfn_jobs.sh
done
parallel --jobs $NPROC < qchem_jobs.sh
sh copy_qchem_wfn_jobs.sh

# Run AIMAll on all outputs
mkdir gamess_aimall
> gamess_aimall_jobs.sh
for INPATH in gamess_output/*.dat
do
    BASENAME=$(basename $INPATH)
    COMBINATION_ID=${BASENAME%.dat}
    WFN_PATH=gamess_aimall/$COMBINATION_ID.wfn
    python ~/repos/qtaim-utilities/extract_wfn.py $INPATH --functional B3LYP > $WFN_PATH
    echo "aimqb.ish -nogui -encomp=4 $COMBINATION_ID.wfn" >> gamess_aimall_jobs.sh
done
cd gamess_aimall
parallel --jobs $HALF_NPROC < ../gamess_aimall_jobs.sh
cd ..
# -a Argument is needed because grep incorrectly infers some files are binary
grep -a -E "Warning! *Significant cumulative integration error." gamess_aimall/*.sum | grep -E -o "[A-Za-z0-9]*_[0-9]*_-?[01]*" > gamess_bad_integration_combination_id.txt

mkdir aimall
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
parallel --jobs $HALF_NPROC < ../aimall_jobs.sh
cd ..
# -a Argument is needed because grep incorrectly infers some files are binary
grep -a -E "Warning! *Significant cumulative integration error." aimall/*.sum | grep -E -o "[A-Za-z0-9]*_[0-9]*_-?[01]*" > bad_integration_combination_id.txt

mkdir gamess_aimall_tbl
for INPATH in gamess_aimall/*.sum
do
    MOL=$(basename $INPATH .sum)
    python ~/repos/qtaim-utilities/parse_sum.py $INPATH gamess_aimall_tbl/$MOL
done

mkdir aimall_tbl
for INPATH in aimall/*.sum
do
    MOL=$(basename $INPATH .sum)
    python ~/repos/qtaim-utilities/parse_sum.py $INPATH aimall_tbl/$MOL
done

# Extract "Lam" values from Q-Chem simulation, which is presumably the potential
grep -E "Lam *-?[0-9.]+" qchem_logs/*.log > lamvals.txt

Rscript simulation_status.R

Rscript charge_energy_tables.R

Rscript smooth_energy.R
Rscript smoothed_energy_plots.R
Rscript softness_table.R
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
# Make animation frames
mkdir -p density_images
for INPATH in density_difference_tables/*.csv.gz
do
    Rscript density_images.R $INPATH
done
# Make a gif animation from each folder of frames
mkdir -p density_animation
for FOLDER in density_images/*
do
    MOLID=$(basename $FOLDER)
    convert -delay 5 -loop 0 $FOLDER/*.png density_animation/$MOLID.gif
done

