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
mkdir qchem_output
python make_qchem_input.py

# Run GAMESS on all input files and save logs
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
parallel --jobs $(nproc) < gamess_jobs.sh
sh copy_jobs.sh
grep UNCONVERGED gamess_logs/*.log | grep -E -o "[A-Za-z0-9]*_[0-9]*_-?[01]*" > unconverged_combination_id.txt

# Run AIMAll on all outputs
mkdir aimall
> aimall_jobs.sh
for INPATH in gamess_output/*.dat
do
    BASENAME=$(basename $INPATH)
    COMBINATION_ID=${BASENAME%.dat}
    WFN_PATH=aimall/$COMBINATION_ID.wfn
    python ~/repos/qtaim-utilities/extract_wfn.py $INPATH --functional B3LYP > $WFN_PATH
    echo "aimqb.ish -nogui -encomp=4 $COMBINATION_ID.wfn" >> aimall_jobs.sh
done
cd aimall
parallel --jobs $(nproc) < ../aimall_jobs.sh
cd ..
# -a Argument is needed because grep incorrectly infers some files are binary
grep -a -E "Warning! *Significant cumulative integration error." aimall/*.sum | grep -E -o "[A-Za-z0-9]*_[0-9]*_-?[01]*" > bad_integration_combination_id.txt

mkdir aimall_tbl
for INPATH in aimall/*.sum
do
    MOL=$(basename $INPATH .sum)
    python ~/repos/qtaim-utilities/parse_sum.py $INPATH aimall_tbl/$MOL
done

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

