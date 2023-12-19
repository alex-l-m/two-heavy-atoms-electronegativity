Rscript element_properties.R
python cccbdb_coords_to_csv.py
Rscript mol_coords.R
mkdir xyz
mkdir gamess_input
mkdir gamess_logs
mkdir gamess_output
python make_gamess_input.py
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
grep UNCONVERGED gamess_logs/*.log | grep -E -o "[A-Z0-9]*_[0-9]*_-?[01]*" > unconverged_combination_id.txt

# Run AIMAll on all outputs
mkdir aimall
> aimall_jobs.sh
for INPATH in gamess_output/*.dat
do
    BASENAME=$(basename $INPATH)
    COMBINATION_ID=${BASENAME%.dat}
    WFN_PATH=aimall/$COMBINATION_ID.wfn
    python ~/repos/qtaim-utilities/extract_wfn.py $INPATH --functional B3LYP > $WFN_PATH
    echo "aimqb.ish -nogui -encomp=2 $COMBINATION_ID.wfn" >> aimall_jobs.sh
done
cd aimall
parallel --jobs $(nproc) < ../aimall_jobs.sh
cd ..
