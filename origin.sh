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
