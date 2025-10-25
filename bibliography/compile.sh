# Run only if "display_annotated_bibliography.aux" does not exist
if [ ! -f display_annotated_bibliography.aux ]; then
    pdflatex display_annotated_bibliography.tex
fi
bibtex display_annotated_bibliography.aux
pdflatex display_annotated_bibliography.tex
