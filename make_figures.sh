# Compile tikz image that is going to be used for a figure

# File containing the tikz image
IMGBASE=da_circuit
TIKZFILE=$IMGBASE.tex

# Compile with pdflatex
pdflatex $TIKZFILE
PDFFILE=$IMGBASE.pdf

# Convert pdf to png with Poppler
PNGFILE=$IMGBASE.png
pdftoppm -png -r 300 $PDFFILE > $PNGFILE

# Assemble figures into a patchwork
Rscript figures.R
