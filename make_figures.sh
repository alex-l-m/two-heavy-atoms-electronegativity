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

# For single panels shown in the main text, copy all panels to show the
# supplements
cp 3-5:zincblende:0_lam_comparison.png manuscript
cp 3-5:zincblende:0_energy_with_nofield.png manuscript
