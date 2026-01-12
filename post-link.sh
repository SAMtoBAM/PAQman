#!/bin/bash

##modify craq paths
sed -i 's|/../src/|/../share/CRAQ/src/|g' "$PREFIX/bin/craq"

##add multithreading to CRAQ samtools depth step
sed -i "s/samtools depth -a/samtools depth -@ \$t -a/" "$PREFIX/src/runLR.sh"
sed -i "s/samtools depth -a/samtools depth -@ \$t -a/" "$PREFIX/src/runSR.sh"

##install ggradar (now doing by conda)
#Rscript -e "devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)"

