#!/bin/bash

##modify craq paths
sed -i 's|/../src/|/../share/CRAQ/src/|g' "$PREFIX/bin/craq"

##install ggradar (now doing by conda)
#Rscript -e "devtools::install_github("ricardo-bion/ggradar", dependencies = TRUE)"

