#!/usr/bin/env bash

OUTPUT_DIR="$(pwd)/dependencies.files"
# dependencies for 
# mamba env export -n qc3c | grep -v '^prefix:' >| ${OUTPUT_DIR}/qc3C.yml
# dependencies for 
mamba env export -n multiqc | grep -v '^prefix:' >| ${OUTPUT_DIR}/multiqc.yml
# dependencies for 
mamba env export -n distiller | grep -v '^prefix:' >| ${OUTPUT_DIR}/distiller.yml
# dependencies for running HiCRep
mamba env export -n hicrep | grep -v '^prefix:' >| ${OUTPUT_DIR}/HiCRep.yml
# dependencies for calling TADs
mamba env export -n TADs | grep -v '^prefix:' >| ${OUTPUT_DIR}/TADs.yml
# cooltools instalation
mamba env export -n cooltools | grep -v '^prefix:' >| ${OUTPUT_DIR}/cooltools.yml
# system dependencies needed for running R
mamba env export -n r | grep -v '^prefix:' >| ${OUTPUT_DIR}/R.yml
# create tsv of all installed R packages
R -q -e "library(magrittr); library(tidyverse); installed.packages() %>% as.data.frame() %>% tibble() %>% select(-c(LibPath)) %>% write_tsv('${OUTPUT_DIR}/R.packages.tsv')"

