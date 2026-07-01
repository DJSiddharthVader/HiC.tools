#!/usr/bin/env bash

OUTPUT_DIR="$(pwd)/dependencies.files"
# list of mamba envs to export to files
mamba_envs=(distiller multiqc hicrep TADs cooltools r)
for env_name in ${mamba_envs[@]}; do
    mamba env export -n ${env_name} | grep -v '^prefix:' >| ${OUTPUT_DIR}/${env_name}.yml
done
# create tsv of all installed R packages
R -q -e "library(magrittr); library(tidyverse); installed.packages() %>% as.data.frame() %>% tibble() %>% select(-c(LibPath)) %>% write_tsv('${OUTPUT_DIR}/R.packages.tsv')"

