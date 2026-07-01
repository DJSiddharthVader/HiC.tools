######################################################################
# Depdendencies
######################################################################
library(here)
here::i_am('scripts/distiller/make.distiller.configs.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR, 'scripts/constants.R'))
    source(file.path(BASE_DIR, 'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

######################################################################
# Make distiller configs for each sample
######################################################################
# Template config file
DISTILLER_SAMPLE_CONFIGS_DIR <- file.path(BASE_DIR, 'sample.configs')
DISTILLER_TEMPLATE_FILE      <- file.path(DISTILLER_SAMPLE_CONFIGS_DIR, 'template.distiller.yml') 
# DISTILLER_TEMPLATE_FILE
distiller.template.txt <- 
    readLines(DISTILLER_TEMPLATE_FILE) %>%
    paste(collapse='\n')
# Sub in sample-specific info to create 1 distiller-nf config file per SampleID
make.config.file <- function(
    results.dir,
    SampleID,
    fastq.R1.filepath,
    fastq.R2.filepath,
    ...){
    glue(
        distiller.template.txt,
        results.dir=results.dir,
        SampleID=SampleID,
        fastq.R1.filepath=fastq.R1.filepath,
        fastq.R2.filepath=fastq.R2.filepath
    )
}
# Individual sample metadata
# sample.metadata.df %>%
load_sample_metadata() %>% 
    add_column(results.dir=RESULTS_DIR) %>% 
    mutate(
        distiller.sample.config.filepath=
            file.path(
                DISTILLER_SAMPLE_CONFIGS_DIR,
                glue('{SampleID}.distiller.yml')
            )
    ) %>% 
    rowwise() %>% 
    mutate(
        sample.config.text=
            make.config.file(
                results.dir=RESULTS_DIR,
                SampleID=SampleID,
                fastq.R1.filepath=fastq.R1.filepath,
                fastq.R2.filepath=fastq.R2.filepath
            )
    ) %>% 
    pwalk(
        function(sample.config.text, distiller.sample.config.filepath, ...){
            cat(
                sample.config.text,
                file=distiller.sample.config.filepath
            )
        }
    )
