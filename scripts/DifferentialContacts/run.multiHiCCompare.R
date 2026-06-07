###################################################
# Depdendencies
###################################################
library(here)
# here::i_am('scripts/DifferentialContacts/run.multiHiCCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    library(BiocParallel)
    library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
    library(magrittr)
    library(tidyverse)
})

###################################################
# Set up all comparisons
###################################################
# All combinations of multiHiCCompare hyper-params to test
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
message(glue('using {parsed.args$threads} core to parallelize'))
register(MulticoreParam(workers=parsed.args$threads * 2 / 4), default=TRUE)
plan(multisession,      workers=parsed.args$threads * 2 / 4)

###################################################
# Generate DAC results for each comparison
###################################################
# 2 group comparison + no covariates -> use exact test
hyper.params.df <- 
    expand_grid(
        resolution=parsed.args$resolutions,
        zero.p=c(0.8),
        A.min=c(5)
    )
# GRanges object with Centro/Telomere regions to filter
data('hg38_cyto') 
# List all separate sets of individual replicates to compare + parameters to run multiHiCComapre
comparisons.df <- 
    ALL_SAMPLE_GROUP_COMPARISONS %>% 
    list_all_sample_group_comparisons(merging='both')
# Run  multiHiCCompare on everything
comparisons.df %>% 
    run_all_multiHiCCompare(    
        hyper.params.df=hyper.params.df,
        remove.regions=hg38_cyto,
        covariates.df=NULL,
        chromosomes=CHROMOSOMES,
        force_redo=parsed.args$force.redo,
        sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC,
        group1_colname='Sample.Group.P1',
        group2_colname='Sample.Group.P2'
    )

