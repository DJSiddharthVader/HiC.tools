################################################################################
# Depdendencies
################################################################################
library(here)
here::i_am('scripts/TADs/run.TADCompare.R')
BASE_DIR <- here()
# BASE_DIR <- '/data/talkowski/Samples/cohesin_project/HiC'
suppressPackageStartupMessages({
    library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

################################################################################
# Handle arguments/parameters
################################################################################
options(scipen=999)
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# TADCompare hypper-parameters
hyper.params.df <- 
    expand_grid(
        resolutions=parsed.args$resolutions,
        normalization=c('weight', 'NONE'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# Load TAD annotations to compare 
TADs.df <- 
    load_all_TAD_results_for_TADCompare() %>%
    filter(resolution %in% parsed.args$resolutions)

################################################################################
# Generate TADCompare results from merged matrices
################################################################################
# List of pairs of merged matrices to compare
comparisons.df <- 
    # List merged matrices
    list_mcool_files() %>%
    mutate(Sample.Group=glue('{Edit}.{Celltype}.{Genotype}')) %>% 
    filter(isMerged == 'Merged') %>% 
    select(Sample.Group, filepath) %>% 
    # add the TAD boundaries to compare between matrices
    left_join(
        TADs.df,
        relationship='many-to-many',
        by=join_by(Sample.Group)
    ) %>% 
    # Specify which comparisons to evaluate
    enumerate_pairwise_comparisons(
        # sample.group.comparisons=comparisons.list,
        sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS,
        pair_grouping_cols=c('resolution', 'chr', 'TAD.method', 'TAD.params'),
        sampleID_col='Sample.Group',
        suffixes=c('Numerator', 'Denominator'),
        delim='.',
        SampleID.fields=c('Edit', 'Celltype', 'Genotype')
    )
    # select()
# used by calls to future_pmap() in functions below
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)
# Run TADCompare on everything
comparisons.df %>% 
    run_all_TADCompare(    
        hyper.params.df=hyper.params.df,
        TADs.df=TADs.df,
        force_redo=parsed.args$force.redo
    )

