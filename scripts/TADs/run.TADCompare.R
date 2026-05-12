################################################################################
# Depdendencies
################################################################################
library(here)
here::i_am('scripts/TADs/run.TADCompare.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.Comparing.TADs.R'))
    library(tidyverse)
    library(magrittr)
    library(purrr)
})

################################################################################
# Handle arguments/parameters
################################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# used by calls to future_pmap() in functions below
if (parsed.args$threads > 1) {
    message(glue('using {parsed.args$threads} core to parallelize'))
    plan(multisession, workers=parsed.args$threads)
} else {
    plan(sequential)
}

################################################################################
# Generate TADCompare results from merged matrices
################################################################################
# TADCompare hypper-parameters
hyper.params.df <- 
    expand_grid(
        # normalization=c('weight', 'NONE'),
        normalization=c('balanced', 'raw'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# List of pairs of merged matrices to compare
# comparisons.df <- 
    # group all replicate matrices by condition
    list_all_mcool_files(merge_status='merged') %>%
    # select(isMerged, Sample.Group, filepath) %>% 
    nest(samples.df=-c(isMerged, Sample.Group)) %>% 
    # nest all individual replicates, so 1 row per sample group e.g. NIPBL.iN.DEL
    # add the TAD boundaries to compare between matrices
    left_join(
        load_all_TAD_results_for_TADCompare(),
        relationship='many-to-many',
        by=join_by(resolution, Sample.Group)
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
# Run TADCompare on everything
comparisons.df %>% 
    run_all_TADCompare(    
        hyper.params.df=hyper.params.df,
        TADs.df=TADs.df,
        force_redo=parsed.args$force.redo
    )

