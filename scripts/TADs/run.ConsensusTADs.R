################################################################################
# Depdendencies
################################################################################
library(here)
here::i_am('scripts/TADs/run.ConsensusTADs.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADCompare.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    library(tidyverse)
    library(magrittr)
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
# Generate ConsensusTAD results
################################################################################
# TADCompare parameters (all defaults)
hyper.params.df <- 
    expand_grid(
        resolution=parsed.args$resolutions,
        normalization=c('balanced', 'raw'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# list all individual HiC replicates per condition
set_up_sample_groups(use_merged=FALSE) %>% 
    # Generate TADCompare results for each condition
    # {.} -> sample.groups.df
    run_all_ConsensusTADs(
        hyper.params.df=hyper.params.df,
        force_redo=parsed.args$force.redo
    )

