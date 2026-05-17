################################################################################
# Depdendencies
################################################################################
library(here)
here::i_am('scripts/TADs/run.ConsensusTADs.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
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
        # normalization=c('balanced', 'raw'),
        normalization=c('balanced'),
        z_thresh=c(3),
        window_size=c(15),
        gap_thresh=c(0.2)
    )
# list all individual HiC replicates per condition
edit.specific.sample.groups <- 
    list_all_mcool_files(merge_status='individual') %>%
    nest(samples.df=-c(isMerged, Sample.Group)) %>% 
    mutate(SampleID=glue('{Sample.Group}.Consensus.Consensus'))
# Group all WT replicates across all edits together (per celltype)
cross.edit.sample.groups <- 
    list_all_mcool_files(merge_status='individual') %>%
    parse_metadata_from_names(info.format='Sample.Group') %>% 
    nest(samples.df=-c(isMerged, Celltype, Genotype)) %>% 
    filter(Genotype %in% c('WT')) %>% 
    mutate(Sample.Group=glue('All.{Celltype}.{Genotype}')) %>% 
    mutate(SampleID=glue('{Sample.Group}.Consensus.Consensus')) %>% 
    select(-c(Celltype, Genotype))
    # Generate TADCompare results for each condition
bind_rows(
    edit.specific.sample.groups,
    cross.edit.sample.groups
    ) %>% 
        # {.} -> sample.groups.df; sample.groups.df
    run_all_ConsensusTADs(
        hyper.params.df=hyper.params.df,
        force_redo=parsed.args$force.redo
    )

