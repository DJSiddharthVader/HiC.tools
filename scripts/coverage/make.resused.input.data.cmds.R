################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/coverage/make.resused.input.data.cmds.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    # library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'coverage/utils.coverage.R'))
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
# All combinations of tool-specific hyper-params to call TADs with
hyper.params.df <- 
    bind_rows(
        # cooltools params
        expand_grid(
            ignore.diags=c(2),
            normalization=c('balanced'),
            track.type=c('genecov'),
            contact.type=c('cis')
        )
    ) %>%  
    cross_join(tibble(resolution=parsed.args$resolutions)) %>% 
    add_column(threads=parsed.args$threads)
# parsed.args$force.redo=TRUE
# parsed.args$force.redo=FALSE

################################################################################
# Generate cmds to call TADs with specified params
################################################################################
# First binify the genome and save bin coords to files
hyper.params.df  %>%
    distinct(resolution) %>% 
    pull(resolution) %>%
    generate_all_genome_binning_cmds(
        cmds.output.filepath=file.path(SAMPLE_QC_DIR, 'generate.all.bins.files.cmds.txt'),
        force_redo=parsed.args$force.redo
    )
# Using the bin files we can compute the bin-wise gene coverage for 
# consistently orienting the PCA vectors computed from the contact matrices
hyper.params.df %>% 
    distinct(track.type) %>% 
    pull(track.type) %>% 
    generate_all_phasing_track_computation_cmds(
        cmds.output.filepath=file.path(SAMPLE_QC_DIR, 'generate.all.phasing.tracks.cmds.txt'),
        force_redo=parsed.args$force.redo
    )
# Generate expected number of contacts for a pairs of bins a given distance away
# i.e. avg number of contacts for all bins along the same TL->BR diagonal band in each contact matrix
hyper.params.df %>% 
    generate_all_distance_expectation_calculation_cmds(
        merge_status='merged',
        cmds.output.filepath=file.path(SAMPLE_QC_DIR, 'generate.all.distance.expectation.cmds.txt'),
        force_redo=parsed.args$force.redo
    )
# Generate binwise marginal contact totals
# i.e. avg number of contacts for all bins along the same TL->BR diagonal band in each contact matrix
hyper.params.df %>% 
    generate_all_marginal_coverage_calculation_cmds(
        merge_status='merged',
        cmds.output.filepath=file.path(SAMPLE_QC_DIR, 'generate.all.marginal.contacts.cmds.txt'),
        force_redo=parsed.args$force.redo
    )

