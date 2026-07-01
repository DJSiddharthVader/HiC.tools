######################################################################
# Depdendencies
######################################################################
library(here)
# here::i_am('scripts/loops/calculate.loop.nesting.levels.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    library(plyranges)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.loading.results.R'))
    source(file.path(SCRIPT_DIR, 'coverage/utils.coverage.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    library(magrittr)
    library(tidyverse)
})
# cli args
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)

######################################################################
# Load loops data 
######################################################################
# load all per condition loop data with/without separation by differential status across comparisons
all.loop.data.df <- 
    load_all_loop_data_for_nesting_analysis()
# Load all genomic bins for counting nesting
all.bins.df <- 
    list_all_genome_bin_files() %>%
    mutate(
        bins.df=
            pmap(
                .l=list(genomic.bins.filepath),
                .f=
                    function(genomic.bins.filepath) {
                        genomic.bins.filepath %>% 
                        read_tsv(
                            show_col_types=FALSE,
                            progress=FALSE
                        ) %>%
                        dplyr::rename('chr'=chrom)
                    }
            )
    ) %>%
    select(-c(genomic.bins.filepath))

######################################################################
# Generate bedtools cmds to calculate loop nesting
######################################################################
# Use plyranges to map loops to bins and summarize contiguous bins sets into segments
# bin.10      111111111122222222233333333333444444444455555555556666666666
# bin 01      123456789012345678901234567890123456789012345678901234567890
# Loop 1      --------------|==================|--------------------------
# Loop 2      --------------|============|--------------------------------
# Loop 3      --------------|=========|-----------------------------------
# Loop 3      --------|=====|---------------------------------------------
# nesting lvl 000000001111114333333333322211111100000000000000000000000000
# map all overlapping loops to each bin they overlap
# summarize loop statistics for each contiguously set of bins at the same nesting lvl
# compute_all_loop_nesting_results(
#     all.loop.data.df=all.loop.data.df,
#     all.bins.df=all.bins.df,
#     force_redo=parsed.args$force.redo
#     # force_redo=TRUE
# )
# # combine all nesting results into a single 
# check_cached_results(
#     results_file=ALL_LOOP_NESTING_RESULTS_FILE,
#     results_fnc=load_all_loop_nesting_results,
#     force_redo=parsed.args$force.redo
#     # force_redo=TRUE
# )

######################################################################
# Calculate rolling correlation of nesting data between conditions
######################################################################
# compute all between condition correlations of binwise nesting lvl data
compute_all_loop_nesting_correlation_results(
    all.loops.df=all.loop.data.df,
    all.bins.df=all.bins.df,
    window.sizes=c(1000),
    force_redo=parsed.args$force.redo
    # force_redo=TRUE
)
# combine correlation results into a single file
check_cached_results(
    results_file=ALL_LOOP_NESTING_CORR_RESULTS_FILE,
    results_fnc=load_all_loop_nesting_correlation_results,
    force_redo=parsed.args$force.redo
    # force_redo=TRUE
)

