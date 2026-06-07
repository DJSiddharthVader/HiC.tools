###################################################
# Depdendencies
###################################################
library(here)
# here::i_am('scripts/loops/calculate.loop.nesting.levels.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    library(plyranges)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
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

###################################################
# Load loops calls
###################################################
# List of all genomic bins at each resolution
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
# Load all loop data to quantify nesting with
all.loops.df <- 
    # load loop results
    check_cached_results(
        results_file=ALL_COOLTOOLS_LOOPS_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_cooltools_dots
    ) %>%
    # Filter and clean up loops
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %>% 
    mutate(
        end=anchor.right + resolution,
        log10.qvalue=-log10(qvalue)
    ) %>% 
    select(-c(qvalue, anchor.right)) %>% 
    dplyr::rename('start'=anchor.left) %>% 
    # anchor.right is the bin start, so change it to bin end to capture that bin in each loop
    nest(
        loops.df=
            c(
                chr, start, end,
                FeatureID,
                count, length, enrichment, log10.qvalue
            )
    )

###################################################
# Generate bedtools cmds to calculate loop nesting
###################################################
# Generate bedops commands to calculate how many loops intersect with each genomic bin e.g.
# bin.10      111111111122222222233333333333444444444455555555556666666666
# bin.01      123456789012345678901234567890123456789012345678901234567890
# Loop 1      --------------|==================|--------------------------
# Loop 2      --------------|============|--------------------------------
# Loop 3      --------------|=========|-----------------------------------
# Loop 3      --------|=========|----------------------------------------
# nesting lvl 000000001111114444433333322211111100000000000000000000000000
# map all overlapping loops to each bin they overlap
# summarize loop statistics for each contiguously set of bins at the same nesting lvl
check_cached_results(
    results_file=ALL_LOOP_NESTING_RESULTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=compute_all_loop_nesting_results,
    all.loops.df=all.loops.df,
    all.bins.df=all.bins.df
)

