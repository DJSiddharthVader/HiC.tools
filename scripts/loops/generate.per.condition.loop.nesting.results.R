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
options(future.globals.maxSize=3.3 * 1024**3)

###################################################
# Load loops data 
###################################################
# Load all genomic bins since nesting is first computed binwise and squashed
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
# load all loops called per condition + loop features
all.per.condition.loop.data.df <- 
    load_per_condition_loop_data_for_nesting_analysis() %>% 
    mutate(anchor.right=anchor.right + resolution) %>% 
    dplyr::rename(
        'start'=anchor.left,
        'end'=anchor.right
    ) %>% 
    nest(
        loops.df=
            c(
                FeatureID,
                chr, start, end, 
                loop.feature, loop.value, 
            )
    )

###################################################
# Compute nesting data for all loops within each condition
###################################################
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
# For each set of loops (differnetial or not) compute the nesting structure and save as 
# table where each row is a genomic segment + nesting level + loop summary stats 
# each segment is genomic range overlapped by the same set of loops
compute_all_loop_nesting_results(
    all.loop.data.df=all.per.condition.loop.data.df,
    all.bins.df=all.bins.df,
    force_redo=parsed.args$force.redo
    # force_redo=TRUE
)
# combine all nesting results into a single file
check_cached_results(
    results_file=,
    results_fnc=load_all_loop_nesting_results,
    force_redo=parsed.args$force.redo
    # force_redo=TRUE
)

###################################################
# calculate segment-wiose correlation of loop nesting between pairs of conditions (loop sets)
###################################################
# get loop data for all pairs of conditions being comparied
all.loop.condition.comparisons.df <- 
    all.per.condition.loop.data.df %>% 
    inner_join(
        .,
        {.},
        suffix=c('.Numerator', '.Denominator'),
        by=join_by(resolution)
    ) %>% 
    inner_join(
       ALL_SAMPLE_MERGED_MATRIX_COMPARISONS,
       by=colnames(ALL_SAMPLE_MERGED_MATRIX_COMPARISONS)
    )
# compute correlation of nesting structure of bins between all pairs of conditions
# correlation is of nesting level per bin, for all sets of adjecent bins with >= 1 loop in either condition
# so given some pair of nesting structures, we compute correaltion across bins separately 
# for each contiguous |~~~~| segment, all bins with x are ignored
# bin.10                11111111|112222222223333333333344444444445|55555|55556666666666|77777
# bin 01                12345678|901234567890123456789012345678901|23456|78901234567890|12345
# nesting condtion del: 00000000|111111443333333332221111110000000|00000|00002222222211|00000
# nesting condtion  wt: 00000000|000000444443333332221112222222222|00000|33333322221111|00000
# segmentid:            xxxxxxxx|            segment 1            |xxxxx|   segment 2  |xxxxx
# correlation bins:     xxxxxxxx|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|xxxxx|~~~~~~~~~~~~~~|xxxxx
# so in this example we compute correlation statistics separately for the two segments (i.e. sets of bins)
compute_all_loop_nesting_correlation_results(
    all.loop.data.df=all.loop.condition.comparisons.df,
    all.bins.df=all.bins.df,
    force_redo=parsed.args$force.redo
    # force_redo=TRUE
)

