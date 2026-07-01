###################################################
# Depdendencies
###################################################
library(here)
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
# options(future.globals.maxSize=2.5 * 1024**3)
options(future.globals.maxSize=4.5 * 1024**3)
# options(future.globals.maxSize=7.5 * 1024**3)

###################################################
# Compute nesting data for all loops within each condition
###################################################
# Load all loops called per condition + loop features
all.per.condition.loop.data.df <- 
    # get all loops per condition and pivot metrics columns
    load_per_condition_loops() %>%
    # clean up columns
    mutate(log10.qvalue=-log10(qvalue)) %>% 
    mutate(anchor.right=anchor.right + resolution) %>% 
    dplyr::rename('start'=anchor.left, 'end'=anchor.right) %>% 
    select(
        resolution, SampleID, 
        FeatureID, chr, start, end,
        enrichment, log10.qvalue
    ) %>% 
    # pivot so I can calculate nesting + stats over loop features
    pivot_longer(
        # c(count, enrichment, log10.qvalue),
        c(enrichment, log10.qvalue),
        names_to='feature',
        values_to='value'
    ) %>% 
    nest(
        loops.df=
            c(
                FeatureID,
                chr, start, end, 
                value, 
            )
    ) %>% 
    # Load lists of genomic bins since nesting is first computed binwise and squashed
    left_join(
        list_all_genome_bin_files() %>%
        mutate(
            bins.df=
                pmap(
                    .l=list(genomic.bins.filepath),
                    .f=
                        function(genomic.bins.filepath, ...){
                            genomic.bins.filepath %>% 
                            read_tsv(show_col_types=FALSE, progress=FALSE) %>%
                            dplyr::rename('chr'=chrom)
                        }
                )
        ) %>%
        select(-c(genomic.bins.filepath)),
        by=join_by(resolution)
    )
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
all.per.condition.loop.data.df %>% 
    # set up output filepath for each set of results
    mutate(
        results_file=
            file.path(
                ALL_LOOP_NESTING_RESULTS_DIR ,
                glue('resolution_{resolution}'),
                glue('feature_{feature}'),
                glue('{SampleID}-binwise.nesting.stats.tsv')
            )
    ) %>% 
    # For each genomic bin count how many loops overlap that bin
    # and summary stats of loop feature (e.g. pvalue) for those overlappign loops
    future_pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=compute_nesting_stats_per_bin,
        force_redo=parsed.args$force.redo,
        # force_redo=TRUE,
        return_data=FALSE,
        .progress=TRUE
    )
# Now squash binwise data into segments and save all results to a single file
check_cached_results(
    results_file=ALL_LOOP_NESTING_RESULTS_FILE,
    results_fnc=squash_all_binwise_nesting_data_into_segments,
    results.files.df=load_nesting_results('binwise.nesting'),
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    return_data=FALSE
)
    
###################################################
# calculate segment-wiose correlation of loop nesting between pairs of conditions (loop sets)
###################################################
# First list all files with the binwise results
# Now get all relevant pairs of conditions with matched binwise nesting data
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
binwise.results.files.df <- 
    load_nesting_results('binwise.nesting') %>% 
    filter(feature == 'log10.qvalue')
ALL_SAMPLE_MERGED_MATRIX_COMPARISONS %>% 
    left_join(
        binwise.results.files.df,
        relationship='many-to-many',
        by=
            join_by(
                SampleID.Numerator == SampleID
            )
    ) %>% 
    left_join(
        binwise.results.files.df,
        relationship='many-to-many',
        suffix=c('.Numerator', '.Denominator'),
        by=
            join_by(
                resolution,
                feature,
                SampleID.Denominator == SampleID
            )
    ) %>% 
    # set up output filepath for each set of results
    mutate(
        results_file=
            file.path(
                ALL_LOOP_NESTING_TESTING_DIR,
                glue('resolution_{resolution}'),
                glue('feature_{feature}'),
                glue('{SampleID.Numerator}-{SampleID.Denominator}-nesting.difference.stats.tsv')
            )
    ) %>% 
    {
        if (!parsed.args$force.redo) {
            filter(., !file.exists(results_file))
        } else {
            .
        }
    } %>% 
    arrange(desc(resolution)) %>% 
        # {.} -> tmp; tmp
    # future_pmap(
    pmap(
        .l=.,
        # .f=compute_nesting_correlation_results,
        .f=check_cached_results,
        results_fnc=compute_nesting_correlation_results,
        # force_redo=parsed.args$force.redo,
        force_redo=TRUE,
        .progress=TRUE
    )
# now load all differences results files into a single clean table
check_cached_results(
    results_file=ALL_LOOP_NESTING_DIFFERENCE_RESULTS_FILE,
    results_fnc=load_all_segmentwise_nesting_difference_results,
    results.files.df=load_nesting_results('nesting.differences'),
    # force_redo=parsed.args$force.redo,
    force_redo=TRUE,
    return_data=FALSE
)

