###################################################
# Depdendencies
###################################################
library(here)
# here::i_am('scripts/loops/run.IDR2D.loops.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    library(magrittr)
    library(tidyverse)
})
# Cli args
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
message(glue('using {parsed.args$threads} core to parallelize'))
    )
plan(multisession, workers=parsed.args$threads)

###################################################
# load loops to compare for each conditions
###################################################
# Prepare loops for comparison between conditions 
# each row is 1 nested set of loop calls per condition + context
nested.loops.df <- 
    # load loop results
    check_cached_results(
        results_file=ALL_COOLTOOLS_LOOPS_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_cooltools_dots
    ) %>%
    # Filter and clean up loops
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %>% 
    standardize_data_cols() %>% 
    # prep columns for input to IDR2D
    mutate(resolution=scale_numbers(resolution, force_numeric=TRUE)) %>% 
    dplyr::rename(
        'start.A'=anchor.left,
        'start.B'=anchor.right
    ) %>% 
    mutate(
        log10.qvalue=-log10(qvalue),
        end.A=start.A + resolution,
        end.B=start.B + resolution,
        chr.A=chr,
        chr.B=chr
    ) %>% 
    select(-c(qvalue)) %>% 
    # nest so one set of loop calls per row (SampleID + res + chr + weight)
    nest(
        loops=
            c(
                FeatureID,
                chr.A, start.A, end.A,
                chr.B, start.B, end.B, 
                count, length, enrichment, log10.qvalue
            )
    )

###################################################
# Run IDR2D across all params
###################################################
# All IDR2D hyper-params to compute
hyper.params.df <- 
    tribble(
        ~metric_colname, ~value_transformation,
        'enrichment',    'identity',  # high enrichment => most important loops 
        'log10.qvalue',  'identity'   # high log10.qval => most important loops
    ) %>% 
    cross_join(
        tibble(
            ambiguity_resolution_method=
                c(
                    "overlap",
                    "midpoint",
                    "value"
                )
        )
    ) %>% 
    cross_join(
        tibble(
            max_gap_bins=
                c(
                    0,
                    1,
                    2,
                    5
                )
        )
    )
# list all pairs of matrices to compare  loops for
comparisons.df <- 
    ALL_SAMPLE_GROUP_COMPARISONS %>% 
    rename_with(~ str_replace(.x, 'Sample.Group', 'SampleID')) %>% 
    mutate(across(everything(), ~ str_replace(.x, '$', '.Merged.Merged')))
# run IDR2D on all comparisons of sample groups + param sets
nested.loops.df %>% 
    run_all_IDR2D_analysis(
        hyper.params.df=hyper.params.df,
        force.redo=parsed.args$force.redo,
        # force.redo=TRUE,
        sample.group.comparisons=comparisons.df,
        # only compare loop call sets with matching param values for these columns
        pair_grouping_cols=
            c(
                'isMerged',
                'method',
                'kernel',
                'type',
                'normalization',
                'resolution',
                'chr'
            )
    )

###################################################
# Combine generated results into single files for downstream analyses
###################################################
# combine all IDR2D results generated with pre-specified hyper-params into a single file
check_cached_results(
    results_file=FILTERED_IDR2D_RESULTS_FILE,
    force_redo=parsed.args$force_redo,
    # force_redo=TRUE,
    results_fnc=load_all_IDR2D_results
)
# count loops by differential status across all conditions + hyper-params, combine into a single file
check_cached_results(
    results_file=ALL_IDR2D_COUNTS_RESULTS_FILE,
    force_redo=parsed.args$force_redo,
    # force_redo=TRUE,
    results_fnc=count_all_IDR2D_results
)

