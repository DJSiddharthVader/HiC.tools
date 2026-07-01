######################################################################
# Depdendencies
######################################################################
library(here)
# here::i_am('scripts/loops/calculate.loop.valency.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(purrr)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.loading.results.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    library(magrittr)
    library(tidyverse)
})

######################################################################
# Set up all comparisons
######################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )

######################################################################
# Calculate Loop Valency for each loop anchor
######################################################################
# 2 group comparison + no covariates -> use exact test
message(glue('using {parsed.args$threads} core to parallelize'))
plan(multisession, workers=parsed.args$threads)
#  each row is 1 nested set of loop calls per condition + params + context
loops.df <- 
    # load loop results
    load_per_condition_loops() %>% 
    mutate(log10.qvalue=-log10(qvalue)) %>% 
    select(-c(qvalue)) %>% 
    nest(
        loops=
            c(
                FeatureID,
                anchor.left,
                anchor.right,
                count,
                length,
                enrichment,
                log10.qvalue
            )
    )
# Also calculate loop valency  i.e. how many loops each anchor is a part of
message('calculating loop valency...')
check_cached_results(
    results_file=ALL_LOOP_VALENCY_RESULTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=calculate_all_loop_valency,
    loops.df=loops.df
)

