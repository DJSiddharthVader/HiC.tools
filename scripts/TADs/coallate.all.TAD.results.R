################################################################################
# Depdendencies
################################################################################
library(here)
here::i_am('scripts/TADs/coallate.all.TAD.results.R')
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
# only include results with certain param values in combined files
inclusion.params.df <- 
    tribble(
        ~TAD.method, ~TAD.params, ~TAD.metric,
        # 'ConsensusTAD', '3#15#0.2',   'Consensus.Score',
        # 'cooltools',    '60#0.66#Li', 'log2_insulation_score',
        'hiTAD',        NA,           'ADI'
    ) %>% 
    cross_join(
        expand_grid(
            resolution=parsed.args$resolutions,
            normalization=c('balanced')
        )
    )

################################################################################
# Combine all TAD results 
################################################################################
# combine all TAD results into a single table with metadata i.e.
# 1 TAD per row i.e. per condition + resolution + TAD Calling parameters
# hiTAD results 
all.hiTAD.TADs.df <- 
    check_cached_results(
        results_file=HITAD_TAD_RESULTS_FILE,
        results_fnc=load_all_hiTAD_TADs,
        # force_redo=TRUE
        force_redo=parsed.args$force.redo
    )
# cooltools results
all.cooltools.TADs.df <- 
    check_cached_results(
        results_file=COOLTOOLS_TAD_RESULTS_FILE,
        results_fnc=load_all_cooltools_insulation_TADs,
        # force_redo=TRUE
        force_redo=parsed.args$force.redo
    )
# ConsensusTAD results
all.ConsensusTAD.TADs.df <- 
    check_cached_results(
        results_file=CONSENSUSTAD_TAD_RESULTS_FILE,
        results_fnc=load_all_ConsensusTAD_TADs,
        # force_redo=TRUE
        force_redo=parsed.args$force.redo
    )
    # all.ConsensusTAD.TADs.df %>% count(TAD.bins) %>% arrange(desc(n))
# Combine all methods together 
bind_rows(
    all.ConsensusTAD.TADs.df,
    all.cooltools.TADs.df,
    all.hiTAD.TADs.df,
) %>% 
inner_join(
    inclusion.params.df,
    by=
        join_by(
            resolution,
            normalization,
            TAD.method, 
            TAD.params,
            TAD.metric
        )
) %>% 
mutate(TAD.idx=glue('{TAD.method}${TAD.params}${TAD.metric}$T{row_number()}')) %>% 
write_tsv(ALL_TAD_RESULTS_FILE)

################################################################################
# Save all calculated binwise scores used to call TADs for each method
################################################################################
# Now pivot combined TADs to boundaries file
ALL_TAD_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    pivot_all_TADs_to_boundaries() %>% 
    write_tsv(ALL_TAD_BOUNDARIES_FILE)

################################################################################
# Combine all binwise scores used to call TADs 
################################################################################
# hiTAD results 
all.hiTAD.scores.df <- 
    check_cached_results(
        results_file=HITAD_SCORE_RESULTS_FILE,
        results_fnc=load_all_hiTAD_DIs,
        # force_redo=TRUE
        force_redo=parsed.args$force.redo
    )
# cooltools results
all.cooltools.scores.df <- 
    check_cached_results(
        results_file=COOLTOOLS_SCORE_RESULTS_FILE,
        results_fnc=load_all_cooltools_insulation_scores,
        # force_redo=TRUE
        force_redo=parsed.args$force.redo
    )
# ConsensusTAD results
all.ConsensusTAD.scores.df <- 
    check_cached_results(
        results_file=CONSENSUSTAD_SCORE_RESULTS_FILE,
        results_fnc=load_all_ConsensusTAD_scores,
        # force_redo=TRUE
        force_redo=parsed.args$force.redo
    )
# Combine all methods together 
bind_rows(
    all.ConsensusTAD.scores.df,
    all.cooltools.scores.df,
    all.hiTAD.scores.df
) %>% 
inner_join(
    inclusion.params.df,
    by=
        join_by(
            resolution,
            normalization,
            TAD.method, 
            TAD.params,
            TAD.metric
        )
) %>% 
select(-c(isConsensusBoundary, window.size, threshold, mfvp)) %>% 
write_tsv(ALL_TAD_SCORES_FILE)

