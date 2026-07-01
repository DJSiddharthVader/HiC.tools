################################################################################
# Depdendencies
################################################################################
library(here)
# here::i_am('scripts/TADs/coallate.all.TADCompare.results.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.Comparing.TADs.R'))
    library(tidyverse)
    library(magrittr)
})

################################################################################
# Handle arguments/parameters
################################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
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

################################################################################
# Parse  + load all TADCompare results
################################################################################
check_cached_results(
    results_file=TADCOMPARE_RESULTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=load_all_TADCompare_results,
    gw.fdr.threshold=0.1,
    nom.threshold=0.05
)

################################################################################
# Parse  + count number of significant differentail TAD Boundaries
################################################################################
check_cached_results(
    results_file=TADCOMPARE_COUNTS_RESULTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=load_correct_count_all_TADCompare_results
)

