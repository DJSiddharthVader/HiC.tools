######################################################################
# Dependencies
######################################################################
library(here)
BASE_DIR <- here()
suppressPackageStartupMessages({
    # library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.loading.results.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    library(tidyverse)
    library(magrittr)
    library(furrr)
})

######################################################################
# Handle arguments/parameters
######################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
# All combinations of tool-specific hyper-params to call TADs with
# parsed.args$resolutions=c(100, 50, 25) * 1e3
# parsed.args$force.redo=TRUE
plan(multisession, workers=length(availableWorkers()))

######################################################################
# Combine + compartmentalize all PC1 data to define compartments
######################################################################
# load all binwise PC1 data into a single file
binwise.df <- 
check_cached_results(
    results_file=ALL_COMPARTMENT_BINWISE_FILE,
    results_fnc=load_all_cooltools_compartment_results,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    resolutions=c(100, 50, 25) * 1e3,
    n.compartment.lvls=20
)
# squash binwise PC1 data into segments with A/B labels
check_cached_results(
    results_file=ALL_COMPARTMENT_REGIONS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    binwise.df=binwise.df,
    results_fnc=squash_all_bins_into_compartments
)
# define bins where a compartment switches from A->B or B->A and
# get position and summary stats i.e. difference in PC1 at switch
check_cached_results(
    results_file=ALL_COMPARTMENT_SWITCHES_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    binwise.df=binwise.df,
    results_fnc=
        function(binwise.df){
            binwise.df %>% 
            filter(does.compartment.switch)
        }
)
# same data but with +/- bin.context bins around it for switch analysis?
# check_cached_results(
#     results_file=ALL_COMPARTMENT_SWITCHES_AND_CONTEXT_FILE,
#     force_redo=parsed.args$force.redo,
#     # force_redo=TRUE,
#     binwise.df=binwise.df,
#     bin.context=10,
#     results_fnc=get_switches_and_context_from_all_bins
# )

