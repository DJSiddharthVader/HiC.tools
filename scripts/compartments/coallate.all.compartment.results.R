######################################################################
# Dependencies
######################################################################
library(here)
here::i_am('scripts/compartments/coallate.all.compartment.results.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    # library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
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

######################################################################
# Combine + compartmentalize all PC1 data to define compartments
######################################################################
# load all binwise PC1 data into a single file
binwise.df <- 
check_cached_results(
    results_file=ALL_COMPARTMENTS_RESULTS_FILE,
    results_fnc=load_all_cooltools_compartment_results,
    force_redo=parsed.args$force.redo,
    resolutions=parsed.args$resolutions,
    n.compartment.lvls.list=c(5, 10)
)

