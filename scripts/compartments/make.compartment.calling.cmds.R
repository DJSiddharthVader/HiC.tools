################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/compartments/make.compartment.calling.cmds.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    # library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
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
# All combinations of tool-specific hyper-params to call TADs with
hyper.params.df <- 
    bind_rows(
        # cooltools params
        expand_grid(
            track.type=c('genecov'),
            normalization=c('balanced'),
            # phasing.track=c('genecov', 'gene'),
            # normalization=c('balanced', 'raw'),
            Comp.method='cooltools'
        )
    ) %>%  
    cross_join(tibble(resolution=parsed.args$resolutions)) %>% 
    add_column(threads=parsed.args$threads)
# parsed.args$force.redo=TRUE

################################################################################
# Generate cmds to call TADs with specified params
################################################################################
# Now generate commands to run cooltools eigs-cis to calculate + orient PC1 from the contact matrix
# We can bin the PC1 data to define compartment type + strength (i.e. Weak A, Strong B etc.)
hyper.params.df %>% 
    generate_all_cooltools_calling_cmds(
        cmds.output.filepath=file.path(COMPARTMENTS_DIR, 'all.compartment.calling.cmds.txt'),
        merge_status='merged',
        force_redo=parsed.args$force.redo
    )

