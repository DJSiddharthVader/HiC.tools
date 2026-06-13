################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/compartments/make.saddle.data.cmds.R')
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
            normalization=c('balanced'),
            # normalization=c('balanced', 'raw'),
            contact.type=c('cis'),
            track.col.name=c('E1'),  # which track score to use for binning 
            expected.col.name=c('balanced.avg.smoothed'),
            qrange=c('0.02 0.98'),  # remove bins with outlier scores 
            n.bins=c(50)
            # n.bins=c(10, 50, 100)
        )
    ) %>%  
    cross_join(tibble(resolution=parsed.args$resolutions)) %>% 
    add_column(threads=parsed.args$threads)
# parsed.args$force.redo=TRUE
# parsed.args$force.redo=FALSE

################################################################################
# Generate cmds to call TADs with specified params
################################################################################
# First binify the genome and save bin coords to files
# Now generate commands to run cooltools eigs-cis to calculate + orient PC1 from the contact matrix
# We can bin the PC1 data to define compartment type + strength (i.e. Weak A, Strong B etc.)
hyper.params.df %>% 
    generate_all_saddle_data_calculation_cmds(
        merge_status='merged',
        force_redo=parsed.args$force.redo
    ) %>% 
    select(cmd) %>% 
    write_tsv(
        file.path(COMPARTMENTS_DIR, 'all.saddle.data.cmds.txt'),
        col_names=FALSE
    )

