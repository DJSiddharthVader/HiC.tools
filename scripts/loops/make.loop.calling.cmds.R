######################################################################
# Dependencies
######################################################################
library(here)
here::i_am('scripts/loops/make.loop.calling.cmds.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    # library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'coverage/utils.coverage.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    library(tidyverse)
    library(magrittr)
})

######################################################################
# Handle arguments/parameters
######################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )

######################################################################
# Generate cmds to call TADs with specified params
######################################################################
# All combinations of tool-specific hyper-params to call TADs with
hyper.params.df <- 
    bind_rows(
        # cooltools params
        expand_grid(
            # normalization=c('balanced', 'raw'),
            expected.col.name=c('balanced.avg.smoothed'),
            normalization=c('balanced'),
            loop.method='cooltools'
        ),
    ) %>%  
    cross_join(tibble(resolution=parsed.args$resolutions)) %>% 
    add_column(threads=parsed.args$threads)
# Now generated commands  to call TADs for each input matrix + param combo
# This will generate 1 txt file with 1 cmd per line, each cmd generates all the output files in 
# a named directory reflecting params/hyper-params used to generate those rsults
# Using the file as input, you can run all the cmds together with GNU parllel and/or via SLURM/SGE
hyper.params.df %>% 
    generate_all_loop_calling_cmds(
        cmds.output.filepath=file.path(LOOPS_DIR, 'all.loop.calling.cmds.txt'),
        force_redo=parsed.args$force.redo,
        merge_status='merged'
    )

