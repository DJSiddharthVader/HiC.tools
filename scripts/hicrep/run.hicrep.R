######################################################################
# Depdendencies
######################################################################
library(here)
here::i_am('scripts/hicrep/run.hicrep.R')
BASE_DIR <- here()
# BASE_DIR <- '/data/talkowski/Samples/WAPL_NIPBL/HiC'
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
    source(file.path(SCRIPT_DIR, 'hicrep/utils.hicrep.R'))
    library(magrittr)
    library(tidyverse)
})

######################################################################
# Set up param combinations
######################################################################
options(scipen=999)
parsed.args <- 
    handle_CLI_args(
        args=c('force', 'resolutions'),
        has.positional=FALSE
    )
# dcHiC hyper-params
hyper.params.df <- 
    expand_grid(
        max.window.size=c(100000, 500000, 1000000, 5000000),
        resolution=parsed.args$resolutions
    ) %>%
    mutate(h=get_smoothing_param(resolution)) %>% 
    cross_join(
        tribble(
            ~is.downsampled.flag, ~is.downsampled,
            '',                  FALSE,
            '--bDownSample ',    TRUE
        )
    )

######################################################################
# Generate cli commands for all pairs of individual matrices
######################################################################
# List all pairs of matrices
list_mcool_files() %>%
    select(SampleID, filepath, isMerged) %>% 
    enumerate_pairwise_comparisons(
        SampleID.fields=c('Edit', 'Celltype', 'Genotype', NA, NA),
        sampleID_col='SampleID',
        suffixes=c('.P1', '.P2'),
        include_merged_col=FALSE,
        sample.group.comparisons=NULL,
        pair_grouping_cols=c('isMerged')
    ) %>% 
    # all pairs of matrices X all sets of hyper params
    cross_join(hyper.params.df) %>% 
    mutate(
        output_dir=
            file.path(
                HICREP_RESULTS_DIR,
                glue("resolution_{resolution}"),
                glue("h_{h}"),
                glue("is.downsampled_{is.downsampled}"),
                glue("window.size_{max.window.size}")
            ),
        output_file=
            file.path(
                output_dir,
                glue("{SampleID.P1}-{SampleID.P2}-hicrep.txt")
            ),
        redundant_file=
            file.path(
                output_dir,
                glue("{SampleID.P2}-{SampleID.P1}-hicrep.txt")
            )
    ) %>%
    mutate(
        across(
            c(filepath.P1, filepath.P2, output_dir, output_file, redundant_file),
            ~ str_replace(.x, BASE_DIR, '.')
        )
    ) %>% 
    {
        if (parsed.args$force.redo) {
            .
        } else {
            filter(., !(file.exists(output_file) | file.exists(redundant_file)))
        }
    } %>% 
    mutate(
        cmd=glue("hicrep {is.downsampled.flag}--dBPMax {max.window.size} --binSize {resolution} --h {h} {filepath.P1} {filepath.P2} {output_file} 2> /dev/null")
    ) %>% 
    select(cmd) %>% 
    write_tsv(
        file.path(HICREP_DIR, 'all.hicrep.cmds.txt'),
        col_names=FALSE
    )

