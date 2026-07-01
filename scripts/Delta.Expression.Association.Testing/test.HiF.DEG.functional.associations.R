######################################################################
# Dependencies
######################################################################
library(here)
here::i_am('scripts/Delta.Expression.Association.Testing/test.HiF.DEG.functional.associations.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.association.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.testing.R'))
    library(tidyverse)
    library(magrittr)
})

######################################################################
# Testing params
######################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('resolutions', 'threads', 'force'),
        has.positional=FALSE
    )
# parsed.args$resolutions <- c(100, 50, 25) * 1e3
# bins within X bins of a HiC feature == feature
HiCFeatureRadius.bins <- 
    c(0, 1, 2, 3, 4) 
test.hyper.params.df <- 
    bind_rows(
        # Fisher test params
        expand_grid(
            HiF.bin.radius=HiCFeatureRadius.bins,
        ) %>%
        add_column(test.type='fisher.test'),
        # T-test params
        tibble(
            HiF.bin.radius=HiCFeatureRadius.bins
        ) %>%
        add_column(test.type='t.test'),
        # Corr test params
        expand_grid(
            max.bin.dist.to.nearest.HiF=c(Inf, 200, 100, 50, 25),
        ) %>%
        add_column(test.type='corr.test'),
    )

######################################################################
# Load Gene & HiF annotations
######################################################################
all.HiF.DEG.mappings.df <- 
    HIF_GENE_ASSOCIATION_MAPPING_DIR %>%
    parse_results_filelist(
        filename.column.name='SampleID',
        suffix='-HiF.Gene.Associations.tsv'
    ) %>%
    mutate(
        HiF.DEG.associations.df=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                show_col_types=FALSE,
                progress=FALSE,
                .progress=TRUE
            )
    ) %>%
        {.} -> tmp; tmp
        tmp %>%
        cross_join(tibble(test.type=c('fisher'))) %>% 
        cross_join(tibble(nearby.dist.thresh.bins=c(1))) %>% 
        # cross_join(tibble(nearby.dist.thresh.bins=c(1, 2, 3))) %>% 
    calculate_all_association_test_results(
        p.corr.group.cols=
            c(
                'association.type',
                'association.source',
                'resolution',
                'HiF.type',
                'HiF.scope'
            )
    ) %>%
    select(association.subtype, resolution, SampleID, test.p.value)

