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
    # list_all_per_condition_gene_association_HiF_mapping_results() %>% 
    HIF_GENE_ASSOCIATION_MAPPING_DIR %>%
    parse_results_filelist(
        suffix='-HiF.Gene.Associations.tsv',
        filename.column.name='SampleID'
    ) %>%
    filter(
        HiF.type %in% c(
            'loop',
            'compartment.region'
        )
    ) %>% 
        {.} -> tmp; tmp;
        tmp %>% 
        head(10) %>% 
    mutate(
        test.results=
            pmap(
                .l=.,
                .f=compute_all_test_results,
                .progress=TRUE
            )
    )
compute_all_test_results <- function(
    filepath,
    HiF.type,
    ...){
    # paste('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()', sep='; ')
    # row.index=1; filepath=tmp$filepath[[row.index]]; association.type=tmp$association.type[[row.index]]; association.subtype=tmp$association.subtype[[row.index]]; association.source=tmp$association.source[[row.index]]; association.strategy=tmp$association.strategy[[row.index]]; resolution=tmp$resolution[[row.index]]; HiF.type=tmp$HiF.type[[row.index]]; HiF.scope=tmp$HiF.scope[[row.index]]; SampleID=tmp$SampleID[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% select(-c(filepath)) %>% t()
    filepath %>% 
    read_tsv(show_col_types=FALSE, progress=FALSE) %>%
        {.} -> tmp2; tmp2
    tmp2 %>% 
    pivot_wider(names_from=feature.y, values_from=value.y)
    tmp2 %>% colnames()
    tmp2 %>% count(!is.na(FeatureID))
    # tmp2 %>% count(feature.x, feature.y)
    tmp2 %>% 
        select(
            # SampleID.Numerator, SampleID.Denominator,
            association.status,
            Sample.Group.Numerator, Sample.Group.Denominator,
            FeatureID, chr, start, end,
            compartment,
            ends_with('.DESeq2')
        ) %>%
        distinct() %>%
        count(padj.DESeq2 < 0.1, association.status, !is.na(FeatureID))
}
        cross_join(tibble(test.type=c('t.test'))) %>% 
        # cross_join(tibble(test.type=c('t.test', 'fisher'))) %>% 
        # cross_join(tibble(nearby.dist.thresh.bins=c(1))) %>% 
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

# all.between.condition.HiF.DEG.mappings.df <- 
#     list_all_bewteen_condition_gene_association_HiF_mapping_results() %>% 
#         {.}

