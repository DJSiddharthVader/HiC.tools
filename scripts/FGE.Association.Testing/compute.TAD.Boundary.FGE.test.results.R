###################################################
# Dependencis
###################################################
library(here)
library(broom)
here::i_am('scripts/functional.enrichment/test.TAD.boundary.CTCF.enrichment.R')
BASE_DIR <- here()
source(file.path(BASE_DIR,   'scripts/constants.R'))
source(file.path(SCRIPT_DIR, 'locations.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'utils.plot.R'))
source(file.path(SCRIPT_DIR, 'functional.enrichment/utils.enrichment.R'))
source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
# source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
library(tidyverse)
library(magrittr)
options(future.globals.maxSize=2.5 * (1024 ** 3))
# plan(sequential)
plan(multisession, workers=N_WORKERS_FOR_PARLLELIZATION)
RESOLUTIONS <- c(100, 50, 25) * 1e3
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )
FORCE_REDO <- parsed.args$force.redo

###################################################
# Join TAD boundaries to bin-wise CTCF data
###################################################
# load TAD boundaries
TAD.boundaries.df <- 
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    filter(resolution %in% RESOLUTIONS) %>% 
    # ignore wierd artifacts idk
    filter(TAD.length > 0) %>%
    # change end of TAD to represent the start position of last bin in the TAD 
    # instead of the end of the last bin
    mutate(end=end - resolution) %>% 
    pivot_TADs_to_boundaries() %>% 
    dplyr::rename('feature.start'=boundary.start) %>%
    # all boundaries are only 1 bin
    mutate(feature.end=feature.start + resolution)
# load bin-wise summary CTCF stats
functional.annotations.binwise.df <- 
    load_specific_overlap_results(
        # force.redo=TRUE,
        feature.type='binwise',
        annotation='CTCF'
    ) %>% 
    filter(resolution %in% RESOLUTIONS) %>% 
    dplyr::rename(
        'chr'=chrom,
        'motif'=annotation.type
    ) %>% 
    select(-c(SampleID)) %>% 
    select(-ends_with(c('.var', '.q25', '.q75', '.min', '.max')))
# map all bins to each TAD
TAD.CTCF.overlaps.df <- 
    TAD.boundaries.df %>% 
    # nest(features.df=-c(resolution, method, TAD.params, Sample.Group, chr)) %>% 
    nest(features.df=-c(resolution, method, TAD.params, Sample.Group)) %>% 
    inner_join(
        functional.annotations.binwise.df %>%
        # nest(bins.df=-c(annotation, feature.type, resolution, motif, chr)),
        nest(bins.df=-c(annotation, feature.type, resolution, motif)),
        relationship='many-to-many',
        by=join_by(resolution)
        # by=join_by(resolution, chr)
    ) %>%
    select(
        annotation, feature.type, motif,
        method, Sample.Group, 
        resolution,
        bins.df, features.df
    )

# TAD.CTCF.overlaps.df %>% select(bins.df, features.df)
###################################################
# Calculate Fisher Enrichment tests
###################################################
check_cached_results(
    results_file=TAD_CTCF_FISHER_ENRICHMENTS_FILE,
    force_redo=FORCE_REDO,
    # force_redo=TRUE,
    results_fnc=calculate_all_feature_CTCF_enrichments,
    overlaps.df=
        TAD.CTCF.overlaps.df %>% 
        cross_join(
            expand_grid(
                HiCFeatureRadius.bins=c(0, 1, 2, 3), # bins within X a feature == feature
                n.CTCF.min.thresh=c(1, 5, 10, 20, 30, 40) # min number of CTCF sites to classify bins
            )
        ),
    test.type='fisher',
    # when correcting pvalues, group tests by these columns, to avoid over-correction
    p.corr.group.cols=
        c(
            # 'Sample.Group',
            # 'chr',
            'resolution',
            'method',
            'motif'
        )
)

###################################################
# Calculate T-tests
###################################################
check_cached_results(
    results_file=TAD_CTCF_TTEST_ENRICHMENTS_FILE,
    force_redo=FORCE_REDO,
    # force_redo=TRUE,
    results_fnc=calculate_all_feature_CTCF_enrichments,
    overlaps.df=TAD.CTCF.overlaps.df,
    test.type='t.test',
    # when correcting pvalues, group tests by these columns, to avoid over-correction
    p.corr.group.cols=
        c(
            # 'Sample.Group',
            # 'chr',
            'resolution',
            'method',
            'motif'
        )
)

