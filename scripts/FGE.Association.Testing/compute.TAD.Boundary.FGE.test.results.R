######################################################################
# Dependencies
######################################################################
library(here)
here::i_am('scripts/FGE.Association.Testing/compute.TAD.Boundary.FGE.test.results.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(broom)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'FGE.Association.Testing/utils.enrichment.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    library(tidyverse)
    library(magrittr)
})
# options(future.globals.maxSize=2.5 * (1024 ** 3))

######################################################################
# Testing params
######################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('resolutions', 'threads', 'force'),
        has.positional=FALSE
    )
parsed.args$resolutions <- c(100, 50, 25) * 1e3
# bins within X bins of a HiC feature == feature
HiCFeatureRadius.bins <- 
    c(0, 1, 2, 3, 4) 
test.hyper.params.df <- 
    bind_rows(
        # Fisher test params
        expand_grid(
            HiF.bin.radius=HiCFeatureRadius.bins,
            FGE.classification.threshold=c(1, 5, 10, 20, 30, 40) # min threshold of FGEs per bin
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
# Load all TAD boundaries
######################################################################
# load TAD boundaries
TAD.boundaries.df <- 
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    # change end of TAD to represent the start position of last bin in the TAD 
    # instead of the end of the last bin
    mutate(end=end - resolution) %>% 
    pivot_all_TADs_to_boundaries() %>% 
    dplyr::rename('feature.start'=boundary.start) %>%
    # all boundaries are only 1 bin
    mutate(feature.end=feature.start + resolution)

######################################################################
# Join TAD boundaries to binwise FGE signals
######################################################################
# map all bins to each TAD
TAD.Boundary.FGE.signal.joint.df <- 
    list_all_binwise_signal_files() %>% 
    inner_join(
        TAD.boundaries.df %>% 
        # nest(features.df=-c(resolution, method, TAD.params, Sample.Group, chr)),
        nest(features.df=-c(resolution, TAD.method, TAD.params, TAD.metric, Sample.Group)),
        relationship='many-to-many',
        by=join_by(resolution)
        # by=join_by(resolution, chr)
    ) %>%
    filter(resolution %in% parsed.args$resolutions)

TAD.Boundary.FGE.signal.joint.df
TAD.CTCF.overlaps.df %>% select(bins.df, features.df)
######################################################################
# Calculate cCRE Fisher Enrichment tests
######################################################################
check_cached_results(
    results_file=TAD_CCRE_FISHER_ENRICHMENTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=calculate_all_feature_CTCF_enrichments,
    overlaps.df=
        TAD.cCRE.overlaps.df %>% 
        cross_join(FISHER_PARAMS),
    test.type='fisher',
    # when correcting pvalues, group tests by these columns, to avoid over-correction
    p.corr.group.cols=
        c(
            # 'Sample.Group',
            # 'chr',
            'resolution',
            'method'
        )
)

######################################################################
# Calculate cCRE T-tests
######################################################################
check_cached_results(
    results_file=TAD_CCRE_TTEST_ENRICHMENTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=calculate_all_feature_CTCF_enrichments,
    overlaps.df=
        TAD.cCRE.overlaps.df %>% 
        cross_join(TTEST_PARAMS),
    test.type='t.test',
    # when correcting pvalues, group tests by these columns, to avoid over-correction
    p.corr.group.cols=
        c(
            # 'Sample.Group',
            # 'chr',
            'resolution',
            'method',
            'cCRE.type'
        )
)


######################################################################
# Join TAD boundaries to bin-wise CTCF data
######################################################################
# map all bins to each TAD
TAD.CTCF.overlaps.df <- 
# load bin-wise summary CTCF stats
    load_specific_overlap_results(
        # force.redo=TRUE,
        feature.type='binwise',
        annotation='CTCF'
    ) %>% 
    filter(resolution %in% RESOLUTIONS) %>% 
    dplyr::rename('motif'=annotation.type) %>% 
    select(-c(SampleID)) %>% 
    select(-ends_with(c('.var', '.q25', '.q75', '.min', '.max'))) %>% 
    # nest(bins.df=-c(annotation, feature.type, resolution, motif, chr)),
    nest(bins.df=-c(annotation, feature.type, resolution, motif)) %>% 
    inner_join(
        TAD.boundaries.df %>% 
        # nest(features.df=-c(resolution, method, TAD.params, Sample.Group, chr)),
        nest(features.df=-c(resolution, method, TAD.params, Sample.Group)),
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

######################################################################
# Calculate CTCF Fisher Enrichment tests
######################################################################
check_cached_results(
    results_file=TAD_CTCF_FISHER_ENRICHMENTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=calculate_all_feature_CTCF_enrichments,
    overlaps.df=
        TAD.CTCF.overlaps.df %>% 
        cross_join(FISHER_PARAMS),
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

######################################################################
# Calculate CTCF T-tests
######################################################################
check_cached_results(
    results_file=TAD_CTCF_TTEST_ENRICHMENTS_FILE,
    force_redo=parsed.args$force.redo,
    # force_redo=TRUE,
    results_fnc=calculate_all_feature_CTCF_enrichments,
    overlaps.df=
        TAD.CTCF.overlaps.df %>% 
        cross_join(TTEST_PARAMS),
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

