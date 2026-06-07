################################################################################
# Load data for HiFs annotated in individual conditions
################################################################################
load_per_condition_TADs <- function(){
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    select(-starts_with('TAD.inner.')) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            TAD.method, TAD.params, TAD.metric,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='TAD')
}

load_per_condition_TAD_Boundaries <- function(){
    ALL_TAD_BOUNDARIES_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    filter(
        TAD.method == 'hiTAD' | 
        (TAD.method != 'hiTAD' & boundary.side == 'start')
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            TAD.method, TAD.params, TAD.metric,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='TAD.Boundary')
}

load_per_condition_loops <- function(force_redo=FALSE){
    check_cached_results(
        results_file=ALL_COOLTOOLS_LOOPS_RESULTS_FILE,
        force_redo=force_redo,
        results_fnc=load_all_cooltools_dots
    ) %>%
    # Filter and clean up loops
    post_process_cooltools_dots_results() %>%
    filter_loop_results() %>% 
    add_column(HiF.type='loop')
}

load_per_condition_loop_anchors <- function(){
    ALL_LOOP_VALENCY_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    add_column(HiF.type='loops.anchor')
}

load_per_condition_loop_nesting <- function(){
    stop('Not Implemented')
    ALL_LOOP_NESTING_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    add_column(HiF.type='loops.nesting')
}

load_per_condition_compartment_regions <- function(){
    ALL_COMPARTMENTS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    merge_bins_into_compartments() %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            Comp.method,
            normalization,
            track.type,
            score.source,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='compartment.region')
}

load_per_condition_compartment_switches <- function(){
    ALL_COMPARTMENTS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    filter(
        !is.na(compartment.switch),
        compartment.switch != 'no.switch',
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            Comp.method,
            normalization,
            track.type,
            score.source,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='compartment.switch')
}

load_specific_per_condition_HiFs <- function(HiF.name){
    case_when(
        HiF.name == 'TAD'                ~ list(load_per_condition_TADs()),
        HiF.name == 'TAD.Boundary'       ~ list(load_per_condition_TAD_Boundaries()),
        HiF.name == 'loop'               ~ list(load_per_condition_loops()),
        # HiF.name == 'loop.anchor'        ~ list(load_per_condition_loop_anchors()),
        # HiF.name == 'loop.nesting'       ~ list(load_per_condition_loop_nesting()),
        HiF.name == 'compartment.region' ~ list(load_per_condition_compartment_regions()),
        HiF.name == 'compartment.switch' ~ list(load_per_condition_compartment_switches()),
        # HiF.name == '' ~ list(load_per_condition_()),
        .unmatched='error'
    ) %>% 
    {.[[1]]}
}

combine_all_per_condition_HiFs <- function(
    HiFs.list,
    association.params.df){
    HiFs.list %>% 
    sapply(
        FUN=
            function(HiF.name, ...) {
                HiF.name %>% 
                load_specific_per_condition_HiFs() %>% 
                dplyr::rename('seqnames'=chr) %>% 
                nest(
                    HiFs.df=
                        -c(
                            HiF.type,
                            resolution,
                            SampleID
                        )
                )
            },
        simplify=FALSE
    ) %>% 
    bind_rows() %>% 
    # make all coordinate tibbles into irange objects to apply plyrange functions
    mutate(HiFs.df=pmap(.l=list(HiFs.df), .f=as_granges))
}

################################################################################
# Load data for differential HiFs annotated by comparing 2 conditions
################################################################################
load_between_condition_TADs <- function(){
    # nest so 1 set of TADs per row
    all.TADs.df <- 
        load_per_condition_TADs() %>% 
        nest(
            TADs.df=
                -c(
                    resolution,
                    TAD.method,
                    TAD.params,
                    TAD.metric,
                    normalization,
                    starts_with('Sample')
                )
        )
    # nest so 1 set of TADs per row
    diff.TAD.boundaries.df <- 
        load_between_condition_TAD_Boundaries() %>% 
        nest(
            diff.boundaries.df=
                -c(
                    resolution,
                    TADCompare.params,
                    TAD.method,
                    TAD.params,
                    normalization, 
                    starts_with('Sample')
                )
        )
    # map all differential boundaries to which TADs they are inside of or a boundary tp
    diff.TAD.boundaries.df %>% 
    mutate(comparison=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}')) %>% 
    # now every row is a set of diff TADs, 2 identical copies per comparison
    pivot_longer(
        c(Sample.Group.Numerator, Sample.Group.Denominator),
        names_to='Sample.Direction',
        names_prefix='Sample.Group.',
        values_to='Sample.Group'
    ) %>%
    # map diff tads to tads called in the numerator + denominator separatelty
    left_join(
        all.TADs.df,
        .,
        by=
            join_by(
                resolution, 
                normalization,
                TAD.method, TAD.params,
                Sample.Group
            )
    ) %>%
    # map diff Tads to TAds they are inside of in each conditions
    mutate(
        diff.TADs.df=
            pmap(
                .l=.,
                .f=
                    function(TADs.df, diff.boundaries.df, resolution, ...){
                        TADs.df %>% 
                        left_join(
                            diff.boundaries.df %>% 
                                select(chr, start, end, Gap.Score, p.adj.gw),
                            suffix=c('.TAD', '.diffBoundary'),
                            by=
                                join_by(
                                    chr,
                                    between(x$start, y$start, y$end)
                                )
                        ) %>%
                        group_by(FeatureID) %>%
                        slice_min(p.adj.gw, n=1) %>% 
                        select(-c(ends_with('.diffBoundary'))) %>% 
                        rename_with(~ str_remove(.x, '.TAD$'))
                    }
            )
    ) %>%
    unnest(diff.TADs.df)
}

load_between_condition_TAD_Boundaries <- function(force_redo=FALSE){
    check_cached_results(
        results_file=TADCOMPARE_RESULTS_FILE,
        force_redo=force_redo,
        # force_redo=TRUE,
        results_fnc=load_all_TADCompare_results
    ) %>% 
    # only keep comparisons for called boundaries
    # filter(isBoundary) %>%
    # select(-c(isBoundary)) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            TADCompare.params,
            # TAD.method, TAD.params, TAD.metric,
            TAD.method, TAD.params,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='TAD.Boundary')
}

load_between_condition_loops <- function(force_redo=FALSE){
    check_cached_results(
        results_file=ALL_IDR2D_RESULTS_FILE,
        force_redo=force_redo,
        # force_redo=TRUE,
        results_fnc=load_all_IDR2D_results
    ) %>% 
    post_process_IDR2D_results() %>% 
    filter_loop_IDR2D_results() %>% 
    add_column(HiF.type='loops')
}

load_between_condition_DIRs <- function(){
    check_cached_results(
        results_file=FILTERED_MULTIHICCOMPARE_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_multiHiCCompare_results,
        sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC,
        sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS,
        # resolutions=c(100, 50, 25, 10) * 1e3,
        resolutions=c(100, 50, 25) * 1e3,
        # resolutions=c(100, 50) * 1e3,
        gw.fdr.threshold=0.1,
        fdr.threshold=0.1,
        nom.threshold=0.05
    ) %>% 
    post_process_multiHiCCompare_results() %>% 
    add_column(HiF.type='DIRs')
}

load_between_condition_DIR_anchors <- function(){
    FILTERED_MULTIHICCOMPARE_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    post_process_multiHiCCompare_results() %>% 
    pivot_longer(
        c(region1 , region2),
        names_to='bin.pair.side',
        values_to='start'
    ) %>% 
    add_column(HiF.type='DIR.anchor')
}       

load_specific_between_condition_HiFs <- function(diff.HiF.name){
    case_when(
        diff.HiF.name == 'TAD'                ~ list(load_between_condition_TADs()),
        diff.HiF.name == 'TAD.Boundary'       ~ list(load_between_condition_TAD_Boundaries()),
        diff.HiF.name == 'loop'               ~ list(load_between_condition_loops()),
        diff.HiF.name == 'DIR'                ~ list(load_between_condition_DIRs()),
        diff.HiF.name == 'DIR.anchor'         ~ list(load_between_condition_DIR_anchors()),
        # HiF.name == 'compartment.region' ~ list(load_between_condition_compartment_regions()),
        # HiF.name == 'compartment.switch' ~ list(load_between_condition_compartment_switches()),
        # HiF.name == '' ~ list(load_per_condition_()),
        .unmatched='error'
    ) %>% 
    {.[[1]]}
}

combine_all_between_condition_HiFs <- function(
    diff.HiFs.list,
    association.params.df){
    diff.HiFs.list %>% 
    sapply(
        FUN=
            function(diff.HiF.name, ...) {
                diff.HiF.name %>% 
                load_specific_between_condition_HiFs() %>% 
                nest(
                    HiFs.df=
                        -c(
                            HiF.type,
                            resolution,
                            SampleID.Numerator, SampleID.Denominator
                        )
                )
            },
        simplify=FALSE
    ) %>% 
    bind_rows() %>% 
    mutate(HiFs.df=pmap(.l=list(HiFs.df), .f=as_iranges)) %>% 
    inner_join(
        association.params.df,
        by=
            join_by(
                HiF.type,
                resolution
            )
    )
}

