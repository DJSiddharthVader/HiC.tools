################################################################################
# Fixed association analaysis parameter sets
################################################################################
HiF.association.strategies.df <- 
    tribble(
        ~HiF.type,            ~association.strategy,
        'DIR.anchors',        'nearby',
        'DIR',                'within',
        'compartment.switch', 'nearby',
        'compartment.region', 'within',
        'loop.nesting',       'within',
        'loop.anchor',        'within',
        'loop',               'within',
        'TAD.Boundary',       'nearby',
        'TAD',                'within'
    )
HiF.association.hyper.params.df <- 
    bind_rows(
        expand_grid(
            fuzzy.matching.threshold.bins=c(0, 1, 2, 3),
            association.strategy='nearby'
        ),
        expand_grid(
            frac.gene.matching.overlap=c(1.0, 0.5, 0.1),
            association.strategy='within'
        )
    ) %>% 
    inner_join(
        HiF.association.strategies.df,
        relationship='many-to-many',
        by=join_by(association.strategy)
    )

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
    # TADCOMPARE_RESULTS_FILE %>% 
    ALL_TAD_BOUNDARIES_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
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

load_per_condition_loops <- function(){
    stop('Not Implemented')
    ALL_COOLTOOLS_LOOPS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='loop')
}

load_per_condition_loop_anchors <- function(){
    stop('Not Implemented')
    ALL_LOOP_VALENCY_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='loops.anchor')
}

load_per_condition_loop_nesting <- function(){
    stop('Not Implemented')
    ALL_LOOP_NESTING_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='loops.nesting')
}

load_per_condition_compartment_regions <- function(){
    stop('Not Implemented')
    ALL_COMPARTMENTS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='compartment.region')
}

load_per_condition_compartment_switches <- function(){
    stop('Not Implemented')
    ALL_COMPARTMENTS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='compartment.switch')
}

load_specific_per_condition_HiFs <- function(HiF.name){
    case_when(
        HiF.name == 'TAD'                ~ list(load_per_condition_TADs()),
        HiF.name == 'TAD.Boundary'       ~ list(load_per_condition_TAD_Boundaries()),
        # HiF.name == 'loop'               ~ list(load_per_condition_loops()),
        # HiF.name == 'loop.anchor'        ~ list(load_per_condition_loop_anchors()),
        # HiF.name == 'loop.nesting'       ~ list(load_per_condition_loop_nesting()),
        # HiF.name == 'compartment.region' ~ list(load_per_condition_compartment_regions()),
        # HiF.name == 'compartment.switch' ~ list(load_per_condition_compartment_switches()),
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
                nest(
                    HiFs.df=
                        -c(
                            HiF.type,
                            resolution,
                            normalization,
                            SampleID
                        )
                )
            },
        simplify=FALSE
    ) %>% 
    bind_rows() %>% 
    # make all coordinate tibbles into irange objects to apply plyrange functions
    mutate(HiFs.df=pmap(.l=list(HiFs.df), .f=as_iranges)) %>% 
    # limit to relevant params
    inner_join(
        association.params.df,
        relationship='many-to-many',
        by=
            join_by(
                HiF.type,
                resolution
            )
    )
}

################################################################################
# Load data for differential HiFs annotated by comparing 2 conditions
################################################################################
load_between_condition_TADs <- function(){
    stop('Not Implemented')
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            TADCompare.params,
            TAD.method, TAD.params, TAD.metric,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='TAD')
}

load_between_condition_TAD_Boundaries <- function(){
    TADCOMPARE_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    # only keep comparisons for called boundaries
    filter(isBoundary) %>%
    select(-c(isBoundary)) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            TADCompare.params,
            TAD.method, TAD.params, TAD.metric,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='TAD.Boundary')
}

load_between_condition_loops <- function(){
    stop('Not Implemented')
    ALL_IDR2D_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='loops')
}

load_between_condition_DIRs <- function(){
    stop('Not Implemented')
    FILTERED_MULTIHICCOMPARE_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='DIRs')
}

load_between_condition_DIR_anchors <- function(){
    stop('Not Implemented')
    FILTERED_MULTIHICCOMPARE_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            normalization,
            chr, start, end
        )
    ) %>% 
    add_column(HiF.type='DIR.anchor')
}       

load_specific_between_condition_HiFs <- function(diff.HiF.name){
    case_when(
        diff.HiF.name == 'TAD'                ~ list(load_between_condition_TADs()),
        diff.HiF.name == 'TAD.Boundary'       ~ list(load_between_condition_TAD_Boundaries()),
        diff.HiF.name == 'loop'               ~ list(load_between_condition_loops()),
        diff.HiF.name == 'DIR'                ~ list(load_between_condition_DIRs()),
        # diff.HiF.name == 'DIR.anchor'         ~ list(load_between_condition_DIR_anchors()),
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

