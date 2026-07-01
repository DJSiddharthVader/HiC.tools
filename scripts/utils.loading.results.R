######################################################################
# Data utils
######################################################################
filter_loop_results <- function(
    results.df,
    q.thresh=0.1){
    results.df %>% 
    filter(type == 'cis') %>% 
    filter(kernel == 'donut') %>% 
    filter(qvalue < q.thresh)
}

######################################################################
# Load mcool files
######################################################################
list_all_mcool_files <- function(
    pattern='.mcool',
    merge_status='merged',
    only_use_included_samples=TRUE, 
    rm_metadata_cols=TRUE,
    ...){
    # pattern='.mcool'; only_use_included_samples=TRUE; rm_metadata_cols=TRUE;
    # List all .mcool files
    MCOOL_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    mutate(MatrixID=str_remove(basename(filepath), pattern)) %>% 
    # parse sample metadata from filename
    parse_metadata_from_names(
        info.format='MatrixID',
        include_merged_col=TRUE,
        keep_id=FALSE
    ) %>% 
    build_name_from_metadata(info.format='SampleID') %>% 
    build_name_from_metadata(info.format='Sample.Group') %>% 
    build_name_from_metadata(info.format='MatrixID') %>% 
    filter(ReadFilter == 'mapq_30') %>% 
    {
        case_when(
            merge_status == 'both'       ~ list(.),
            merge_status == 'merged'     ~ list(filter(., isMerged == 'Merged')),
            merge_status == 'individual' ~ list(filter(., isMerged == 'Individual')),
            .unmatched='error'
        )
    } %>% 
    {.[[1]]} %>% 
    {
        if (only_use_included_samples){
            filter_included_samples(df=.)
        } else {
            .
        }
    } %>% 
    {
        if (rm_metadata_cols){
            select(., !all_of(intersect(colnames(.), ALL.METADATA.FIELDS)))
        } else {
            .
        }
    }
}

load_mcool_file <- function(
    filepath,
    resolution,
    normalization,
    range1="",
    range2="",
    cis=TRUE,
    type='df',
    include_ends=FALSE,
    ...){
    if (normalization %in% c('weight', 'balanced', 'ICE')){
        normalization <- "weight"
    } else {
        normalization <- "NONE"
    }
    if (range2 == "") {
        range2 <- range1
    }
    filepath %>% 
    File(resolution=resolution) %>% 
    fetch(
        range1=range1,
        range2=range2,
        normalization=normalization,
        type=type,
        join=TRUE,
        query_type='UCSC'
    ) %>% 
    # format column names
    {
        if (type == 'df') {
            if (cis) {
                as_tibble(.) %>%
                filter(chrom1 == chrom2) %>% 
                {
                    if (include_ends){
                        dplyr::rename(
                            .,
                            'chr.A'=chrom1,
                            'start.A'=start1,
                            'end.A'=end1,
                            'chr.B'=chrom2,
                            'start.B'=start2,
                            'end.B'=end2,
                            'IF'=count
                        )
                    } else {
                        select(
                            .,
                            -c(
                                chrom2,
                                end1,
                                end2
                            )
                        ) %>% 
                        dplyr::rename(
                            'chr'=chrom1,
                            'range1'=start1,
                            'range2'=start2,
                            'IF'=count
                        )
                    }
                }
            } else {
                as_tibble(.) %>%
                dplyr::rename(
                    'chr1'=chrom1,
                    'chr2'=chrom2,
                    'range1'=start1,
                    'range2'=start2,
                    'IF'=count
                )
            }
        } else {
            as.matrix(.)
        }
    }
}

load_mcool_files <- function(
    hic.params.df,
    merge_status='merged',
    regions.df=NULL,
    range1s=NULL,
    range2s=NULL,
    progress=TRUE,
    ...){
    # hic.params.df=expand.grid(resolution=c(100) * 1e3, normalization='NONE'); pattern='.mcool'; regions.df=NULL; range1s=NULL; range2s=NULL; progress=TRUE; keep_metadata_columns=FALSE;
    # Define all genomic regions to load contacts 
    regions.df <- 
        {
            if (is.null(regions.df)) {
                # Get the whole genome
                if ((is.null(range1s)) & (is.null(range2s))) {
                    tibble()
                    # tibble(
                    #     range1=CHROMOSOMES,
                    #     range2=CHROMOSOMES
                    # )
                # get all contacts within all regions in range1s
                } else if (is.null(range2s)) {
                    tibble(
                        range1=range1s,
                        range2=range1s
                    )
                # get all contacts between all pairs of regions only (not intra-region contacts)
                } else {
                    expand_grid(
                        range1=range1s,
                        range2=range2s
                    )
                }
            # Just load the specified regions (1 region per row: chr, start, end)
            } else {
                regions.df
            }
        }
    # List all regions for all samples
    list_all_mcool_files(merge_status=merge_status) %>% 
    join_all_rows(regions.df) %>% 
    # must contain 2 columns: resolution (int) and normalization (passed to load_mcool_file())
    # all samples will be loaded per each pair of resolution+normalization listed i.e. per row in hic.params.df
    cross_join(hic.params.df) %>% 
    # Load contacts if specified or just return sample metadata + filepaths + regions
    mutate(
        .,
        contacts=
            purrr::pmap(
                .l=.,
                .f=load_mcool_file,
                .progress=progress
            )
    ) %>% 
    select(-c(filepath, range1, range2)) %>% 
    unnest(contacts)
}

######################################################################
# Load data for HiFs annotated in individual conditions
######################################################################
load_per_condition_TADs <- function(){
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    select(-starts_with('TAD.inner.')) %>% 
    unite(
        TAD.params,
        sep='#',
        remove=FALSE,
        c(normalization, TAD.method, TAD.params)
        # c(normalization, TAD.method, TAD.params, TAD.metric)
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            TAD.params,
            chr, start, end
        )
    )
}

load_per_condition_TAD_Boundaries <- function(){
    load_per_condition_TADs() %>% 
    mutate(end=start + resolution) %>% 
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
    select(
        resolution,
        TAD.method,
        SampleID, Sample.Group,
        FeatureID, chr, start, end,
        TAD.length, TAD.start.score
    )
}

load_per_condition_loops <- function(force_redo=FALSE){
    check_cached_results(
        results_file=ALL_LOOP_RESULTS_FILE,
        force_redo=force_redo,
        results_fnc=load_all_cooltools_dots
    ) %>%
    # Filter and clean up loops
    filter_loop_results() %>% 
    unite(
        Loop.params,
        sep='#',
        remove=FALSE,
        c(method, type, normalization, kernel)
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            Loop.params,
            chr, anchor.left, anchor.right
        )
    )
}

load_per_condition_loop_anchors <- function(){
    ALL_LOOP_VALENCY_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            method, type,
            normalization,
            kernel,
            chr, anchor.position
        )
    )
}

load_per_condition_loop_nesting <- function(force_redo=FALSE){
    ALL_LOOP_NESTING_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE, progress=FALSE) %>%
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            chr, start, end
        )
    )
}

load_per_condition_compartment_regions <- function(force_redo=FALSE){
    ALL_COMPARTMENT_REGIONS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
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
    )
}

load_per_condition_compartment_switches <- function(){
    ALL_COMPARTMENT_SWITCHES_FILE %>% 
    # ALL_COMPARTMENT_SWITCHES_AND_CONTEXT_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
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
    )
}

load_specific_per_condition_HiFs <- function(HiF.name){
    case_when(
        HiF.name == 'TAD'                ~ 
            load_per_condition_TADs() %>% 
            unite(
                TAD.Params,
                sep='#',
                c(normalization, TAD.metric)
            ) %>% 
            dplyr::rename('length'=TAD.length) %>% 
            select(
                resolution, 
                TAD.Params,
                SampleID,
                FeatureID, chr,start, end,
                length, TAD.start.score
            ) %>% 
            pivot_longer(
                c(TAD.start.score),
                names_to='feature',
                values_to='value'
            ) %>% 
            list(),
        # HiF.name == 'TAD.Boundary'       ~ 
        #     load_per_condition_TAD_Boundaries() %>% 
        #     list(),
        HiF.name == 'loop'               ~ 
            load_per_condition_loops() %>%
            unite(
                Loop.Params,
                sep='#',
                c(method, type, normalization, kernel)
            ) %>% 
            dplyr::rename('start'=anchor.left, 'end'=anchor.right) %>% 
            select(
                resolution, 
                Loop.Params,
                SampleID,
                FeatureID, chr, start, end,
                length, enrichment, qvalue
            ) %>% 
            pivot_longer(
                c(enrichment, qvalue),
                names_to='feature',
                values_to='value'
            ) %>% 
            list(),
        # HiF.name == 'loop.anchor'        ~ 
        #     load_per_condition_loop_anchors() %>% 
        #     list(),
        HiF.name == 'loop.nesting'       ~ 
            load_per_condition_loop_nesting() %>% 
            list(),
        HiF.name == 'compartment.region' ~
            load_per_condition_compartment_regions() %>% 
            select(
                resolution, 
                compartment.params,
                SampleID,
                FeatureID, chr,start, end, length,
                compartment, score_mean, score_total
            ) %>% 
            pivot_longer(
                c(score_mean, score_total),
                names_to='feature',
                values_to='value'
            ) %>% 
            list(),
        # HiF.name == 'compartment.switch' ~ 
        #     load_per_condition_compartment_switches() %>% 
        #     list(),
        .unmatched='error'
    ) %>% 
    {.[[1]]} %>% 
    add_column(HiF.type=HiF.name)
}

combine_all_per_condition_HiFs <- function(HiFs.list){
    HiFs.list %>% 
    sapply(
        FUN=
            function(HiF.name, ...) {
                HiF.name %>% 
                load_specific_per_condition_HiFs() %>% 
                { if ('chr' %in% colnames(.)) { dplyr::rename(., 'seqnames'=chr) } else { . } } %>% 
                nest(
                    HiFs.df=
                        -c(
                            HiF.type,
                            resolution,
                            SampleID
                        )
                ) %>%
                # make all coordinate tibbles into irange objects to apply plyrange functions
                mutate(HiFs.df=pmap(.l=list(HiFs.df), .f=as_granges))
            },
        simplify=FALSE
    ) %>% 
    bind_rows()
}

######################################################################
# Load data for differential HiFs annotated by comparing 2 conditions
######################################################################
load_between_condition_TADs <- function(){
    # nest so 1 set of TADs per row
    all.TADs.df <- 
        load_per_condition_TADs() %>% 
        select(
            resolution,
            Sample.Group,
            TAD.params,
            chr, start, end, TAD.length,
            TAD.start.score
        ) %>% 
        dplyr::rename('length'=TAD.length, 'seqnames'=chr) %>% 
        nest(
            TADs.df=
                -c(
                    resolution,
                    TAD.params,
                    starts_with('Sample')
                )
        )
    # nest so 1 set of TADs per row
    diff.TAD.boundaries.df <- 
        load_between_condition_TAD_Boundaries() %>% 
        dplyr::rename('seqnames'=chr) %>% 
        select(
            resolution, TAD.params, TADCompare.params, 
            Sample.Group.Numerator, Sample.Group.Denominator,
            seqnames, start, end, Gap.Score, p.adj.gw
        ) %>% 
        nest(
            diff.boundaries.df=
                -c(
                    resolution,
                    TAD.params,
                    TADCompare.params,
                    starts_with('Sample')
                )
        )
    # map all differential boundaries to which TADs they are inside of or a boundary tp
    # diff.TAD.boundaries.df; all.TADs.df
    diff.TAD.boundaries.df %>% 
    mutate(comparison=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}')) %>% 
    # now every row is a set of diff TADs, 2 identical copies per comparison
    pivot_longer(
        c(Sample.Group.Numerator, Sample.Group.Denominator),
        names_to='Sample.Direction',
        names_prefix='Sample.Group.',
        values_to='Sample.Group'
    ) %>%
    separate_wider_delim(
        comparison,
        delim=' vs ',
        names=c('Sample.Group.Numerator', 'Sample.Group.Denominator')
    ) %>% 
    # map diff tads to tads called in the numerator + denominator separatelty
    left_join(
        all.TADs.df,
        .,
        by=
            join_by(
                resolution, 
                TAD.params,
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
                        # TADs.df %>% 
                        all.TADs.df$TADs.df[[1]] %>% 
                        as_granges() %>% 
                        mutate(
                            diff.boundaries.within=
                                count_overlaps(
                                    .,
                                    # diff.boundaries.df %>% 
                                    diff.TAD.boundaries.df$diff.boundaries.df[[1]] %>% 
                                        as_granges(),
                                    # minoverlap=resolution
                                    minoverlap=resolution
                                )
                        ) %>%
                        as_tibble()
                    }
            )
    ) %>%
    select(-c(TADs.df, diff.boundaries.df)) %>% 
    unnest(diff.TADs.df) %>%
    mutate(is.differential=diff.boundaries.within >= 1) %>% 
    select(-c(strand, width))
}

load_between_condition_TAD_Boundaries <- function(force_redo=FALSE){
    check_cached_results(
        results_file=TADCOMPARE_RESULTS_FILE,
        force_redo=force_redo,
        # force_redo=TRUE,
        results_fnc=load_all_TADCompare_results,
        gw.fdr.threshold=0.1,
        nom.threshold=0.05,
    ) %>% 
    # only keep comparisons for called boundaries
    # filter(isBoundary) %>%
    # select(-c(isBoundary)) %>% 
    unite(
        TAD.params,
        sep='#',
        c(normalization, TAD.method, TAD.params)
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            TADCompare.params,
            TAD.params,
            chr, start, end
        )
    )
}

load_between_condition_loops <- function(force_redo=FALSE){
    check_cached_results(
        results_file=FILTERED_IDR2D_RESULTS_FILE,
        force_redo=force_redo,
        results_fnc=load_all_IDR2D_results
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            method, type,
            normalization,
            kernel,
            chr, anchor.left, anchor.right
        )
    )
}

load_between_condition_loop_nesting <- function(){
    ALL_LOOP_NESTING_DIFFERENCE_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE, progress=FALSE) %>% 
    # ALL_IDR2D_NESTING_RESULTS_FILE %>% 
    # read_tsv(show_col_types=FALSE, progress=FALSE) %>% 
    # ALL_IDR2D_NESTING_DIFFERENCE_RESULTS_FILE %>% 
    # read_tsv(show_col_types=FALSE, progress=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            chr, start,end 
        )
    )
}

load_loop_nesting_by_reproducibility <-  function(){
    ALL_IDR2D_NESTING_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE, progress=FALSE) %>% 
    # ALL_IDR2D_NESTING_DIFFERENCE_RESULTS_FILE %>% 
    # read_tsv(show_col_types=FALSE, progress=FALSE) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            IDR2D.Params,
            chr, start,end 
        )
    )
}

# load_bewteen_compartment_regions() <- function(){
# }

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
    # Create a unique ID for each bin tested, to check overlaps across experiments
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            zero.p, A.min, merged,
            chr, region1, region2
        )
    )
}

load_between_condition_DIR_anchors <- function(){
    FILTERED_MULTIHICCOMPARE_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    pivot_longer(
        c(region1 , region2),
        names_to='bin.pair.side',
        values_to='start'
    ) %>% 
    unite(
        FeatureID,
        sep='#',
        remove=FALSE,
        c(
            zero.p, A.min, merged,
            chr, start, bin.pair.side
        )
    )
}       

load_specific_between_condition_HiFs <- function(diff.HiF.name){
    case_when(
        diff.HiF.name == 'TAD'~
            load_between_condition_TADs() %>% 
            dplyr::rename(
                'Numerator'=Sample.Group.Numerator,
                'Denominator'= Sample.Group.Denominator
            ) %>% 
            pivot_longer(
                c(TAD.start.score, diff.boundaries.within),
                names_to='feature',
                values_to='value'
            ) %>% 
            list(),
        # diff.HiF.name == 'TAD.Boundary'       ~ list(load_between_condition_TAD_Boundaries()),
        diff.HiF.name == 'loop' ~
            load_between_condition_loops() %>% 
            dplyr::rename('start'=anchor.left, 'end'=anchor.right) %>% 
            mutate(length=end - start) %>% 
            dplyr::rename(
                'Numerator'=SampleID.Numerator,
                'Denominator'= SampleID.Denominator
            ) %>% 
            select(
                resolution, 
                IDR2D.Params,
                Numerator, Denominator,
                chr, start, end, length,
                loop.status,
                diff.value, diff.rank, IDR
            ) %>% 
            pivot_longer(
                c(diff.value, diff.rank, IDR),
                names_to='feature',
                values_to='value'
            ) %>% 
            list(),
        diff.HiF.name == 'loop.nesting' ~ 
            load_between_condition_loop_nesting() %>% 
            filter(!str_detect(feature.stat, 'min|var|enrichment')) %>% 
            filter(!str_detect(diff.stat, 'var')) %>% 
            mutate(feature=glue('{feature.stat}_{diff.stat}')) %>% 
            mutate(length=n.bins * resolution) %>% 
            dplyr::rename(
                'Numerator'=SampleID.Numerator,
                'Denominator'= SampleID.Denominator
            ) %>% 
            select(
                resolution, 
                # SampleID.Numerator, SampleID.Denominator,
                Numerator, Denominator,
                SegmentID, chr, start, end, length,
                feature.stat, diff.stat, feature, value
            ) %>% 
            list(),
        diff.HiF.name == 'DIR' ~ 
            load_between_condition_DIRs() %>% 
            filter(merged == 'Individual') %>% 
            dplyr::rename(
                'start'=region1,
                'end'=region2,
                'length'=distance.bp
            ) %>% 
            dplyr::rename(
                'Numerator'=Sample.Group.Numerator,
                'Denominator'=Sample.Group.Denominator 
            ) %>% 
            select(
                resolution, merged,
                Numerator, Denominator,
                # Sample.Group.Numerator, Sample.Group.Denominator,
                FeatureID, chr, start, end, length,
                logFC, log.p.adj.gw
            ) %>% 
            list(),
        # HiF.name == 'compartment.region' ~ 
        #     load_between_condition_compartment_regions() %>% 
        #     list(),
        .unmatched='error'
    ) %>% 
    {.[[1]]} %>% 
    add_column(HiF.type=diff.HiF.name)
}

combine_all_between_condition_HiFs <- function(diff.HiFs.list){
    diff.HiFs.list %>% 
    sapply(
        FUN=
            function(diff.HiF.name, ...) {
                diff.HiF.name %>% 
                load_specific_between_condition_HiFs() %>% 
                { if ('chr' %in% colnames(.)) { dplyr::rename(., 'seqnames'=chr) } else { . } } %>% 
                # make all coordinate tibbles into irange objects to apply plyrange functions
                nest(
                    HiFs.df=
                        -c(
                            HiF.type,
                            resolution,
                            Numerator, Denominator
                        )
                ) %>% 
                mutate(HiFs.df=pmap(.l=list(HiFs.df), .f=as_granges))
            },
        simplify=FALSE
    ) %>% 
    bind_rows()
}

