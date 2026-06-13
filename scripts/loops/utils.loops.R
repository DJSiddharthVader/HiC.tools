library(stringi)
library(furrr)
library(idr2d)
# library(plyranges)

###################################################
# cooltools dots
###################################################
generate_cooltools_dots_calling_cmds <- function(
    threads,
    normalization,
    resolution,
    SampleID,
    mcool.filepath,
    distance.expectation.filepath,
    expected.col.name,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("resolution_{resolution}"),
            glue("normalization_{normalization}")
        )
    # Create filenames
    expected.uri    <- glue("{distance.expectation.filepath}::{expected.col.name}")
    mcool.uri       <- glue("{mcool.filepath}::resolutions/{resolution}")
    output.filepath <- glue("{output_dir}/{SampleID}-dots.tsv")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    main.cmd  <- glue("cooltools dots {weight_flag} --nproc {threads} --output {output.filepath} {mcool.uri} {expected.uri}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=output.filepath,
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    main.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_loop_calling_cmds <- function(
    hyper.params.df,
    cmds.output.filepath=NULL,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # merge_status='merged';
    list_all_distance_expectation_files() %>% 
    # filter(type == 'cis') %>% 
    filter(contact.type == 'cis') %>% 
    inner_join(
        hyper.params.df,
        by=
            join_by(
                normalization,
                resolution
            )
    ) %>% 
    # list contacts matrices for all samples to generate compartments for
    inner_join(
        list_all_mcool_files(merge_status=merge_status) %>%
        dplyr::rename('mcool.filepath'=filepath),
        by=join_by(MatrixID)
    ) %>% 
    mutate(output_dir=file.path(LOOP_RESULTS_DIR, glue("loop.method_{loop.method}"))) %>% 
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(loop.method, ...) {
                        case_when(
                            loop.method == 'cooltools' ~ generate_cooltools_dots_calling_cmds(...),
                            .unmatched='error'
                        )
                    },
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    save_cmds_to_file(
        cmds.output.filepath=cmds.output.filepath,
        force_redo=force_redo
    )
}

list_all_cooltools_dots_results <- function(resolutions=NULL){
    LOOP_RESULTS_DIR %>% 
    parse_results_filelist(
        filename.column.name='MatrixID',
        suffix='-dots.tsv'
    ) %>%
    convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } 
}

load_cooltools_dots <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        show_col_types=FALSE,
        progress=FALSE
    )
}

load_all_cooltools_dots <- function(resolutions=NULL){
    list_all_cooltools_dots_results(resolutions=resolutions) %>% 
    mutate(
        loops=
            # pmap(
            future_pmap(
                .,
                load_cooltools_dots,
                .progress=TRUE
            )
    ) %>%
    unnest(loops) %>% 
    mutate(length=start2 - start1) %>% 
    # pivot so each lloop as one kernel estimate per row
    pivot_longer(
        starts_with('la_exp.'),
        names_to='metric',
        names_prefix='la_exp.',
        values_to='value'
    ) %>% 
    separate_wider_delim(
        metric,
        delim='.',
        names=c('kernel', 'metric')
    ) %>% 
    pivot_wider(
        names_from='metric',
        values_from='value'
    )
}

post_process_cooltools_dots_results <- function(results.df) {
    results.df %>%
    dplyr::rename(
        'chr'=chrom1,
        'anchor.left'=start1,
        'anchor.right'=start2,
        'normalization'=weight,
        'enrichment'=value,
        'qvalue'=qval
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
    ) %>% 
    select(-c(ends_with(c('1', '2')), cstart1 ,cstart2, c_label, c_size ,region, filepath))
}

filter_loop_results <- function(
    results.df,
    q.thresh=0.1){
    results.df %>% 
    filter(kernel == 'donut') %>% 
    filter(type == 'cis') %>% 
    filter(normalization == 'balanced') %>% 
    filter(kernel == 'donut') %>% 
    filter(qvalue < q.thresh)
}

###################################################
# IDR2D Analysis
###################################################
tidy_IDR2D_sided_results <- function(
    results.obj,
    metric_colname,
    ...){
    results.obj %>% 
    as_tibble() %>% 
    mutate(
        # idr=ifelse(is.na(idr), -1, idr), # to help identify rep-exclusive loops
        diff.value=value - rep_value,
        diff.rank=rank - rep_rank
    ) %>% 
    dplyr::rename(
        'chr'=chr_a,
        'anchor.left'=start_a,
        'anchor.right'=start_b,
        'IDR'=idr
    ) %>% 
    select(
        chr, anchor.left, anchor.right,
        diff.value, diff.rank,
        IDR
    )
}

tidy_IDR2D_results <- function(
    results,
    metric_colname,
    ...){
    # all loops from numerator
    tidy.loops.Numerator <- 
        tidy_IDR2D_sided_results(
            results$rep1_df,
            metric_colname
        )
    # all loops from denominator
    tidy.loops.Denominator <- 
        tidy_IDR2D_sided_results(
            results$rep2_df,
            metric_colname
        ) %>%
        # consistent sign so -ve => rep1 > rep2 for both sets of results
        mutate(
            across(
                starts_with('diff.'),
                ~ -.x
            )
        )
    # loops detected in both condtions with IDR stats
    tidy.common.loops <- 
        bind_rows(
            tidy.loops.Numerator,
            tidy.loops.Denominator
        ) %>%
        filter(!is.na(IDR)) %>% 
        distinct()
    # Bind all loops together into a single df
    # create column indicating which loops are reproduced or exclusive between conditions
    bind_rows(
        tidy.common.loops %>% 
        add_column(loop.type='detected.in.both'),
        # loops only detected in the numerator
        anti_join(
            tidy.loops.Numerator,
            tidy.common.loops
        ) %>% 
        mutate(IDR=1) %>% 
        add_column(loop.type='Numerator.only'),
        # loops only detected in the denominator
        anti_join(
            tidy.loops.Denominator,
            tidy.common.loops
        ) %>% 
        mutate(IDR=1) %>% 
        add_column(loop.type='Denominator.only'),
    ) %>%
    arrange(chr, anchor.left, anchor.right)
}

run_IDR2D_analysis <- function(
    loops.Numerator,
    loops.Denominator,
    metric_colname,
    value_transformation,
    ambiguity_resolution_method,
    max_gap,
    ...){
    # paste('row.index=1;', paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row.index]]", collapse='; '), ';t(head(tmp, 1))')
    loops.Numerator <- 
        loops.Numerator %>% 
        mutate(across(c(chr_A, chr_B), as.character)) %>% 
        mutate(across(c(start_A, start_B, end_A, end_B), as.integer)) %>% 
        select(
            chr_A, start_A, end_A,
            chr_B, start_B, end_B, 
            !!sym(metric_colname)
        ) %>%
        as.data.frame()
    loops.Denominator <-  
        loops.Denominator %>% 
        mutate(across(c(chr_A, chr_B), as.character)) %>% 
        mutate(across(c(start_A, start_B, end_A, end_B), as.integer)) %>% 
        select(
            chr_A, start_A, end_A,
            chr_B, start_B, end_B, 
            !!sym(metric_colname)
        ) %>%
        as.data.frame()
    # run IDR2D to define replicable loops between conditions
    estimate_idr2d(
        loops.Numerator,
        loops.Denominator, 
        value_transformation=value_transformation,
        ambiguity_resolution_method=ambiguity_resolution_method,
        max_gap=ifelse(max_gap < 1, -1L, max_gap)
    ) %>% 
    # tidy up results into nice tabular format
    tidy_IDR2D_results(
        loops.Numerator,
        loops.Denominator,
        metric_colname
    )
}

run_all_IDR2D_analysis <- function(
    all.loop.data.df,
    hyper.params.df,
    force.redo,
    sample.group.comparisons,
    pair_grouping_cols,
    ...){
    # force.redo=TRUE; sample.group.comparisons=comparisons.df; pair_grouping_cols=c('isMerged', 'method', 'kernel', 'type', 'normalization', 'resolution', 'chr')
    sample.group.comparisons %>% 
    left_join(
        all.loop.data.df,
        relationship='many-to-many',
        by=join_by(SampleID.Numerator == SampleID),
    ) %>% 
    left_join(
        all.loop.data.df %>% 
        dplyr::rename('SampleID.Denominator'=SampleID),
        relationship='many-to-many',
        suffix=c('.Numerator', '.Denominator'),
        by=c(pair_grouping_cols, 'SampleID.Denominator'),
    ) %>% 
    # Evaluate all comparisons for all combinations of specified parameters
    cross_join(hyper.params.df) %>% 
    mutate(
        max_gap=resolution * max_gap_bins,
        output_dir=
            file.path(
                LOOPS_IDR2D_DIR,
                glue('method_{method}'),
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('normalization_{normalization}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                glue('metric_{metric_colname}'),
                glue('resolve.method_{ambiguity_resolution_method}'),
                glue('max.gap_{max_gap}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID.Numerator}_vs_{SampleID.Denominator}-IDR2D.tsv')
            )
    ) %>% 
    filter(chr != 'chrY') %>% 
    arrange(desc(chr)) %>% 
    {
        if (!force.redo) {
            filter(., !file.exists(results_file))
        } else{
            .
        }
    } %T>% 
    {
        message('Generating the following results files')
        print(
            mutate(
                ., 
                comparison=glue('{SampleID.Numerator}-{SampleID.Denominator}')
            ) %>% 
            dplyr::count(
                # comparison
                metric_colname, ambiguity_resolution_method,
                resolution, max_gap_bins
            ) %>%
            dplyr::rename(
                'Metric'=metric_colname,
                # 'V.T'=value_transformation,
                'A.R.M'=ambiguity_resolution_method,
                'Max Gap'=max_gap_bins,
            )
        )
    } %>%
    future_pmap(
        .l=.,
        .f=check_cached_results,
        force_redo=force.redo,
        return_data=FALSE,
        results_fnc=run_IDR2D_analysis,
        .progress=TRUE
    )
}

list_all_IDR2D_results <- function(
    resolutions=NULL,
    ...){
    LOOPS_IDR2D_DIR %>% 
    parse_results_filelist(
        suffix='-IDR2D.tsv',
        filename.column.name='Sample.Group.Pair'
    ) %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } %>% 
    separate_wider_delim(
        Sample.Group.Pair,
        delim='_vs_',
        names=c('SampleID.Numerator', 'SampleID.Denominator')
    )
}

load_IDR2D_results <- function(filepath, ...){
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
    mutate(
        loop.status=
            case_when(
                loop.type == 'detected.in.both' & IDR >= 0.1 & IDR <= 1 ~ 'Irreproducible',
                loop.type == 'detected.in.both' & IDR < 0.1             ~ 'IDR < 0.1',
                loop.type == 'Numerator.only'                           ~ 'Numerator.only',
                loop.type == 'Denominator.only'                         ~ 'Denominator.only',
                TRUE                                                    ~ NA
            )
    ) %>%
    select(-c(loop.type))
}

load_all_IDR2D_results <- function(resolutions=NULL){
    list_all_IDR2D_results(resolutions=resolutions) %>% 
    mutate(max.gap.bins=max.gap / resolution) %>% 
    filter_loop_IDR2D_results() %>% 
    mutate(
        idr2d=
            # pmap(
            future_pmap(
                .l=.,
                .f=load_IDR2D_results,
                .progress=TRUE
            )
    ) %>%
    unnest(idr2d) %>% 
    select(-c(filepath, max.gap, region)) %>% 
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
    ) %>% 
    unite(
        IDR2D.Params,
        sep='#',
        remove=FALSE,
        c(metric, resolve.method, max.gap.bins)
    ) %>% 
    relocate(
        resolution,
        metric, resolve.method, max.gap.bins,
        IDR2D.Params,
        method, type, normalization, kernel,
        FeatureID,
        chr, anchor.left, anchor.right,
        diff.value, diff.rank, IDR, loop.status
    )
}

filter_loop_IDR2D_results <- function(results.df){
    results.df %>% 
    filter(type == 'cis') %>% 
    filter(normalization == 'balanced') %>% 
    filter(metric == 'log10.qvalue') %>% 
    filter(resolve.method == 'value') %>%
    filter(max.gap.bins == 5)
}

post_process_IDR2D_results <- function(results.df){
    results.df %>% 
    mutate(
        # max.gap.bins.int=max.gap / resolution,
        max.gap.bins=
            fct_reorder(
                glue('{max.gap.bins} bins'),
                max.gap.bins
            ),
        loop.status=
            factor(
                loop.status,
                levels=
                    c(
                        'Irreproducible',
                        'IDR < 0.1',
                        'Numerator.only',
                        'Denominator.only'
                    )
            )
    )
}

count_all_IDR2D_results <- function(resolutions=NULL){
    list_all_IDR2D_results(resolutions=resolutions) %>% 
    mutate(
        idr2d=
            # pmap(
            future_pmap(
                .l=.,
                .f=
                    function(filepath, ...){ 
                        filepath %>% 
                        load_IDR2D_results() %>% 
                        count(chr, loop.status)
                    },
                .progress=TRUE
            )
    ) %>%
    unnest(idr2d) %>% 
    mutate(max.gap.bins=max.gap / resolution) %>% 
    select(-c(filepath, max.gap, region))
}

join_loop_and_IDR2D_resilts <- function() {
    # Load all loop data to quantify nesting with
    all.loop.data.df <- 
        load_per_condition_loops() %>%
        mutate(log10.qvalue=-log10(qvalue)) %>% 
        pivot_longer(
            c(count, enrichment, log10.qvalue),
            names_to='loop.feature',
            values_to='loop.value'
        ) %>% 
        select(
            resolution, SampleID, 
            FeatureID,
            loop.feature, loop.value
        )
    # load differential loop results and map stats to loops
    load_between_condition_loops() %>%
    select(
        resolution, IDR2D.Params, 
        SampleID.Numerator, SampleID.Denominator,
        FeatureID, chr, anchor.left, anchor.right,
        loop.status
    ) %>%
    mutate(Comparison=glue('{SampleID.Numerator} vs {SampleID.Denominator}')) %>% 
    # pivot so I can easily nest so that
    # loops are per condition + have differential status per comparison
    pivot_longer(
        c(SampleID.Numerator, SampleID.Denominator),
        names_to='Comparison.side',
        names_prefix='SampleID.',
        values_to='SampleID'
    ) %>% 
    filter(
        (loop.status == 'Numerator.only'   & Comparison.side == 'Numerator'  ) |
        (loop.status == 'Denominator.only' & Comparison.side == 'Denominator') |
        (!loop.status %in% c('Numerator.only', 'Denominator.only'))
    ) %>% 
    # add condition specific individual loop data for loop per comparison
    inner_join(
        all.loop.data.df,
        relationship='many-to-many',
        by=
            join_by(
                resolution, 
                FeatureID,
                SampleID
            )
    )
}

###################################################
# Valency Analysis
###################################################
calculate_loop_valency <- function(df) {
    df %>%
    pivot_longer(
        starts_with('anchor.'),
        names_to='anchor.side',
        names_prefix='anchor.',
        values_to='anchor.position'
    ) %>% 
    group_by(anchor.position) %>% 
    dplyr::summarize(
        valency=dplyr::n(),
        across(
            .cols=c(count, length, enrichment, log10.qvalue),
            .fns=list('mean'=mean, 'min'=min, 'max'=max, 'median'=median, 'var'=var, 'total'=sum),
            .names="metric_{.col}_{.fn}"
        )
    )
}

calculate_all_loop_valency <- function(
    loops.df,
    ...){
    # calculate how many loops each anchor is a part of
    loops.df %>%
    mutate(
        valency.results=
            pmap(
                .l=list(loops),
                .f=calculate_loop_valency,
                .progress=TRUE
            )
    ) %>%
    unnest(valency.results) %>%
    select(-c(loops))
}

###################################################
# Loop Nesting Results
###################################################
# input data for calculating nesting per condition + stratified by differential loop status
load_all_loop_data_for_nesting_analysis <- function(){
    # get all loops per condition and pivot metrics columns
    all.per.condition.loop.data.df <- 
        load_per_condition_loops() %>%
        mutate(log10.qvalue=-log10(qvalue)) %>% 
        select(-c(length, count)) %>% 
        pivot_longer(
            # c(length, count, enrichment, log10.qvalue),
            c(enrichment, log10.qvalue),
            names_to='loop.feature',
            values_to='loop.value'
        ) %>% 
        select(
            resolution, SampleID, 
            FeatureID, chr, anchor.left, anchor.right,
            loop.feature, loop.value
        )
    # load differential loop results and map stats to loops
    all.between.condition.loop.data.df <- 
        load_between_condition_loops() %>%
        select(
            resolution, IDR2D.Params, 
            SampleID.Numerator, SampleID.Denominator,
            FeatureID, chr, anchor.left, anchor.right,
            loop.status
        ) %>%
        mutate(Comparison=glue('{SampleID.Numerator} vs {SampleID.Denominator}')) %>% 
        # pivot so I can easily nest so that
        # loops are per condition + have differential status per comparison
        pivot_longer(
            c(SampleID.Numerator, SampleID.Denominator),
            names_to='Comparison.side',
            names_prefix='SampleID.',
            values_to='SampleID'
        ) %>% 
        filter(
            (loop.status == 'Numerator.only'   & Comparison.side == 'Numerator'  ) |
            (loop.status == 'Denominator.only' & Comparison.side == 'Denominator') |
            (!loop.status %in% c('Numerator.only', 'Denominator.only'))
        ) %>% 
        # add condition specific individual loop data for loop per comparison
        inner_join(
            all.per.condition.loop.data.df %>% 
                select(-c(chr, anchor.left, anchor.right)),
            relationship='many-to-many',
            by=
                join_by(
                    resolution, 
                    FeatureID,
                    SampleID
                )
        )
    # bind together since they are analyzed the same way, produces separate output files per differential loop status
    # and marginalzied over all loops
    all.per.condition.loop.data.df %>% 
    add_column(
        loop.status='all',
        Comparison.side='None',
        Comparison='None'
    ) %>% 
    bind_rows(all.between.condition.loop.data.df) %>% 
    # change end to be end of last bin, not start, since we are treating loops as segments to calc nesting
    # anchor.right is the bin start, so change it to bin end to capture that bin in each loop
    mutate(anchor.right=anchor.right + resolution) %>% 
    dplyr::rename(
        'start'=anchor.left,
        'end'=anchor.right
    ) %>% 
    nest(
        loops.df=
            c(
                FeatureID,
                chr, start, end, 
                loop.feature, loop.value, 
            )
    )
}
# compute all nesting stats i.e. summary stats over loops per bin, squashed into segments of continuous bins
compute_nesting_stats_per_bin <- function(
    loops.df,
    bins.df,
    ...){
    # paste('row.index=1', paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row.index]]", collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()', sep='; ')
    # row.index=1; resolution=tmp$resolution[[row.index]]; SampleID=tmp$SampleID[[row.index]]; loop.status=tmp$loop.status[[row.index]]; Comparison.side=tmp$Comparison.side[[row.index]]; Comparison=tmp$Comparison[[row.index]]; IDR2D.Params=tmp$IDR2D.Params[[row.index]]; loops.df=tmp$loops.df[[row.index]]; bins.df=tmp$bins.df[[row.index]]; results_file=tmp$results_file[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% t()
    # compute summary stats across all loops overlapping each bin
    # then collapse run of contiguous bins together (i.e. all bins within the exact same set of loops)
    bins.df %>%
    dplyr::rename('seqnames'=chr) %>%
    as_granges() %>% 
    # map all individual bins to all loops overlapping them
    join_overlap_inner_within(
        loops.df %>% 
        dplyr::rename('seqnames'=chr) %>%
        as_granges(),
        minoverlap=1L,
        suffix=c('.bin', '.loop')
    ) %>%
    as_tibble() %>% 
    select(-c(strand, width)) %>% 
    dplyr::rename('chr'=seqnames) %>% 
    # for each genomic bin + loop feature 
    group_by(
        chr, start, end,
        loop.feature
    ) %>%
    # calculate loop feature summary statistics + count nesting lvl (how many loops overlap each bin)
    summarize(
        nesting.lvl=n(),
        across(
            .cols=c(loop.value),
            .fn=
                list(
                    'var'=var,
                    'min'=min,
                    'mean'=mean,
                    'max'=max,
                    'total'=sum
                ),
            .names="metric_{.fn}"
        )
    ) %>%
    ungroup() %>% 
    # pivot so now each row is a single nested segment and all the stats are separate columns
    pivot_wider(
        names_from=loop.feature,
        names_glue='{.value}_{loop.feature}',
        values_from=starts_with('metric_')
    )
}

compute_loop_nesting_results <- function(
    loops.df,
    bins.df,
    ...){
    # For each genomic bin count how many loops overlap that bin
    # and summary stats of loop feature (e.g. pvalue) for those overlappign loops
    compute_nesting_stats_per_bin(
        loops.df,
        bins.df
    ) %>%
    # Now collapse all sets of contiguous bins at the same nesting lvl into segments
    # so each row is a segment (start-end) instead of a single bin
    # all bins squashed into the same segment are this way are overlapped by 
    # the same set of loops they all have the same summary stats
    group_by(across(c(starts_with('metric_'), 'nesting.lvl', 'chr'))) %>% 
    # Take the leftmost start and righmost end across all bins
    # per group as the segment start/end and just save the segment coords
    summarize(
        start=min(start),
        end=max(end)
    ) %>%
    ungroup() %>% 
    arrange(chr, start, end) %>% 
    relocate(chr, start, end, nesting.lvl)
}

compute_all_loop_nesting_results <- function(
    all.loop.data.df,
    all.bins.df,
    force_redo=FALSE){
    # for every set of loops join the set of all genomic bins at the same resolution
    all.loop.data.df %>%
    left_join(
        all.bins.df,
        by=join_by(resolution)
    ) %>% 
    # set up output filepath for each set of results
    mutate(
        # Comparison=str_replace(Comparison, ' ', '_'),
        # loop.status_=str_replace(loop.status, ' ', '_'),
        results_file=
            file.path(
                ALL_LOOP_NESTING_RESULTS_DIR,
                glue('resolution_{resolution}'),
                glue('IDR2D.Params_{IDR2D.Params}'),
                glue('Comparison_{Comparison}'),
                glue('Comparison.side_{Comparison.side}'),
                glue('loop.status_{loop.status}'),
                glue('{SampleID}-loop.nesting.results.tsv')
            )
    ) %>% 
        {.} -> tmp; tmp
        tmp %>% head(3) %>% 
    # compute and save all nested segments for each set of loop data
    future_pmap(
    # pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=compute_loop_nesting_results,
        force_redo=force_redo,
        .progress=TRUE
    )
}
# list and load all nesting results
list_all_loop_nesting_results <- function(){
    ALL_LOOP_NESTING_RESULTS_DIR %>% 
    parse_results_filelist(
        filename.column.name='SampleID',
        suffix='-loop.nesting.results.tsv'
    )
}

load_all_loop_nesting_results <- function(){
    list_all_loop_nesting_results() %>%
    mutate(
        nesting.results=
            pmap(
                 .l=list(filepath),
                 .f=read_tsv,
                 show_col_types=FALSE
             )
    ) %>%
    unnest(nesting.results)
}

post_process_loop_nesting_result <- function(results.df){
    results.df %>% 
    filter(kernel == 'donut') %>% 
    filter(normalization == 'balanced') %>% 
    mutate(length=end - start) %>% 
    pivot_longer(
        starts_with('metric_'),
        names_to='tmp',
        values_to='value'
    ) %>% 
    separate_wider_delim(
        tmp,
        delim='_',
        names=c(NA, 'loop.stat', 'loop.feature')
    )
}

###################################################
# Loop Nesting Correlation Analysis
###################################################
# compute rolling correlation of nesting stats across bins between a pair of conditions
compute_nesting_stat_rolling_correlations <- function(
    comparison.df,
    window.size){
    comparison.df %>%
    group_by(loop.feature) %>% 
    arrange(chr, start, end) %>% 
    cross_join(tibble(corr.metric=c('pearson', 'spearman', 'kendall'))) %>% 
    mutate(
        rolling.corr=
            cor(
                lead(loop.value.Numerator,   n=window.size),
                lead(loop.value.Denominator, n=window.size),
                method=corr.metric
            )
    ) %>%
    pivot_longer(everthing(), names_to='corr.stat', values_to='corr.value')
}

compute_nesting_correlation_results <- function(
    loops.df.Numerator,
    loops.df.Denominator,
    bins.df,
    window.size,
    ...){
    # paste('row.index=1;', paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row.index]]", collapse='; '), ';t(head(tmp, 1))')
    # row.index=1; resolution=tmp$resolution[[row.index]]; SampleID.Numerator=tmp$SampleID.Numerator[[row.index]]; loop.status=tmp$loop.status[[row.index]]; Comparison.side=tmp$Comparison.side[[row.index]]; Comparison=tmp$Comparison[[row.index]]; IDR2D.Params=tmp$IDR2D.Params[[row.index]]; loops.df.Numerator=tmp$loops.df.Numerator[[row.index]]; SampleID.Denominator=tmp$SampleID.Denominator[[row.index]]; loops.df.Denominator=tmp$loops.df.Denominator[[row.index]]; bins.df=tmp$bins.df[[row.index]]; window.size=tmp$window.size[[row.index]]; results_file=tmp$results_file[[row.index]] ;t(head(tmp, 1))
    # for each bin get nesting stats for numerator + denominator paired together
    comparison.df <- 
        # for each genomic bin, compute summary stats over all loops overlapping said bin 
        compute_nesting_stats_per_bin(
            loops.df.Numerator,
            bins.df
        ) %>% 
        pivot_longer(
            c(starts_with('metric_'), nesting.lvl),
            names_to='loop.feature',
            values_to='loop.value'
        ) %>% 
        # join numerator and denominator nesting stats 
        full_join(
            # for each genomic bin, compute summary stats over all loops overlapping said bin 
            compute_nesting_stats_per_bin(
                loops.df.Denominator,
                bins.df
            ) %>% 
            pivot_longer(
                c(starts_with('metric_'), nesting.lvl),
                names_to='loop.feature',
                values_to='loop.value'
            ),
            suffix=c('.Numerator', '.Denominator'),
            by=join_by(chr, start, end, loop.feature)
        ) %>%
        # set nesting stats to 0 if any bin is only overlapped by loops in either numerator or denominator
        filter(grepl(paste('nesting.lvl', 'mean', 'max', 'total', sep='|'), loop.feature)) %>% 
        mutate(
            across(
                c(loop.value.Numerator, loop.value.Denominator),
                ~ ifelse(is.na(.x), 0, .x)
            )
        )
    # compute correlation of each metric across all bins per chr
    per.chr.stats <- 
        comparison.df %>% 
        group_by(chr, loop.feature) %>% 
        compute_nesting_stat_rolling_correlations(window.size=window.size)
    # compute correlation of each metric genome-wide
    gw.stats <- 
        comparison.df %>% 
        group_by(loop.feature) %>% 
        compute_nesting_stat_rolling_correlations(window.size=window.size)
    # combine results
    bind_rows(
        per.chr.stats,
        gw.stats %>% add_column(chr='GW')
    )
}

compute_all_loop_nesting_correlation_results <- function(
    all.loop.data.df,
    all.bins.df,
    window.sizes,
    force_redo=FALSE){
    all.loop.data.df %>%
    inner_join(
        .,
        {.},
        suffix=c('.Numerator', '.Denominator'),
        by=
            join_by(
                resolution,
                IDR2D.Params,
                Comparison,
                Comparison.side,
                loop.status
            )
    ) %>%
    inner_join(
        ALL_SAMPLE_MERGED_MATRIX_COMPARISONS,
        by=colnames(ALL_SAMPLE_MERGED_MATRIX_COMPARISONS)
    ) %>% 
    left_join(
        all.bins.df,
        by=join_by(resolution)
    ) %>% 
    cross_join(tibble(window.size=window.sizes)) %>% 
    mutate(
        results_file=
            file.path(
                ALL_LOOP_NESTING_CORR_RESULTS_DIR,
                glue('resolution_{resolution}'),
                glue('window.size_{window.size}'),
                glue('IDR2D.Params_{IDR2D.Params}'),
                glue('loop.status_{loop.status}'),
                glue('{Comparison}-loop.nesting.correlation.results.tsv')
            )
    ) %>% 
        {.} -> tmp; tmp
        tmp %>% head(1) %>% 
    future_pmap(
        .l=.,
        .f=check_cached_results,
        force_redo=force_redo,
        results_fnc=compute_nesting_correlation_results,
        .progress=TRUE
    )
}
# list and load all nesting results
list_all_loop_nesting_correlation_results <- function(){
    ALL_LOOP_NESTING_CORR_RESULTS_DIR %>% 
    parse_results_filelist(
        filename.column.name='Comparison',
        suffix='-loop.nesting.correlation.results.tsv'
    )
}

load_all_loop_nesting_correlation_results <- function(){
    list_all_loop_nesting_correlation_results() %>%
    mutate(
        nesting.results=
            pmap(
                 .l=.,
                 .f=read_tsv,
                 show_col_types=FALSE
             )
    ) %>%
    unnest(nesting.results)
}

