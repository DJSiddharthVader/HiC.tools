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
calculate_all_loop_valency <- function(
    loops.df,
    ...){
    # calculate how many loops each anchor is a part of
    loops.df %>%
    mutate(
        valency.results=
            pmap(
                .l=list(loops),
                .f=
                    function(df) {
                        df %>%
                        pivot_longer(
                            starts_with('anchor.'),
                            names_to='anchor.side',
                            names_prefix='anchor.',
                            values_to='anchor.position'
                        ) %>% 
                        group_by(anchor.position) %>% 
                        dplyr::summarize(
                            across(
                                .cols=c(count, length, enrichment, log10.qvalue),
                                .fns=list('mean'=mean, 'min'=min, 'max'=max, 'median'=median, 'var'=var),
                                .names="{.col}-{.fn}"
                            ),
                            valency=dplyr::n()
                        ) %>% 
                        pivot_longer(
                            -c(
                               anchor.position, valency
                            ),
                            names_to='tmp',
                            values_to='value',
                        ) %>%
                        separate_wider_delim(
                            tmp,
                            delim='-',
                            names=c('loop.metric', 'loop.stat')
                        ) %>%
                        pivot_wider(
                            names_from=loop.stat,
                            values_from=value
                        )
                    },
            .progress=TRUE
        )
    ) %>%
    unnest(valency.results) %>%
    select(-c(loops))
}

###################################################
# Nesting Analysis
###################################################
define_nested_loop_regions_and_compute_summary_stats <- function(
    loops.df,
    bins.df,
    resolution,
    ...){
    # paste('row.index=1', paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row.index]]", collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()', sep='; ')
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
    # pivot so each loop feature is on its own row
    pivot_longer(
        -c(seqnames, start, end, FeatureID),
        names_to='loop.feature',
        values_to='loop.value'
    ) %>%
    # for each bin + loop feature 
    group_by(
        seqnames, start, end,
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
                    'mean'=mean,
                    'min'=min,
                    'max'=max,
                    'sum'=sum
                ),
            .names="metric_{.fn}"
        )
    ) %>%
    ungroup() %>% 
    pivot_wider(
        names_from=loop.feature,
        names_glue='{.value}_{loop.feature}',
        values_from=starts_with('metric_')
    ) %>% 
    # Now collapse all sets of contiguous bins at the same nesting lvl into segments
    # so each row is a segment (start-end) instead of a single bin
    # all bins squashed into the same segment are this way overlapped by the same set of loops
    # so all summary stats per nesting segment are the same
    group_by(
        across(
            c(
                starts_with('metric_'), 
                'nesting.lvl',
                'seqnames'
            )
        )
    ) %>%
    # Take the leftmost start and righmost end for all groups of contiguous bins (i.e. segments) that
    # by definition have the same loop nesting stats, since they are covered by the same set of loops
    summarize(
        start=min(start),
        end=max(end)
    ) %>%
    ungroup() %>%
    dplyr::rename('chr'=seqnames) %>% 
    arrange(chr, start) %>% 
    relocate(chr, start, end, nesting.lvl)
}

compute_all_loop_nesting_results <- function(
    all.loops.df,
    all.bins.df){
    all.loops.df %>%
    left_join(
        all.bins.df,
        by=join_by(resolution)
    ) %>% 
    mutate(
        nesting.results=
            future_pmap(
            # pmap(
                .l=.,
                .f=define_nested_loop_regions_and_compute_summary_stats,
                .progress=TRUE
            )
    ) %>%
    unnest(nesting.results) %>% 
    select(-c(loops.df, bins.df))
    # mutate(
    #     results_file=
    #         file.path(
    #             ALL_LOOP_NESTING_RESULTS_DIR,
    #             glue('resolution_{resolution}'),
    #             glue("normalization_{normalization}"),
    #             glue("kernel_{kernel}"),
    #             glue('{SampleID}-loop.nesting.results.tsv'),
    #         )
    # ) %>% 
    # future_pmap(
    #     .l=.,
    #     .f=check_cached_results,
    #     results_fnc=define_nested_loop_regions_and_compute_summary_stats,
    #     .progress=TRUE
    # )
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
        names=c(NA, 'stat', 'loop.feature')
    )
}

###################################################
# Expression Integration Analysis
###################################################
join_expr_and_IDR2D_results <- function(
    idr2d.results.df,
    expression.df,
    ...){
    # samples.avail <- expression.df %>% colnames() %>% grep('16p.', ., value=TRUE)
    samples.avail <- unique(expression.df$SampleID)
    idr2d.results.df %>% 
    filter(
           SampleID.Numerator %in% samples.avail,
           SampleID.Denominator %in% samples.avail
    ) %>% 
    inner_join(
        expression.df,
        by=
            join_by(
                chr,
                between(y$start, x$anchor.left, x$anchor.right),
                between(y$end,   x$anchor.left, x$anchor.right)
            )
    )
}

calc_expr_loop_ztest <- function(
    idr2d.results.df,
    expression.df,
    ...){
    idr2d.results.df %>% 
    inner_join(
        expression.df,
        by=
            join_by(
                SampleID.Numerator == SampleID,
                chr,
                between(y$start, x$anchor.left, x$anchor.right),
                between(y$end,   x$anchor.left, x$anchor.right)
            )
    ) %>% 
    inner_join(
        expression.df,
        suffix=c('.Numerator', '.Denominator'),
        by=
            join_by(
                SampleID.Denominator == SampleID,
                chr,
                start, end,
                symbol, EnsemblID
            )
    ) %>%
    # for each gene compute pvalue if mean expression is different between conditions
    add_column(
        n.rna.replicates.Numerator=6,
        n.rna.replicates.Denominator=6
    ) %>% 
    mutate(
        TPM.se.Numerator=TPM.sd.Numerator**2 / n.rna.replicates.Numerator,
        TPM.se.Denominator=TPM.sd.Denominator**2 / n.rna.replicates.Denominator,
        expr.Z=(TPM.mean.Denominator - TPM.mean.Numerator) / sqrt(TPM.se.Numerator + TPM.se.Denominator),
        expr.p=2 * (1 - pnorm(abs(expr.Z)))
    ) %>%
    # adjust p-values genome-wide
    group_by(
        normalization, resolution, kernel,
        resolve.method, metric,
        SampleID.Numerator, SampleID.Denominator
    ) %>% 
    mutate(expr.p.adjust=p.adjust(expr.p, method='BH')) %>% 
    ungroup() %>% 
    select(
        -c(
            n.rna.replicates.Numerator, n.rna.replicates.Denominator,
            TPM.se.Numerator, TPM.se.Denominator,
            # TPM.sd.Numerator, TPM.sd.Denominator,
            # TPM.mean.Numerator, TPM.mean.Denominator,
            expr.Z, expr.p
        )
    )
}

