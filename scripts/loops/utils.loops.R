library(stringi)
library(furrr)
library(idr2d)
library(broom)
library(twosamples)
library(infer)
library(lsa) # for cosine distance
library(plyranges)

######################################################################
# cooltools dots
######################################################################
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
    filter(normalization == 'balanced') %>% 
        # head(3) %>% 
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
    ) %>% 
    dplyr::rename(
        'chr'=chrom1,
        'anchor.left'=start1,
        'anchor.right'=start2,
        'enrichment'=value,
        'qvalue'=qval
    ) %>% 
    select(-c(ends_with(c('1', '2')), cstart1 ,cstart2, c_label, c_size ,region, filepath))
}

######################################################################
# IDR2D Analysis
######################################################################
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
    # paste('row.index=1;', paste0(colnames(tmp2), "=tmp2$", colnames(tmp2), "[[row.index]]", collapse='; '), ';t(head(tmp2, 1))')
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

######################################################################
# Valency Analysis
######################################################################
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

######################################################################
# Loop Nesting Results
######################################################################
compute_nesting_stats_per_bin <- function(
    loops.df,
    bins.df,
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
    dplyr::rename('chr'=seqnames) %>% 
    # for each genomic bin + loop feature 
    group_by(chr, start, end) %>%
    # calculate loop feature summary statistics + count nesting lvl (how many loops overlap each bin)
    summarize(
        nesting.lvl=n(),
        across(
            .cols=c(value),
            .fn=
                list(
                    'var'=var,
                    'min'=min,
                    'mean'=mean,
                    'max'=max,
                    'total'=sum
                ),
            .names="stat_{.fn}"
        )
    ) %>%
    ungroup()
}

squash_binwise_nesting_data_into_segments <- function(
    filepath,
    ...){
    # filepath with binwise nesting results
    filepath %>%
    read_tsv(show_col_types=FALSE) %>% 
    # Now collapse all sets of contiguous bins at the same nesting lvl into segments
    # so each row is a segment (start-end) instead of a single bin
    # all bins squashed into the same segment are this way are overlapped by 
    # the same set of loops they all have the same summary stats
    group_by(across(-c(start, end))) %>% 
    # Take the leftmost start and righmost end across all bins
    # per group as the segment start/end and just save the segment coords
    summarize(
        start=min(start),
        end=max(end)
    ) %>%
    ungroup() %>% 
    arrange(chr, start, end) %>% 
    relocate(chr, start, end)
}

squash_all_binwise_nesting_data_into_segments <- function(
    results.files.df,
    ...){
    results.files.df %>% 
    mutate(
        segments=
            future_pmap(
                 .l=.,
                 .f=squash_binwise_nesting_data_into_segments,
                 .progress=TRUE
            )
    ) %>% 
    unnest(segments) %>% 
    pivot_longer(
        starts_with('stat_'),
        names_to='stat',
        # names_prefix='stat_',
        values_to='value'
    ) %>% 
    unite(
        'feature_stat',
        sep='_',
        c(stat, feature)
    ) %>% 
    pivot_wider(
        names_from=feature_stat,
        values_from=value
    ) %>% 
    pivot_longer(
        c(nesting.lvl, starts_with('stat_')),
        names_to='feature',
        names_prefix='stat_',
        values_to='value'
    ) %>% 
    filter(!is.na(value)) %>% 
    select(-c(filepath))
}

######################################################################
# Loop Nesting Comparison Analysis 
######################################################################
load_binwise_nesting_data <- function(filepath){
    filepath %>% 
    read_tsv(show_col_types=FALSE, progress=FALSE) %>% 
    pivot_longer(
        c(starts_with('stat_'), nesting.lvl),
        names_to='stat',
        names_prefix='stat_',
        values_to='value'
    )
}

compute_nesting_summary_stats_per_segment <- function(
    segment.nesting.comparison.df,
    alternative='two.sided',
    ...){
    summary.stats <- 
        segment.nesting.comparison.df %>%
        arrange(chr, start, end) %>% 
        mutate(FC=value.Numerator / value.Denominator) %>% 
        summarize(
            # n.bins=n(),
            n.bins.larger.in.numerator=sum(FC > 1, na.rm=TRUE),
            cosine.dist=  cosine(value.Numerator, value.Denominator)[[1]],
            corr_pearson= cor(value.Numerator,    value.Denominator, method='pearson'),
            corr_kendall= cor(value.Numerator,    value.Denominator, method='kendall'),
            corr_spearman=cor(value.Numerator,    value.Denominator, method='spearman'),
            across(
                .cols=c(FC),
                .fns=
                    list(
                        'mean'= partial(mean, na.rm=TRUE),
                        'max'=  partial(max,  na.rm=TRUE),
                        'var'=  partial(var,  na.rm=TRUE),
                        'total'=partial(sum,  na.rm=TRUE)
                    ),
                .names="nest.FC_{.fn}"
            )
        )
    test.results <- 
        segment.nesting.comparison.df %>%
        arrange(chr, start, end) %>% 
        summarize(
            # test_AD=
            #     ad_test(
            #         value.Numerator,
            #         value.Denominator,
            #     ) %>% 
            #     as_tibble_row() %>% 
            #     mutate(across(everything(), as.numeric)) %>% 
            #     dplyr::rename('statistic'=`Test Stat`, 'p.value'=`P-Value`),
            test_KS=
                ks.test(
                    value.Numerator,
                    value.Denominator,
                    alternative=alternative
                ) %>% tidy(),
            test_Wilcox=
                tryCatch( 
                    {
                        wilcox.test(
                            value.Numerator,
                            value.Denominator,
                            paired=TRUE,
                            alternative=alternative
                        ) %>% 
                        tidy()
                    },
                    error=function(e) { tibble_row() }
                ),
            test_Welch=
                tryCatch( 
                    {
                        t.test(
                            value.Numerator,
                            value.Denominator,
                            paired=TRUE,
                            alternative=alternative
                        ) %>% 
                        tidy()
                    },
                    error=function(e) { tibble_row() }
                ),
            test_Sign=
                binom.test(
                    x=sum(value.Numerator > value.Denominator),
                    n=n(),
                    alternative=alternative
                ) %>% tidy(),
            # test_Kendall=
            #     cor.test(
            #         value.Numerator,
            #         value.Denominator,
            #         method='kendall',
            #         alternative=alternative
            #     ) %>% tidy(),
            # test_Spearman=
            #     cor.test(
            #         value.Numerator,
            #         value.Denominator,
            #         method='spearman',
            #         alternative=alternative
            #     ) %>% tidy(),
            test_Pearson=
                cor.test(
                    value.Numerator,
                    value.Denominator,
                    method='pearson',
                    alternative=alternative
                ) %>% tidy()
        ) %>%
        pivot_longer(
            starts_with('test_'),
            names_to='test.type',
            names_prefix='test_',
            values_to='test.results'
        ) %>% 
        unnest(test.results) %>%
        select(test.type, p.value) %>% 
        pivot_wider(names_from=test.type, names_prefix='p.value_', values_from=p.value)
    # bind everything together in tidy format
    bind_cols(
        summary.stats,
        test.results
    ) %>%
    pivot_longer(
        everything(),
        names_to='diff.stat',
        values_to='value'
    ) %>%
    add_column(n.bins=nrow(segment.nesting.comparison.df))
}

compute_nesting_correlation_results <- function(
    filepath.Numerator,
    filepath.Denominator,
    ...){
    # paste('row.index=8', paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row.index]]", collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()', sep='; ')
    # row.index=8; SampleID.Numerator=tmp$SampleID.Numerator[[row.index]]; SampleID.Denominator=tmp$SampleID.Denominator[[row.index]]; filepath.Numerator=tmp$filepath.Numerator[[row.index]]; resolution=tmp$resolution[[row.index]]; feature=tmp$feature[[row.index]]; filepath.Denominator=tmp$filepath.Denominator[[row.index]]; results_file=tmp$results_file[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% t()
    # for each bin get nesting stats for numerator + denominator paired together
    comparisons.df <- 
        # join numerator and denominator nesting stats 
        full_join(
            filepath.Numerator   %>% load_binwise_nesting_data(),
            filepath.Denominator %>% load_binwise_nesting_data(),
            suffix=c('.Numerator', '.Denominator'),
            by=join_by(chr, start, end, stat)
        ) %>%
        # set nesting stats to 0 if any bin is only overlapped by loops in either numerator or denominator
        mutate(
            across(
                c(value.Numerator, value.Denominator),
                ~ ifelse(is.na(.x), 0, .x)
            )
        ) %>% 
        # convert to granges object for computing overlaps 
        dplyr::rename('seqnames'=chr) %>% 
        as_granges()
    # squash bins into contiguous segments, since only bins ith >= 1 overlapping loop are in this 
    # list, contiguous segments must have >= 1 loop overlapping in either condition 
    # and breaks between segments must have 0 loops overlapping in both conditions
    segments.df <- 
        comparisons.df %>%
        select(-c(everything())) %>% 
        # squash all contiguous bins into segments (1 segment per row)
        # to group bins for computing summary/correlation stats between conditions per segment
        reduce_ranges() %>% 
        # clean up column names
        as_tibble() %>% 
        unite('SegmentID', sep='#', c(seqnames, start, end), remove=FALSE) %>% 
        as_granges()
    # now group bins by which contiguous segment they are within for computing summart stats
    comparisons.df %>% 
    join_overlap_inner_within(segments.df) %>% 
    as_tibble() %>%
    dplyr::rename('chr'=seqnames) %>%
    select(-c(strand, width)) %>% 
    nest(
        segment.nesting.comparison.df=
            c(
                chr, start, end,
                value.Numerator, value.Denominator
            )
    ) %>% 
    # compute different/correlation stats per segment
    mutate(
        segment.stats.df=
            future_pmap(
                 .l=.,
                 .f=compute_nesting_summary_stats_per_segment,
                 .progress=FALSE
            )
    ) %>%
    select(-c(segment.nesting.comparison.df)) %>% 
    unnest(segment.stats.df) %>% 
    separate_wider_delim(
        SegmentID,
        delim='#',
        names=c('chr', 'start', 'end'),
        cols_remove=FALSE
    )
}

load_all_segmentwise_nesting_difference_results <- function(
    results.files.df,
    ...){
    results.files.df %>% 
    mutate(
        diff.results=
            future_pmap(
                 .l=list(filepath),
                 .f=read_tsv,
                 show_col_types=FALSE,
                 progress=FALSE,
                 .progress=TRUE
            )
    ) %>% 
    unnest(diff.results) %>% 
    unite(
        'feature.stat',
        sep='_',
        c(stat, feature)
    ) %>% 
    filter(!is.na(value)) %>% 
    select(-c(filepath))
}

######################################################################
# List loop nesting resulst files
######################################################################
# Nesting analysis across loops stratified by differential status
load_nesting_results <- function(results.type){
    if (results.type == 'binwise.nesting'){
        ALL_LOOP_NESTING_RESULTS_DIR %>%
        parse_results_filelist(
            filename.column.name='SampleID',
            suffix='-binwise.nesting.stats.tsv'
        )
    } else if (results.type == 'nesting.differences'){
        ALL_LOOP_NESTING_DIFFERENCE_DIR %>% 
        parse_results_filelist(
            filename.column.name='Comparison',
            suffix='-nesting.difference.stats.tsv'
        ) %>% 
        separate_wider_delim(
            Comparison,
            delim='-',
            names=c('SampleID.Numerator', 'SampleID.Denominator')
        )
    } else if (results.type == 'binwise.nesting.by.reproducibility'){
        ALL_IDR2D_NESTING_RESULTS_DIR %>%
        parse_results_filelist(
            filename.column.name='Comparison',
            suffix='-binwise.nesting.stats.tsv'
        ) %>% 
        separate_wider_delim(
            Comparison,
            delim='-',
            names=c('SampleID.Numerator', 'SampleID.Denominator')
        )
    } else if (results.type == 'nesting.differences.by.reproducibility'){
        ALL_IDR2D_NESTING_DIFFERENCE_DIR %>% 
        parse_results_filelist(
            filename.column.name='meta.comparison',
            suffix='-nesting.difference.stats.tsv'
        ) %>% 
        separate_wider_delim(
            meta.comparison,
            delim='-',
            names=c('Comparison.Numerator', 'Comparison.Denominator')
        )
    } else {
        stop(glue('Invalid results.type: {results.type}'))
    }
}

