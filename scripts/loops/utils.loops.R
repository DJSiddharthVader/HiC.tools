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
    merge_status='merged';
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
    parse_results_filelist(suffix='-dots.tsv') %>%
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
        # {.} -> tmp; tmp
        # tmp %>% select(SampleID, loops)
    unnest(loops) %>% 
    rename_with(
        ~ stri_replace_all(
            .x, 
            regex='la_exp.([a-z]+).(value|qval)',
            '$2.$1'
        )
    ) %>% 
    select(-c(filepath))
}

post_process_cooltools_dots_results <- function(results.df) {
    results.df %>%
    dplyr::rename(
        'chr'=chrom1,
        'anchor.left'=start1,
        'anchor.right'=start2
    ) %>% 
    pivot_longer(
        c(
            starts_with('value.'),
            starts_with('qval.')
        ),
        names_to='statistic',
        values_to='value'
    ) %>% 
    separate_wider_delim(
        statistic,
        delim='.',
        names=c('statistic', 'kernel')
    ) %>% 
    pivot_wider(
        names_from=statistic,
        values_from=value
    ) %>%
    mutate(
        SampleID=str_replace_all(SampleID, '.Merged.Merged', ''),
        log10.qval=-log10(qval),
        length=anchor.right - anchor.left
    ) %>% 
    dplyr::rename('enrichment'=value) %>% 
    dplyr::select(
        type,
        weight,
        resolution,
        Edit,
        Celltype,
        Genotype,
        # CloneID,
        # TechRepID,
        # ReadFilter,
        # isMerged,
        SampleID,
        chr,
        anchor.left,
        anchor.right,
        count,
        # c_label,
        # c_size,
        length,
        kernel,
        enrichment,
        log10.qval
    )
}

filter_loop_results <- function(
    results.df,
    q.thresh=LOOP_QVALUE_THRESHOLD){
    results.df %>% 
    filter(type == 'cis') %>% 
    filter(weight == 'balanced') %>% 
    filter(kernel == 'donut') %>% 
    filter(log10.qval >= -log10(q.thresh))
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
    # all loops from rep1
    reproducible.loops.P1 <- 
        tidy_IDR2D_sided_results(
            results$rep1_df,
            metric_colname
        )
    # all loops from rep2
    reproducible.loops.P2 <- 
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
    # combine all loop results from both replicates
    # loops detected in both
    bind_rows(
        reproducible.loops.P1 %>% filter(!is.na(IDR)),
        reproducible.loops.P2 %>% filter(!is.na(IDR)),
    ) %>% 
    # loops detected in exactly one of the replicates
    bind_rows(
        reproducible.loops.P1 %>% 
            filter(is.na(IDR)) %>%
            add_column(loop.type='P1.only'),
        reproducible.loops.P2 %>%
            filter(is.na(IDR)) %>%
            add_column(loop.type='P2.only'),
    ) %>%
    distinct(pick(-c('loop.type')), .keep_all=TRUE) %>% 
    # create column indicating which loops are reproduced between conditions
    mutate(
        loop.type=
            case_when(
                is.na(loop.type) & !is.na(IDR) ~ 'detected.in.both',
                loop.type == 'P1.only'         ~ loop.type,
                loop.type == 'P2.only'         ~ loop.type,
                TRUE                           ~ NA
            )
    )
}

run_IDR2D_analysis <- function(
    loops.P1,
    loops.P2,
    metric_colname,
    value_transformation,
    ambiguity_resolution_method,
    max_gap,
    ...){
    # paste0(colnames(tmp), "=tmp$", colnames(tmp), "[[row_index]]", collapse='; ')
    loops.P1 <- 
        loops.P1 %>% 
        mutate(across(c(chr.A, chr.B), as.character)) %>% 
        mutate(across(c(start.A, start.B, end.A, end.B), as.integer)) %>% 
        select(
            chr.A, start.A, end.A,
            chr.B, start.B, end.B, 
            !!sym(metric_colname)
        ) %>%
        as.data.frame()
    loops.P2 <-  
        loops.P2 %>% 
        mutate(across(c(chr.A, chr.B), as.character)) %>% 
        mutate(across(c(start.A, start.B, end.A, end.B), as.integer)) %>% 
        select(
            chr.A, start.A, end.A,
            chr.B, start.B, end.B, 
            !!sym(metric_colname)
        ) %>%
        as.data.frame()
    # run IDR2D to define replicable loops between conditions
    estimate_idr2d(
        loops.P1,
        loops.P2, 
        value_transformation=value_transformation,
        ambiguity_resolution_method=ambiguity_resolution_method,
        max_gap=ifelse(max_gap < 1, -1L, max_gap)
    ) %>% 
    # tidy up results into nice tabular format
    tidy_IDR2D_results(
        loops.P1,
        loops.P2,
        metric_colname
    )
}

run_all_IDR2D_analysis <- function(
    nested.loops.df,
    hyper.params.df,
    force.redo,
    sample.group.comparisons,
    pair_grouping_cols,
    SampleID.fields,
    sampleID_col='SampleID',
    suffixes=c('.P1', '.P2'),
    ...){
    # force.redo=parsed.args$force.redo; sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS %>% rename( 'SampleID.P1'=Sample.Group.Numerator, 'SampleID.P2'=Sample.Group.Denominator); suffixes=c('.P1', '.P2'); pair_grouping_cols=c('kernel', 'type', 'weight', 'resolution', 'chr'); sampleID_col='SampleID'; SampleID.fields=c(NA, 'Celltype', 'Genotype')
    # list + format metadata for all specified sample groups to compare
    nested.loops.df %>% 
    enumerate_pairwise_comparisons(
        sample.group.comparisons=sample.group.comparisons,
        pair_grouping_cols=pair_grouping_cols,
        sampleID_col=sampleID_col,
        suffixes=suffixes,
        SampleID.fields=SampleID.fields,
        include_merged_col=FALSE
    ) %>% 
    # Evaluate all comparisons for all combinations of specified parameters
    cross_join(hyper.params.df) %>% 
    mutate(
        max_gap=resolution * max_gap_bins,
        output_dir=
            file.path(
                LOOPS_IDR2D_DIR,
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                glue('metric_{metric_colname}'),
                glue('resolve.method_{ambiguity_resolution_method}'),
                glue('max.gap_{max_gap}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID.P1}_vs_{SampleID.P2}-IDR2D.tsv')
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
                comparison=glue('{SampleID.P1}-{SampleID.P2}')
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
        # {.} -> tmp; tmp
        # tmp %>% 
    future_pmap(
        .l=.,
        .f=check_cached_results,
        force_redo=force.redo,
        return_data=FALSE,
        results_fnc=run_IDR2D_analysis,
        .progress=TRUE
    )
}

load_IDR2D_results <- function(filepath, ...){
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
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
    ) %>% 
    extract_all_sample_pair_metadata(
        SampleID.cols=c('SampleID.Numerator', 'SampleID.Denominator'),
        SampleID.fields=c('Edit', 'Celltype', 'Genotype'),
        suffixes=c('Numerator', 'Denominator')
    )
}

load_all_IDR2D_results <- function(resolutions=NULL){
    list_all_IDR2D_results(resolutions=resolutions) %>% 
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
    select(-c(filepath))
}

post_process_IDR2D_results <- function(results.df){
    results.df %>% 
    mutate(
        max.gap.bins.int=max.gap / resolution,
        max.gap=
            fct_reorder(
                glue('{scale_numbers(max.gap, force_chr=TRUE)}'),
                max.gap.bins.int
            ),
        max.gap.bins=
            fct_reorder(
                glue('{max.gap.bins.int} bins'),
                max.gap.bins.int
            ),
        is.loop.shared=
            case_when(
                loop.type == 'P1.only' ~ 'P1.only',
                loop.type == 'P2.only' ~ 'P2.only',
                IDR <= 0.1             ~ 'IDR < 0.1',
                IDR <= 1               ~ 'Irreproducible',
                TRUE                   ~ NA
            ) %>%
            factor(
                levels=
                    c(
                        'IDR < 0.1',
                        'Irreproducible',
                        'P1.only',
                        'P2.only'
                    )
            )
    ) %>%
    select(
        -c(
            # metric, 
            # resolve.method,
            # weight,
            # kernel,
           # max.gap,
           # max.gap.bins.int,
            # max.gap.bins,
            loop.type,
            region
        )
    ) %>%
    relocate(
        type, weight, kernel,
        resolution,
        # metric, resolve.method, max.gap, max.gap.bins, max.gap.bins.int
        metric, resolve.method, max.gap,
        chr, anchor.left, anchor.right,
        diff.value, diff.rank, IDR, is.loop.shared
    )
}

filter_loop_IDR2D_results <- function(results.df){
    results.df %>% 
    filter(kernel == 'donut') %>% 
    filter(type == 'cis') %>% 
    filter(weight == 'balanced') %>% 
    filter(metric == 'log10.qval') %>% 
    filter(resolve.method == 'value') %>%
    filter(max.gap.bins.int == 5)
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
                                .cols=c(count, length, enrichment, log10.qval),
                                .fns=list('mean'=mean, 'min'=min, 'max'=max, 'median'=median),
                                .names="{.col}-{.fn}"
                            ),
                            valency=dplyr::n()
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
           SampleID.P1 %in% samples.avail,
           SampleID.P2 %in% samples.avail
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
                SampleID.P1 == SampleID,
                chr,
                between(y$start, x$anchor.left, x$anchor.right),
                between(y$end,   x$anchor.left, x$anchor.right)
            )
    ) %>% 
    inner_join(
        expression.df,
        suffix=c('.P1', '.P2'),
        by=
            join_by(
                SampleID.P2 == SampleID,
                chr,
                start, end,
                symbol, EnsemblID
            )
    ) %>%
    # for each gene compute pvalue if mean expression is different between conditions
    add_column(
        n.rna.replicates.P1=6,
        n.rna.replicates.P2=6
    ) %>% 
    mutate(
        TPM.se.P1=TPM.sd.P1**2 / n.rna.replicates.P1,
        TPM.se.P2=TPM.sd.P2**2 / n.rna.replicates.P2,
        expr.Z=(TPM.mean.P2 - TPM.mean.P1) / sqrt(TPM.se.P1 + TPM.se.P2),
        expr.p=2 * (1 - pnorm(abs(expr.Z)))
    ) %>%
    # adjust p-values genome-wide
    group_by(
        weight, resolution, kernel,
        resolve.method, metric,
        SampleID.P1, SampleID.P2
    ) %>% 
    mutate(expr.p.adjust=p.adjust(expr.p, method='BH')) %>% 
    ungroup() %>% 
    select(
        -c(
            n.rna.replicates.P1, n.rna.replicates.P2,
            TPM.se.P1, TPM.se.P2,
            # TPM.sd.P1, TPM.sd.P2,
            # TPM.mean.P1, TPM.mean.P2,
            expr.Z, expr.p
        )
    )
}

