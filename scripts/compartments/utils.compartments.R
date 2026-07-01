######################################################################
# Depdendencies
######################################################################
library(lsa)
library(broom)

######################################################################
# Generate Cooltools Results
######################################################################
generate_cooltools_calling_cmd <- function(
    threads,
    normalization,
    resolution,
    MatrixID,
    mcool.filepath,
    track.type,
    phasing.track.filepath,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("track.type_{track.type}"),
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}")
        )
    # Create filepaths
    input.filepath     <- glue("{mcool.filepath}::resolutions/{resolution}")
    compartment.prefix <- glue("{output_dir}/{MatrixID}-")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd       <- glue("mkdir -p {output_dir}")
    compartment.cmd <- glue("cooltools eigs-cis --phasing-track {phasing.track.filepath} --n-eigs 3 {weight_flag}  -o {compartment.prefix} {input.filepath}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=glue("{compartment.prefix}.cis.vecs.tsv"),
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    compartment.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_cooltools_calling_cmds <- function(
    hyper.params.df,
    cmds.output.filepath=NULL,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # list all hyper-params and corresponding input files (i.e. phasing track files)
    # to call compartments with
    GENOME_TRACK_FILES_DIR %>%
    parse_results_filelist(suffix='-genome.track.tsv') %>%
    dplyr::rename('phasing.track.filepath'=filepath) %>% 
    inner_join(
        hyper.params.df,
        by=join_by(track.type, resolution)
    ) %>% 
    # list contacts matrices for all samples to generate compartments for
    cross_join(
        list_all_mcool_files(merge_status=merge_status) %>%
        dplyr::rename('mcool.filepath'=filepath),
    ) %>% 
    mutate(output_dir=file.path(COMPARTMENTS_RESULTS_DIR, glue("Comp.method_{Comp.method}"))) %>% 
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(Comp.method, ...) {
                        case_when(
                            Comp.method == 'cooltools' ~ generate_cooltools_calling_cmd(...),
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

######################################################################
# Load cooltools results
######################################################################
quantize_compartment_scores <- function(
    scores.df,
    n.compartment.lvls,
    ...){
    # scores.df=tmp2$scores.df[[1]]; n.compartment.lvls=tmp2$n.compartment.lvls[[1]]
    scores.df %>%
    # spearate annotation for each separate PC* being analyzed
    pivot_longer(
        -c(start, end),
        names_to='score.source',
        values_to='score'
    ) %>% 
    group_by(score.source) %>% 
    arrange(start, end) %>% 
    # remove bins with no PC* signal
    # filter(!is.na(lag(score)) | !is.na(score) | !is.na(lead(score))) %>% 
    filter(!is.na(score)) %>% 
    mutate(
        # calculate PC1 difference between each bin and the next bin,
        # more important/different switches should have larger differences?
        score.change=lead(score) - score,
        # binnify genomic bins by abs PC1 score
        strength.lvl=
            cut(
                x=abs(score),
                breaks=n.compartment.lvls,
                labels=seq(1, n.compartment.lvls)
            ) %>%
            as.integer(),
        # difference in PC1 bin between adjacent genomic bins
        strength.lvl.change=
            case_when(
                is.na(strength.lvl)       ~ NA,
                is.na(lead(strength.lvl)) ~ 0 - strength.lvl,
                is.na(lag(strength.lvl))  ~     strength.lvl,
                .default=lead(strength.lvl) - strength.lvl
            ),
        # classify each bin as A or B compartments by PC1 sign (already oriented by cooltools eigs-cis)
        compartment=ifelse(score > 0, 'A', 'B'),
        # note all bins where we shift from A->B or B->A compartments or no switch
        compartment.change=
            case_when(
                is.na(compartment)       ~ NA,
                is.na(lead(compartment)) ~ NA,
                .default=glue('{compartment}->{lead(compartment)}')
            ),
        # integrate strength + switch annotations for downstream analyses
        does.compartment.switch=compartment != lead(compartment),
        compartment.switch.type=
            case_when(
                 is.na(does.compartment.switch)                     ~ NA,
                 does.compartment.switch                            ~ 'switch',
                 !does.compartment.switch & strength.lvl.change > 0 ~ 'stronger',
                 !does.compartment.switch & strength.lvl.change < 0 ~ 'weaker',
                 .default='no.switch'
            )
    )
}

load_cooltools_compartment_results <- function(
    filepath,
    # n.compartment.lvls.list,
    n.compartment.lvls,
    ...){
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
    dplyr::rename('chr'=chrom) %>% 
    # select(chr, start, end, E1, E2, E3) %>% 
    select(chr, start, end, E1) %>% 
    nest(scores.df=-c(chr)) %>% 
    add_column(n.compartment.lvls=n.compartment.lvls) %>% 
    # cross_join(tibble(n.compartment.lvls=n.compartment.lvls.list)) %>% 
    mutate(
        compartment.labels=
            pmap(
                .l=.,
                .f=quantize_compartment_scores
            )
    ) %>% 
    select(-c(scores.df)) %>% 
    unnest(compartment.labels)
}

load_all_cooltools_compartment_results <- function(
    resolutions=NULL,
    n.compartment.lvls=20){
    COMPARTMENTS_RESULTS_DIR %>%
    parse_results_filelist(
        filename.column.name='MatrixID',
        suffix='-.cis.vecs.tsv'
    ) %>%  
    convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } %>% 
    mutate(
        compartments=
            future_pmap(
                .l=.,
                .f=load_cooltools_compartment_results,
                 n.compartment.lvls=n.compartment.lvls,
                .progress=TRUE
            )
    ) %>%
    unnest(compartments) %>% 
    unite(
        'compartment.params',
        sep='#',
        c(
            Comp.method,
            track.type,
            normalization,
            score.source
        ),
        remove=FALSE
    ) %>% 
    select(-c(filepath)) 
}

squash_bins_into_compartments <- function(compartments, ...){
    compartments %>% 
    {
        if ('chr' %in% colnames(.)) {
            dplyr::rename(., 'seqnames'=chr)
        } else {
            .
        }
    } %>% 
    as_granges() %>% 
    reduce_ranges(
        score_var=var(score),
        score_total=sum(score),
        score_mean=mean(score),
        score_median=median(score),
        score_min=min(score),
        socre_max=max(score),
        n.bins=n()
    ) %>% 
    as_tibble() %>% 
    dplyr::rename('chr'=seqnames)
}

squash_all_bins_into_compartments <- function(
    binwise.df,
    ...){
    binwise.df %>% 
    dplyr::rename('seqnames'=chr) %>% 
    nest(
        compartments=
            c(
                seqnames, start, end,
                score, score.change,
                strength.lvl, strength.lvl.change,
                # compartment,
                compartment.change, 
                does.compartment.switch, compartment.switch.type
            )
    ) %>% 
    mutate(
        compartments=
            future_pmap(
                .l=.,
                .f=squash_bins_into_compartments,
                .progress=TRUE
            )
    ) %>% 
    unnest(compartments) %>% 
    select(-c(strand, width)) %>% 
    mutate(length=end - start) %>% 
    arrange(chr, start, end, compartment)
}

get_switches_and_context_from_all_bins <- function(
    binwise.df,
    bin.context=10,
    ...){
    # map each bin to its nearest switch (up or downstream) and 
    # keep all bins +/- bin.context bins within a switch in either direction
    binwise.df %>% 
    dplyr::rename('seqnames'=chr) %>%
    nest(
        bins=
            c(
                seqnames, start, end,
                score, score.change,
                strength.lvl, strength.lvl.change,
                compartment, compartment.change, 
                does.compartment.switch, compartment.switch.type
            )
    ) %>% 
    mutate(bins=pmap(.l=list(bins), as_granges)) %>% 
        # {.} -> tmp2; tmp2
        # bins=tmp2$bins[[1]]; resolution=tmp2$resolution[[1]]
        # tmp2 %>% head(5) %>% 
    mutate(
        switch.and.context=
            future_pmap(
                bin.context=bin.context,
                .progress=TRUE,
                .l=.,
                .f=
                    function(bins, resolution, ...){
                        switches <- 
                            bins %>% 
                            filter(does.compartment.switch) %>%
                            select(!everything()) %>%
                            mutate(switch.idx=row_number())
                        bind_ranges(
                            join_nearest_upstream(
                                bins, 
                                switches,
                                distance=TRUE
                            ) %>%
                            mutate(side='upstream'),
                            join_nearest_downstream(
                                bins, 
                                switches,
                                distance=TRUE
                            ) %>% 
                            mutate(side='downstream')
                        ) %>% 
                        as_tibble() %>%
                        filter((distance / resolution) <= bin.context) %>% 
                        mutate(
                            distance=
                                case_when(
                                    distance == 0        ~ 0,
                                    side == 'upstream'   ~ 0 - (distance + 1),
                                    side == 'downstream' ~      distance + 1
                                )
                        ) %>% 
                        dplyr::rename('distance.to.nearest.switch'=distance)
                    }
            )
    ) %>%
    select(-c(bins)) %>% 
    unnest(switch.and.context)
}

######################################################################
# Generate Saddle Plot Data
######################################################################
generate_saddle_data_calculation_cmds <- function(
    normalization,
    resolution,
    contact.type,
    n.bins,
    qrange,
    track.col.name,
    expected.col.name,
    MatrixID,
    mcool.filepath,
    track.filepath,
    expected.path,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}"),
            glue("contat.type_{contact.type}"),
            glue("n.bins_{n.bins}"),
            glue("expectation.metric_{expected.col.name}"),
            glue("saddle.metric_{track.col.name}")
        )
    # Create filepaths
    mcool.uri     <- glue("{mcool.filepath}::resolutions/{resolution}")
    track.uri     <- glue("{track.filepath}::{track.col.name}")
    expected.uri  <- glue("{expected.path}::{expected.col.name}")
    output.prefix <- glue("{output_dir}/{MatrixID}-")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd       <- glue("mkdir -p {output_dir}")
    saddle.cmd <- glue("cooltools saddle --qrange {qrange} --strength {weight_flag} -t {contact.type} --n-bins {n.bins} -o {output.prefix} {mcool.uri} {track.uri} {expected.path}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=glue("{output.prefix}signals.tsv"),
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    saddle.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_saddle_data_calculation_cmds <- function(
    hyper.params.df,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # merge_status='merged'; force_redo=FALSE;
    # list all binwise eigenvector results files
    compartment.results.files.df <- 
        COMPARTMENTS_RESULTS_DIR %>%
        parse_results_filelist(
            filename.column.name='MatrixID',
            suffix='.cis.vecs.tsv'
        ) %>% 
        convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
        dplyr::rename('track.filepath'=filepath)
    # list all distance expected contact files
    distance.expectation.results.files.df <- 
        DISTANCE_EXPECTED_CONTACTS_DIR %>% 
        parse_results_filelist(
            filename.column.name='MatrixID',
            suffix='-expected.tsv'
        ) %>% 
        convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
        dplyr::rename('expected.path'=filepath)
    # Contact matrix files
    matrix.files.df <- 
        list_all_mcool_files(merge_status=merge_status) %>%
        dplyr::rename('mcool.filepath'=filepath)
    # Map all inputs together by matching param valuies
    compartment.results.files.df %>% 
    inner_join(
        matrix.files.df,
        by=join_by(isMerged, Sample.Group, SampleID)
    ) %>% 
    inner_join(
        distance.expectation.results.files.df,
        by=join_by(isMerged, Sample.Group, SampleID, normalization, resolution)
    ) %>% 
    # list contacts matrices for all samples to generate compartments for
    # Map any other hyper-params to sets of related files
    inner_join(
        hyper.params.df,
        by=
            join_by(
                normalization,
                resolution,
                contact.type
            )
    ) %>% 
    add_column(output_dir=COMPARTMENT_SADDLE_FILES_DIR) %>% 
    # build commands from relevant params + input files
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(Comp.method, ...) {
                        case_when(
                            Comp.method == 'cooltools' ~ generate_saddle_data_calculation_cmds(...),
                            .unmatched='error'
                        )
                    },
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    # Only include cmds generating outputfiles that dont exist
    {
        if (!force_redo) {
            filter(., !file.exists(output.filepath))
        } else {
            .
        }
    }
}

######################################################################
# Compare Compartment scores
######################################################################
calculate_binwise_difference_stats <- function(
    segment.paired.binwise.data ,
    alternative='two.sided',
    # seed=9,
    # reps=10000,
    ...){
    # segment.paired.binwise.data=tmp2$segment.paired.binwise.data[[5]]
        # segment.paired.binwise.data
    summary.stats <- 
        segment.paired.binwise.data %>% 
        # mutate(FC=score.Numerator / score.Denominator) %>% 
        mutate(diff=score.Numerator - score.Denominator) %>% 
        summarize(
            n.bins=n(),
            n.bins.diff.compartment=sum(!is.compartment.matched),
            cosine.dist=  cosine(score.Numerator, score.Denominator)[[1]],
            corr.pearson= cor(score.Numerator,    score.Denominator, method='pearson'),
            corr.kendall= cor(score.Numerator,    score.Denominator, method='kendall'),
            corr.spearman=cor(score.Numerator,    score.Denominator, method='spearman'),
            n.bins.larger.in.numerator=sum(diff > 1),
            across(
                # .cols=c(FC), .names="score.FC_{.fn}"
                .cols=c(diff), .names="diff_{.fn}",
                .fns=
                    list(
                        'max'=max,
                        'mean'=mean,
                        'var'=var,
                        'total'=sum
                    )
            )
        )
    test.results <- 
        segment.paired.binwise.data %>% 
        summarize(
            test_KS=
                ks.test(
                    score.Numerator,
                    score.Denominator,
                    alternative=alternative
                ) %>% tidy(),
            test_Wilcox=
                tryCatch( 
                    {
                        wilcox.test(
                            score.Numerator,
                            score.Denominator,
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
                            score.Numerator,
                            score.Denominator,
                            paired=TRUE,
                            alternative=alternative
                        ) %>% 
                        tidy()
                    },
                    error=function(e) { tibble_row() }
                ),
            test_Sign=
                binom.test(
                    x=sum(score.Numerator > score.Denominator),
                    n=n(),
                    alternative=alternative
                ) %>% tidy(),
            test_Pearson=
                tryCatch( 
                    {
                        cor.test(
                            score.Numerator,
                            score.Denominator,
                            method='pearson',
                            alternative=alternative
                        ) %>% tidy()
                    },
                    error=function(e) { tibble_row() }
                )
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
        -c(n.bins),
        names_to='feature',
        values_to='value'
    )
}

calculate_binwise_difference_stats_per_segment <- function(
    segments.Numerator,
    segments.Denominator,
    resolution,
    ...){
    # paste('row.index=1', paste0(colnames(segments.df), '=segments.df$', colnames(segments.df), '[[row.index]]', collapse='; '), sep='; ')
    # group sets of adjacent bins into contiguous segments, each segment has a gap of >= 1 bin 
    segments <- 
        bind_ranges(
            segments.Numerator,
            segments.Denominator
        ) %>% 
        filter(!is.na(score)) %>% 
        # squash all contiguous bins into segments (1 segment per row)
        # to group bins for computing summary/correlation stats between conditions per segment
        reduce_ranges() %>% 
        # clean up column names
        as_tibble() %>% 
        unite('SegmentID', sep='#', c(seqnames, start, end), remove=FALSE) %>% 
        as_granges()
    join_overlap_left(
        segments.Denominator,
        segments.Numerator,
        minoverlap=resolution,
        suffix=c('.Numerator', '.Denominator')
    ) %>%
    join_overlap_inner_within(segments) %>% 
    as_tibble() %>%
    dplyr::rename('chr'=seqnames) %>%
    mutate(is.compartment.matched=(compartment.Numerator == compartment.Denominator)) %>% 
    select(
        SegmentID,
        compartment.Denominator, is.compartment.matched, 
        chr, start, end,
        score.Numerator, score.Denominator,
    ) %>% 
    nest(
        segment.paired.binwise.data=
            c(
                is.compartment.matched, 
                chr, start, end,
                score.Numerator, score.Denominator,
            )
    ) %>% 
    # only keep compartments annotaed in Denominator, longer than 1 bin
    filter(!is.na(compartment.Denominator)) %>% 
    filter(pmap(.l=list(segment.paired.binwise.data), .f=nrow) > 1) %>% 
    # compute different/correlation stats per segment
    mutate(
        segment.comparisons=
            # future_pmap(
            pmap(
                 # .progress=TRUE,
                 .l=.,
                 .f=calculate_binwise_difference_stats,
                 .progress=FALSE
            )
    ) %>%
    select(-c(segment.paired.binwise.data)) %>% 
    unnest(segment.comparisons)
}

calculate_all_binwise_difference_stats <- function(
    segments.df,
    ...){
    segments.df %>% 
    mutate(
        segment.test.results=
            future_pmap(
                .l=.,
                .f=calculate_binwise_difference_stats_per_segment,
                .progress=TRUE
            )
    ) %>% 
    select(-c(starts_with('segments.'))) %>%
    unnest(segment.test.results) %>%
    select(-c(chr)) %>% 
    separate_wider_delim(
        SegmentID,
        delim='#',
        names=c('chr', 'start', 'end'),
        cols_remove=FALSE
    )
}

compute_switch_differences_between_conditions <- function(
    switches.Numerator,
    switches.Denominator,
    ...){
    # paste('row.index=1', paste0(colnames(switches.df), '=switches.df$', colnames(switches.df), '[[row.index]]', collapse='; '), sep='; ')
    # row.index=1; Sample.Group.Numerator=switches.df$Sample.Group.Numerator[[row.index]]; Sample.Group.Denominator=switches.df$Sample.Group.Denominator[[row.index]]; resolution=switches.df$resolution[[row.index]]; compartment.params=switches.df$compartment.params[[row.index]]; chr=switches.df$chr[[row.index]]; switches.Numerator=switches.df$switches.Numerator[[row.index]]; switches.Denominator=switches.df$switches.Denominator[[row.index]]
    # switches.df %>% head(row.index) %>% tail(1) %>% t()
    # a swithc is a bin where the compartment changes in the next adjacent bin 
    # i.e. A -> B or B -> A 
    # map nearest numerator swithc to each denominator switch
    join_nearest(
        switches.Denominator,
        switches.Numerator,
        suffix=c('___Denominator', '___Numerator'),
        distance=TRUE
    ) %>%
    as_tibble() %>% 
    mutate(distance=ifelse(distance == 0, 0, distance+1)) %>% 
    # pivot columns and compute difference in switch stats between switches 
    pivot_longer(
        ends_with(c('___Denominator', '___Numerator')),
        names_to='feature',
        values_to='value'
    ) %>% 
    separate_wider_delim(
        feature,
        delim='___',
        names=c('feature', 'side')
    ) %>%
    pivot_wider(
        names_from=side,
        names_prefix='value.',
        values_from=value
    ) %>%
    # FC as in fold-change for each feature
    mutate(diff=value.Numerator - value.Denominator) %>%
    mutate(FC=value.Numerator / value.Denominator) %>%
    dplyr::rename('chr'=seqnames) %>% 
    select(-c(width, strand, starts_with('value.')))
}

match_all_nearest_switches_between_conditions <- function(
    switches.df,
    ...){
    switches.df %>% 
    mutate(
        switch.matches=
            future_pmap(
                .l=.,
                .f=compute_switch_differences_between_conditions,
                .progress=TRUE
            )
    ) %>% 
    select(-c(starts_with('switches.'))) %>%
    unnest(switch.matches)
}

