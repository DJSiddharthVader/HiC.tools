################################################################################
# Dependencies
################################################################################
# library(tidyverse)
library(stringi)
library(glue)
library(furrr)
library(TADCompare)

################################################################################
# Utilities
################################################################################
convert_boundaries_to_TADs <- function(
    bins.df,
    include.tails=FALSE,
    start.col.name='start',
    end.col.name='end',
    ...){
    # bins.df=tmp2$bins.df[[1]]
    # convert list of TAD boundaries to start/end format
    # every boundary is a start+end excpet for the first and last ones
    bins.df %>% 
    filter(is.boundary) %>%
    pull(bin.start) %>%
    {
        if (include.tails) {
            c(
                bins.df %>% pull(bin.start) %>% min(), # first bin on chr
                .,                                     # all annotated boundaries
                bins.df %>% pull(bin.start) %>% max()  # last bin on chr
            )
        } else {
            .
        }
    } %>% 
    tibble(boundaries=.) %>%
    rename(!!start.col.name := boundaries) %>% 
    mutate(!!end.col.name := lead(!!sym(start.col.name))) %>% 
    filter(!is.na(start), !is.na(end))
}

convert_TADs_to_boundaries <- function(
    # convert list of TAD boundaries to start/end format
    # every bnoundary is a start+end excpet for the first and last ones
    TAD.starts,
    TAD.ends,
    boundary.col.name='boundary',
    ...){
    c(TAD.starts, TAD.ends) %>%
    unique() %>%
    sort() %>% 
    tibble(boundary=.) %>% 
    rename(boundary.col.name=boundary)
}

compute_TAD_stats <- function(
    grouped.df,
    resolution){
    # compute summary stats for each TAD for all bins inside each TAD
    # the score is whatever score is used to decide if a bin is a TAD or not e.g. Insulation, Directionality
    grouped.df %>% 
    # must have columns start, end, bin.start, bin.end, score
    summarize(
        TAD.start.score=sum(ifelse(start == bin.start, bin.score, 0)),
        TAD.end.score=sum(ifelse(end == bin.end, bin.score, 0)),
        across(
            .cols=c(bin.score),
            .fns=
                list(
                    'min'=   partial(min,    na.rm=TRUE),
                    'q25'=   partial(stats::quantile, probs=0.25, na.rm=TRUE),
                    'mean'=  partial(mean,   na.rm=TRUE),
                    'median'=partial(median, na.rm=TRUE),
                    'q75'=   partial(stats::quantile, probs=0.75, na.rm=TRUE),
                    'max'=   partial(max,    na.rm=TRUE),
                    'total'= partial(sum,    na.rm=TRUE)
                ),
            .names="TAD.inner.{.fn}"
        )
    ) %>%
    mutate(
        TAD.length=end - start,
        TAD.bins=TAD.length / resolution
    ) %>% 
    ungroup()
}

################################################################################
# Standardize + Coallate results across methods
################################################################################
generate_all_TAD_Calling_cmds <- function(
    hyper.params.df,
    output_cmd_file,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # resolutions=parsed.args$resolutions; force_redo=FALSE; merge_status='merged';
    list_all_mcool_files(merge_status=merge_status) %>%
    dplyr::rename('mcool.filepath'=filepath) %>% 
    cross_join(hyper.params.df) %>% 
    mutate(output_dir=file.path(TAD_RESULTS_DIR, glue("TAD.method_{TAD.method}"))) %>% 
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(TAD.method, ...) {
                        case_when(
                            TAD.method == 'cooltools' ~ generate_cooltools_calling_cmd(...),
                            TAD.method == 'hiTAD'     ~ generate_hiTAD_calling_cmd(...),
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
            filter(., !file.exists(output.file))
        } else {
            .
        }
    }
}

load_all_TAD_score_results <- function(
    scores.dfs.,
    resolutions=NULL){
    scores.dfs %>% 
    bind_rows() %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    }
}

load_all_TAD_results <- function(
    TAD.dfs.,
    resolutions=NULL){
    # Bind evertything together
    TAD.results.dfs %>% 
    bind_rows() %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    }
}

post_process_all_TAD_results <- function(results.df){
    results.df %>% 
    filter(TAD.length < 10 * 1e6) %>% 
    mutate(
        TAD.size.band=
            case_when(
                length > 5e6 ~ '>  5Mb',
                length > 2e6 ~ '>  2Mb',
                length > 1e6 ~ '>  1Mb',
                length > 5e5 ~ '>  500Kb',
                TRUE         ~ '<= 500Kb'
            ) %>%
            factor(levels=c('<= 500Kb', '>  500Kb', '>  1Mb', '>  2Mb', '>  5Mb'))
    )
}

pivot_all_TADs_to_boundaries <- function(results.df){
    # Pivot so the bins that start and end at the boundary position are both included
    results.df %>% 
    dplyr::select(-starts_with('TAD.inner.')) %>% 
    mutate(TAD.index=row_number()) %>% 
    dplyr::rename(
        "start.boundary"=start,
        "end.boundary"=end,
        "start.score"=TAD.start.score,
        "end.score"=TAD.end.score
    ) %>% 
    pivot_longer(
        c(start.boundary, end.boundary, start.score, end.score),
        names_to='stat',
        values_to='value'
    ) %>%
    separate_wider_delim(
        stat,
        delim='.',
        names=c('boundary.side', 'value.type')
    ) %>% 
    pivot_wider(
        names_from=value.type,
        values_from=value
    ) %>% 
    dplyr::rename(
        "boundary.start"=boundary,
        "boundary.score"=score
    )
}

################################################################################
# Process cooltools TADs + Boundaries
################################################################################
generate_cooltools_calling_cmd <- function(
    threads,
    mfvp,
    threshold,
    window.sizes.str,
    normalization,
    resolution,
    SampleID,
    mcool.filepath,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("resolution_{resolution}"),
            glue("normalization_{normalization}"),
            glue("threshold_{threshold}"),
            glue("mfvp_{mfvp}")
        )
    # Create filenames
    mcool.uri    <- glue("{mcool.filepath}::resolutions/{resolution}")
    TAD.filename <- glue("{output_dir}/{SampleID}-TADs.tsv")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    TAD.cmd   <- glue("cooltools insulation {weight_flag} --verbose --nproc {threads} --append-raw-scores --window-pixels --min-frac-valid-pixels {mfvp} --threshold {threshold} --output {TAD.filename} {mcool.uri} {window.sizes.str}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.file=TAD.filename,
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    TAD.cmd
                ),
                collapse='; '
            )
    )
}

list_all_cooltools_insulation_results_files <- function(){
    COOLTOOLS_TAD_RESULTS_DIR %>% 
    parse_results_filelist(suffix='-TAD.tsv') %>%
    add_column(TAD.method='cooltools') %>% 
    convert_MatrixID_to_SampleID_and_SampleGroup()
}

load_cooltools_insulation_TADs <- function(
    filepath,
    resolution,
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; ')
    Insulation <- 
        filepath %>% 
        read_tsv(
            progress=FALSE,
            show_col_types=FALSE
        ) %>%
        pivot_longer(
            ends_with('000'),
            names_to='TAD.stat',
            values_to='value'
        ) %>%
        # Need to reverse for separate_wider_delim to work since some stat names have _ in them
        mutate(TAD.stat=stri_reverse(TAD.stat)) %>% 
        separate_wider_delim(
            TAD.stat,
            delim='_',
            names=c(
                'window.size',
                'stat'
            ),
            too_many='merge'
        ) %>%
        mutate(
            across(
                c(stat, window.size),
                stri_reverse
            )
        ) %>% 
        pivot_wider(
            names_from=stat,
            values_from=value
        ) %>%
        mutate(is.boundary=as.logical(is_boundary)) %>% 
        filter(!is_bad_bin) %>% 
        select(-c(is_bad_bin, region, is_boundary)) %>% 
        rename(
            'chr'=chrom,
            'bin.start'=start,
            'bin.end'=end
        )
    # Get chr boundaries to map first/last TADs
    TADs <- 
        Insulation %>%
        # select(chr, window.size, bin.start, is.boundary) %>% 
        nest(bins.df=-c(chr, window.size)) %>% 
        mutate(
            TADs=
                pmap(
                    .l=.,
                    .f=convert_boundaries_to_TADs,
                    boundary.indicator.col='is.boundary',
                    .progress=FALSE
                )
        ) %>% 
        unnest(TADs) %>%
        select(window.size, chr, start, end)
    # Join bin-wise scores to TAD intervals, compute metric summary stats over all bins per TAD 
    Insulation %>% 
    select(-c(is.boundary)) %>% 
    # pivot bin-wise stats to tidy format for computing summary stats
    pivot_longer(
        c(
            log2_insulation_score,
            n_valid_pixels,
            sum_counts,
            sum_balanced,
            boundary_strength
        ),
        names_to='TAD.metric',
        values_to='bin.score',
    ) %>% 
    filter(
        TAD.metric %in% c(
            # 'sum_balanced',
            # 'boundary_strength',
            'log2_insulation_score'
        )
    ) %>% 
    # Map all bins to which TAD they are inside of
    right_join(
        TADs,
        suffix=c('.Insulation','.TAD'),
        by=
            join_by(
                chr,
                window.size,
                within(x$bin.start, x$bin.end, y$start, y$end)
            )
    ) %>% 
    filter(!is.na(start), !is.na(end)) %>% 
    # for each TAD compute summary stats over the bin-wise scores
    group_by(
        window.size, TAD.metric,
        chr, start, end
    ) %>% 
    relocate(window.size, TAD.metric, chr, start, end, bin.start, bin.end, bin.score) %>% 
    compute_TAD_stats(resolution=resolution)
}

load_all_cooltools_insulation_TADs <- function(){
    list_all_cooltools_insulation_results_files() %>% 
        # {.} -> tmp
    mutate(
        insulation=
            # pmap(
            future_pmap(
                .,
                load_cooltools_insulation_TADs,
                .progress=TRUE
            )
    ) %>%
    unnest(insulation) %>% 
    mutate(
        window.size.bins=window.size / resolution,
        # window.size=scale_numbers(window.size, force_numeric=TRUE),
    ) %>% 
    unite(
        'TAD.params',
        window.size.bins, mfvp, threshold,
        sep='#',
        remove=TRUE
    ) %>% 
    select(
        -c(
            window.size,
            filepath
        )
    )
}

load_cooltools_insulation_scores <- function(
    filepath,
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; ')
    filepath %>% 
    read_tsv(
        progress=FALSE,
        show_col_types=FALSE
    ) %>%
    pivot_longer(
        ends_with('000'),
        names_to='TAD.stat',
        values_to='value'
    ) %>%
    # Need to reverse for separate_wider_delim to work since some stat names have _ in them
    mutate(TAD.stat=stri_reverse(TAD.stat)) %>% 
    separate_wider_delim(
        TAD.stat,
        delim='_',
        names=c(
            'window.size',
            'stat'
        ),
        too_many='merge'
    ) %>%
    mutate(
        across(
            c(stat, window.size),
            stri_reverse
        )
    ) %>% 
    pivot_wider(
        names_from=stat,
        values_from=value
    ) %>%
    mutate(is.boundary=as.logical(is_boundary)) %>% 
    filter(!is_bad_bin) %>% 
    select(-c(is_bad_bin, region, is_boundary)) %>% 
    rename(
        'chr'=chrom,
        'bin.start'=start,
        'bin.end'=end
    ) %>% 
    # pivot bin-wise stats to tidy format for computing summary stats
    pivot_longer(
        c(
            log2_insulation_score,
            n_valid_pixels,
            sum_counts,
            sum_balanced,
            boundary_strength
        ),
        names_to='TAD.metric',
        values_to='bin.score',
    ) %>% 
    filter(
        TAD.metric %in% c(
            'sum_balanced',
            'boundary_strength',
            'log2_insulation_score'
        )
    )
}

load_all_cooltools_insulation_scores <- function(){
    list_all_cooltools_insulation_results_files() %>% 
        head(5) %>% 
        # {.} -> tmp
    mutate(
        insulation=
            # pmap(
            future_pmap(
                .,
                load_cooltools_insulation_scores,
                .progress=TRUE
            )
    ) %>%
    unnest(insulation) %>% 
    select(-c(filepath))
}

################################################################################
# Generate + Process hiTAD TADs + Boundaries
################################################################################
generate_hiTAD_calling_cmd <- function(
    threads,
    normalization,
    resolution,
    SampleID,
    mcool.filepath,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("resolution_{resolution}"),
            glue("normalization_{normalization}")
        )
    # Create filenames
    input.filename <- glue("{mcool.filepath}::resolutions/{resolution}")
    log.filename   <- glue("{output_dir}/{SampleID}-hitad.log")
    DI.filename    <- glue("{output_dir}/{SampleID}-DI.tsv")
    TAD.filename   <- glue("{output_dir}/{SampleID}-TADs.tsv")
    # Compose command to generate TAD for this set of inputs + params
    weight <- 
        case_when(
            normalization == 'balanced' ~ 'weight',
            normalization == 'raw'      ~ 'RAW',
            .unmatched="error"
        )
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    TAD.cmd   <- glue("domaincaller --cpu-core {threads} --weight-col {weight} --logFile {log.filename} --DI-output {DI.filename} --output {TAD.filename} --uri {input.filename}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.file=TAD.filename,
        cmd=
            paste0(
                c(
                    mkdir.cmd,
                    TAD.cmd
                ),
                collapse='; '
            )
    )
}

list_all_hiTAD_results_files <- function(){
    HITAD_TAD_RESULTS_DIR %>% 
    parse_results_filelist(suffix='.tsv') %>%
    # hiTAD only expects balanced matrices as input 
    add_column(TAD.method='hiTAD') %>% 
    separate_wider_delim(
        filename,
        delim='-',
        names=c('MatrixID', 'feature.type')
    ) %>% 
    pivot_wider(
        names_from=feature.type,
        names_glue="{feature.type}.filepath",
        values_from=filepath
    ) %>% 
    convert_MatrixID_to_SampleID_and_SampleGroup()
}

load_hiTAD_TADs <- function(
    TAD.filepath,
    DI.filepath,
    resolution,
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; '); row.index=1
    # Load bin-wise Adaptibe Directinality scores
    DIs <- 
        read_tsv(
            DI.filepath,
            show_col_types=FALSE,
            progress=FALSE,
            col_names=
                c(
                    'chr',
                    'bin.start',
                    'bin.end',
                    'bin.score'
                )
        )
    # load actual TAboundaries
    TADs <- 
        read_tsv(
            TAD.filepath,
            show_col_types=FALSE,
            progress=FALSE,
            col_names=
                c(
                    'chr',
                    'start',
                    'end'
                )
        )
    # Annotatte ADI summary stats for each TAD
    TADs %>% 
    # Map scores to TADs
    left_join(
        DIs,
        suffix=c('.TAD', '.DI'),
        by=
            join_by(
                chr,
                within(y$bin.start, y$bin.end, x$start, x$end)
            )
    ) %>% 
    add_column(TAD.metric='ADI') %>% 
    # for each TAD compute summary stats over the bin-wise scores
    group_by(TAD.metric, chr, start, end) %>% 
    compute_TAD_stats(resolution=resolution)
}

load_all_hiTAD_TADs <- function(){
    list_all_hiTAD_results_files() %>% 
    mutate(
        TADs=
            # pmap(
            future_pmap(
                .,
                load_hiTAD_TADs,
                .progress=TRUE
            )
    ) %>%
    unnest(TADs) %>% 
    add_column(TAD.params=NA) %>% 
    select(-c(ends_with('.filepath')))
}

load_hiTAD_DIs <- function(
    TAD.filepath,
    DI.filepath,
    resolution,
    ...) {
    DIs <- 
        read_tsv(
            DI.filepath,
            show_col_types=FALSE,
            progress=FALSE,
            col_names=
                c(
                    'chr',
                    'bin.start',
                    'bin.end',
                    'bin.score'
                )
        )
    TADs <- 
        read_tsv(
            TAD.filepath,
            show_col_types=FALSE,
            progress=FALSE,
            col_names=
                c(
                    'chr',
                    'TAD.start',
                    'TAD.end'
                )
        )
    # Annotatte ADI summary stats for each TAD
    # Map scores to TADs
    DIs %>% 
    left_join(
        TADs,
        by=
            join_by(
                chr,
                within(x$bin.start, x$bin.end, y$TAD.start, y$TAD.end)
            )
    ) %>% 
    mutate(
        boundary.type=
            case_when(
                TAD.start == bin.start              ~ 'TAD Start',
                TAD.end == bin.end                  ~ 'TAD End',
                !is.na(TAD.start) & !is.na(TAD.end) ~ 'TAD Interior',
                TRUE                                ~ 'Not inside TAD'
            ),
        is.boundary=boundary.type %in% c('TAD Start', 'TAD End'),
    ) %>% 
    add_column(TAD.metric='ADI')
}

load_all_hiTAD_DIs <- function(){
    list_all_hiTAD_results_files() %>% 
    mutate(
        scores=
            # pmap(
            future_pmap(
                .,
                load_hiTAD_DIs,
                .progress=TRUE
            )
    ) %>%
    unnest(scores) %>% 
    add_column(TAD.params=NA) %>% 
    select(-c(ends_with('.filepath')))
}

################################################################################
# Generate + Process ConsensusTAD TADs
################################################################################
TADCompare_load_matrix <- function(
    filepath,
    ...){
    load_mcool_file(
        filepath,
        type='df',
        cis=TRUE,
        ...
    ) %>% 
    select(c(range1, range2, IF))
}

run_ConsensusTAD <- function(
    samples.df,
    resolution,
    normalization,
    range1,
    range2,
    z_thresh,
    window_size,
    gap_thresh,
    ...){
    # paste('row.index=1;', paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; '))
    sampleID.mapping <- 
        samples.df %>%
        mutate(
            og.sample=as.character(glue('Sample {row_number()}')),
            SampleID=as.character(glue('{SampleID}.score'))
        ) %>%
        select(SampleID, og.sample) %>% 
        deframe()
    # Generate ConsensusTAD results
    consensus.results <- 
        samples.df %>%
        # named list of all Hi-C replicate matrices to use
        mutate(
            contact.matrix=
                pmap(
                    .l=.,
                    .f=TADCompare_load_matrix,
                    resolution=resolution,
                    normalization=normalization,
                    range1=range1,
                    range2=range2,
                    .progress=TRUE
                    # .progress=FALSE
                )
        ) %>%
        select(SampleID, contact.matrix) %>%
        deframe() %>% 
        # pass named list of input matrics to ConsensusTAD
        ConsensusTADs(
            resolution=resolution,
            z_thresh=z_thresh,
            window_size=window_size,
            gap_thresh=gap_thresh
        )
    # Save data for all bins + consensus boundaries
    consensus.results$All_Regions %>% 
    tibble() %>%
    # Specify which bins were actually called as boundaries
    left_join(
        consensus.results$Consensus %>% 
        tibble() %>% 
        select(Coordinate) %>% 
        add_column(isConsensusBoundary=TRUE),
        by=join_by(Coordinate)
    ) %>% 
    # rename score columns to match SampleIDs
    mutate(isConsensusBoundary=ifelse(is.na(isConsensusBoundary), FALSE, isConsensusBoundary)) %>% 
    dplyr::rename(
        all_of(sampleID.mapping),
        'bin.start'=Coordinate,
        'consensus.score'=Consensus_Score
    )
}

run_all_ConsensusTADs <- function(
    sample.groups.df,
    hyper.params.df,
    chromosomes=CHROMOSOMES,
    force_redo=TRUE,
    ...){
    # chromosomes=CHROMOSOMES; force_redo=TRUE;
    sample.groups.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    cross_join(tibble(chr=chromosomes)) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output files accordingly
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                CONSENSUSTAD_TAD_RESULTS_DIR,
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                glue('normalization_{normalization}'),
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{SampleID}-ConsensusTADs.tsv')
            )
    ) %>% 
    {
        if (!force_redo) {
            filter(., !(file.exists(results_file)))
        } else{
            .
        }
    } %T>% 
    {
        message('Generating the following results files')
        print(
            dplyr::count(
                .,
                z_thresh,
                window_size,
                gap_thresh,
                resolution,
                Sample.Group
            )
        )
    } %>%
    arrange(desc(resolution)) %>% 
        # {.} -> tmp; tmp
        # tmp %>% head(3) %>% 
    # pmap(
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_ConsensusTAD
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_ConsensusTAD,
                    # all extra args also passed as input to run_ConsensusTADs() by pmap
                    ...  # passed from the call run_all_multiHiCCompare()
                )
            },
        ...,  # passed from the call to this wrapper function by pmap()
        .progress=TRUE
    )
}

list_all_ConsensusTAD_TADs <- function(){
    CONSENSUSTAD_TAD_RESULTS_DIR %>% 
    parse_results_filelist(
        filename.column='Sample.Group',
        suffix='-ConsensusTADs.tsv'
    ) %>%
    add_column(method='ConsensusTAD')
}

load_ConsensusTAD_TADs <- function(
    filepath,
    resolution,
    ...){
    # row_index=1; paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', sep='', collapse='; ')
    # row.index=51; filepath=tmp$filepath[[row_index]]; z.thresh=tmp$z.thresh[[row_index]]; window.size=tmp$window.size[[row_index]]; gap.thresh=tmp$gap.thresh[[row_index]]; resolution=tmp$resolution[[row_index]]; region=tmp$region[[row_index]]; Sample.Group=tmp$Sample.Group[[row_index]]; method=tmp$method[[row_index]]
    # Load all the bin-wise + sample-wise scores TAD scores
    scores.df <- 
        filepath %>% 
        read_tsv(
            show_col_types=FALSE,
            progress=FALSE,
        ) %>%
        # to uniquely identify rows (boundaries) for pivoting
        mutate(
            idx=row_number(),
            across(ends_with('.score'), as.numeric)
        ) %>% 
        # has 1 column for every sample's individual TAD score + consensus score
        pivot_longer(
            ends_with('.score'),
            names_to='SampleID',
            values_to='score'
        ) %>%
        # only keep consensus score for each bin
        filter(SampleID == 'consensus.score') %>%
        # filter(!is.na(score)) %>% 
        pivot_wider(
            names_from=SampleID,
            values_from=score
        ) %>%
        dplyr::rename('bin.score'=consensus.score) %>% 
        mutate(bin.end=bin.start + resolution) %>% 
        select(-c(idx))
    # Now get list of boundaries and transform to TAD start:end pairs
    TADs.df <- 
        scores.df %>% 
        dplyr::rename('is.boundary'=isConsensusBoundary) %>% 
        select(bin.start, is.boundary) %>% 
        nest(bins.df=c(bin.start, is.boundary)) %>% 
        mutate(
            TADs=
                pmap(
                    .l=.,
                    .f=convert_boundaries_to_TADs,
                    .progress=FALSE
                )
        ) %>% 
        unnest(TADs) %>%
        select(start, end)
    # scores.df %>% nrow()
    # TADs.df %>% nrow()
    # Now map score data for each TAD together and compute summary stats
    TADs.df %>% 
    # Map scores to TADs
    left_join(
        scores.df,
        suffix=c('.TAD', '.bin'),
        by=
            join_by(
                within(y$bin.start, y$bin.end, x$start, x$end)
            )
    ) %>% 
    # for each TAD compute summary stats over the bin-wise scores
    add_column(TAD.metric='Consensus.Score') %>% 
    group_by(start, end) %>% 
    compute_TAD_stats(resolution=resolution)
}

load_all_ConsensusTAD_TADs <- function(){
    list_all_ConsensusTAD_TADs() %>% 
    mutate(
        TADs=
            # pmap(
            future_pmap(
                .,
                load_ConsensusTAD_TADs,
                .progress=TRUE
            )
    ) %>%
    unnest(TADs) %>% 
    unite(
        'TAD.params',
        sep='#',
        z.thresh,
        window.size,
        gap.thresh
    ) %>% 
    dplyr::rename('chr'=region) %>% 
    select(-c(filepath))
}

