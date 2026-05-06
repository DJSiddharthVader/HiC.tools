################################################################################
# Dependencies
################################################################################
library(furrr)
library(TADCompare)

################################################################################
# Calculate Measure of Concordance between sets of TADs
################################################################################
calculate_MoC <- function(
    TADs.P1,
    TADs.P2,
    ...){
    # TADs.P1, TADs.P2 are both tibbles with 
    # the following 3 columns: start, end, length
    # add indices to track all pairs of TADs
    # MoC normalization constant
    nTADs.P1 <- nrow(TADs.P1)
    nTADs.P2 <- nrow(TADs.P2)
    norm_const <- 1 / (sqrt(nTADs.P1 * nTADs.P2) - 1)
    # nTADs.P1; nTADs.P2; norm_const;
    # Now find all overlapping pairs of TADs
    inner_join(
        TADs.P1 %>% mutate(idx=row_number()),
        TADs.P2 %>% mutate(idx=row_number()),
        suffix=c('.P1', '.P2'),
        by=join_by(overlaps(x$start, x$end, y$start, y$end))
    ) %>% 
    # https://link.springer.com/article/10.1186/s13059-018-1596-9#Sec9
    # "Assessment of TAD calller performance"
    rowwise() %>% 
    mutate( 
        # F_ij^2 / (P_i * Q_j), 1 pair of TADs per row
        intersection=min(end.P1, end.P2) - max(start.P1, start.P2),
        moc.inner=((intersection**2) / (TAD.length.P1 * TAD.length.P2))
    ) %>%
    # group_by(idx.P1) %>% slice_max(moc.inner) %>%  ungroup() %>% 
    # group_by(idx.P2) %>% slice_max(moc.inner) %>% 
    # When multiple TADs overlap, count only the most overlapping match
    # group_by(idx.P1) %>%
    # slice_max(moc.inner)
    ungroup() %>% 
    summarize(
        n.Overlaps=n(),
        n.TADs.P1=length(unique(idx.P1)),
        n.TADs.P2=length(unique(idx.P2)),
        MoC=(sum(moc.inner) - 1) / (sqrt(n.TADs.P1 * n.TADs.P2) - 1)
        # MoC=(sum(moc.inner) - 1) * norm_const
    )
}

calculate_all_MoCs <- function(
    nested.TADs.df,
    suffixes=NULL,
    delim='.',
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', sep='', collapse='; ')
    # suffixes=c('Numerator', 'Denominator'); delim='.'
    nested.TADs.df %>% 
    enumerate_pairwise_comparisons(
        delim=delim,
        suffixes=c('P1', 'P2'),
        ...
    ) %>% 
    # Finally compute all MoCs for all listed pairs of annotations
    mutate(
        MoCs=
            # pmap(
            future_pmap(
                .l=.,
                .f=calculate_MoC,
                .progress=TRUE
            )
    ) %>%
    select(-c(TADs.P1, TADs.P2)) %>%
    unnest(MoCs) %>%
    {
        if (!is.null(suffixes) & length(suffixes) == 2) {
            rename_with(
                ., 
                .cols=ends_with('P1'),
                ~str_replace(.x, 'P1$', suffixes[[1]])
            ) %>% 
            rename_with(
                .cols=ends_with('P2'),
                ~str_replace(.x, 'P2$', suffixes[[2]])
            )
        } else {
            .
        }
    }
}

################################################################################
# Generate TADCompare results
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

load_all_TAD_results_for_TADCompare <- function(
    force.redo=FALSE,
    force_redo_sub=FALSE){
    # hiTAD TAD results
    all.TADs.df <- 
        check_cached_results(
            results_file=ALL_TAD_RESULTS_FILE,
            force_redo=force.redo, force_redo_sub=force_redo_sub,
            results_fnc=load_all_TAD_results
        ) %>% 
        select(
            resolution,
            Sample.Group,
            method, TAD.params,
            chr, start, end, TAD.length
        ) %>% 
        mutate(chr.copy=chr) %>% 
        # mutate(length=end - start) %>% 
        nest(TADs=c(chr, start, end, TAD.length)) %>% 
        dplyr::rename(
            'chr'=chr.copy,
            'TAD.method'=method
        )
    # default TADCompare method estimates TADs itself, include nothing
    spectralTAD.TADs.df <- 
        expand_grid(
            Sample.Group=unique(all.TADs.df$Sample.Group),
            chr=CHROMOSOMES,
            resolution=unique(all.TADs.df$resolution)
        ) %>% 
        add_column(
            TADs=NULL, # will be estimated by TADCompare
            TAD.params=NULL,
            TAD.method='spectralTAD'
        )
    # Bind everything together
    bind_rows(
        all.TADs.df,
        spectralTAD.TADs.df
    ) %>%
    unite(
        'TAD.set.index',
        sep='~',
        remove=FALSE,
        c(
          TAD.method,
          TAD.params,
          resolution
        )
    ) %>% 
    dplyr::select(-c(TAD.set.index)) %>% 
    dplyr::rename('pre_tads'=TADs)
}

run_TADCompare <- function(
    filepath.Numerator,
    Sample.Group.Numerator,
    pre_tads.Numerator,
    filepath.Denominator,
    Sample.Group.Denominator,
    pre_tads.Denominator,
    resolution,
    normalization,
    range1,
    range2,
    z_thresh,
    window_size,
    gap_thresh,
    ...){
    # paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', collapse='; ')
    # chr1 @ 10Kb -> 24896x24896 matrix -> 40Gb is enough
    # Run TADCompare on the 2 matrices being compared
    matrix.numerator <-
        TADCompare_load_matrix(
            filepath.Numerator,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2
        )
    matrix.denominator <-
        TADCompare_load_matrix(
            filepath.Denominator,
            resolution=resolution,
            normalization=normalization,
            range1=range1,
            range2=range2
        )
    pre_tads <- 
        if (is.null(pre_tads.Numerator)) {
            NULL
        } else {
            list(pre_tads.Numerator, pre_tads.Denominator)
        }
    tad.compare.results <- 
        TADCompare(
            matrix.numerator,
            matrix.denominator,
            resolution=resolution,
            z_thresh=z_thresh,
            window_size=window_size,
            gap_thresh=gap_thresh,
            pre_tads=pre_tads
        )
    # Format results to include boundary+gap scores for all bins + differential annotations
    # tad.compare.results$Boundary_Scores %>% as_tibble()
    # tad.compare.results$TAD_Frame %>% as_tibble()
    tad.compare.results$Boundary_Scores %>% 
    as_tibble() %>%
    full_join(
        tad.compare.results$TAD_Frame %>%
        as_tibble() %>% 
        add_column(isTADBoundary=TRUE),
        suffix=c('.All', '.TADs'),
        by=join_by(Boundary)
    ) %>% 
        # {.} -> tcr; tcr
        # tcr %>% count(isTADBoundary, Differential.All, Differential.TADs, Type.All, Type.TADs)
        # tcr %>% 
    mutate(
        isTADBoundary=ifelse(is.na(isTADBoundary), FALSE, isTADBoundary),
        Differential=
            case_when(
                is.na(Differential.TADs) ~ Differential.All,
                TRUE                     ~ Differential.TADs
            ),
        is.Differential=!grepl('Non-Differential', Differential),
        Type=
            case_when(
                is.na(Type.TADs) ~ Type.All,
                TRUE             ~ Type.TADs
            ),
        Enriched.Condition=
            case_when(
                Enriched_In.TADs == 'Matrix 1' ~ Sample.Group.Numerator,
                Enriched_In.TADs == 'Matrix 2' ~ Sample.Group.Denominator,
                Enriched_In.All  == 'Matrix 1' ~ Sample.Group.Numerator,
                Enriched_In.All  == 'Matrix 2' ~ Sample.Group.Denominator,
                TRUE                      ~ NA
            ),
        TAD_Score1=
            case_when(
                is.na(TAD_Score1.TADs) ~ TAD_Score1.All,
                TRUE                   ~ TAD_Score1.TADs
            ),
        TAD_Score2=
            case_when(
                is.na(TAD_Score2.TADs) ~ TAD_Score2.All,
                TRUE                   ~ TAD_Score2.TADs
            ),
        Gap_Score=
            case_when(
                is.na(Gap_Score.TADs) ~ Gap_Score.All,
                TRUE                  ~ Gap_Score.TADs
            )
    ) %>%
    dplyr::rename(
        'TAD.Score.Numerator'=TAD_Score1,
        'TAD.Score.Denominator'=TAD_Score2,
        'TAD.isDifferential'=Differential,
        'TAD.Difference.Type'=Type
    ) %>% 
    rename_with(~ str_replace_all(.x, '_', '.')) %>% 
    dplyr::select(-c(ends_with('.All'), ends_with('.TADs')))
}

run_all_TADCompare <- function(
    comparisons.df,
    hyper.params.df,
    force_redo=FALSE,
    ...){
    # force_redo=TRUE;
    comparisons.df %>% 
    # for each comparison list all paramter combinations
    cross_join(hyper.params.df) %>% 
    # Create nested directory structure listing all relevant analysis parameters
    # Name output file as {numerator}_vs_{denominator}-*.tsv
    mutate(
        range1=chr, range2=chr,
        output_dir=
            file.path(
                TADCOMPARE_DIR,
                glue('z.thresh_{z_thresh}'),
                glue('window.size_{window_size}'),
                glue('gap.thresh_{gap_thresh}'),
                glue('TAD.method_{TAD.method}'),
                glue('TAD.params_{TAD.params}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}'),
                # glue('resolution.type_{resolution.type}'),
                glue('region_{chr}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}-TADCompare.tsv')
            )
    ) %>% 
    # filter(chr != 'chrY') %>% 
    arrange(desc(resolution), desc(chr)) %>% 
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
                # z_thresh,
                # window_size,
                # gap_thresh,
                # TAD.params,
                TAD.method, 
                resolution,
                Sample.Group.Numerator,
                Sample.Group.Denominator
            )
        )
    } %>%
        # {.} -> tmp; tmp
        # tmp %>% 
    # pmap(
    future_pmap(
        .l=.,
        .f= # Need this wrapper to pass ... arguments to run_multiHiCCompare
            function(results_file, ...){ 
                check_cached_results(
                    results_file=results_file,
                    force_redo=force_redo,
                    return_data=FALSE,
                    results_fnc=run_TADCompare,
                    # all args also passed as input arguments to run_all*() by pmap
                    ...  # passed from the call to this wrapper()
                )
            },
        ...,  # passed from the call to from run_all_TADCompare
        .progress=TRUE
    )
}

################################################################################
# Load TADCompare results
################################################################################
load_TADCompare_results <- function(
    filepath,
    boundaries_only=TRUE,
    ...){
    # row_index=1; filepath=tmp$filepaths[[row_index]][[1]];
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
        # {.} -> ltr.tmp; ltr.tmp
        # ltr.tmp %>% count(isTADBoundary, is.Differential, TAD.Difference.Type)
    {
        if (boundaries_only) {
            filter(., isTADBoundary)
        } else {
            .
        }
    } %>% 
    filter(!is.na(TAD.Difference.Type))
}

load_and_correct_TADCompare_results <- function(
    filepaths,
    nom.threshold,
    # fdr.threshold,
    gw.fdr.threshold,
    ...){
    # filepaths=tmp$filepaths[[1]]
    filepaths %>% 
    mutate(
        results=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                id='tmpID',
                show_col_types=FALSE,
                progress=FALSE
            )
    ) %>% 
    unnest(results) %>% 
    # The score provided by TADCompare is functionallt a z-score distributed at N(0,1)
    # so we can compute a regular p-value to decide if a TAD's boundary score is significantly different
    # between conditions
    # See Section 2.7 here
    # https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00158/full 
    # calculate pvalue from Gap Score (a z-score) calcualted by TADCompare
    mutate(p.value=2 * pnorm(abs(Gap.Score), lower.tail=FALSE)) %>% 
    mutate(p.adj.gw=p.adjust(p.value, method='BH')) %>% 
    ungroup() %>% 
    # only keep sufficiently sifnificant siginicant differences
    filter(
        p.adj.gw < gw.fdr.threshold,
        # p.adj    < fdr.threshold,
        p.value  < nom.threshold
    ) %>%
    select(-c(tmpID))
}

list_all_TADCompare_results <- function(){
    # Get a list of all results files
    TADCOMPARE_DIR %>% 
    parse_results_filelist(
        suffix='-TADCompare.tsv',
        filename.column.name='pair.name'
    ) %>% 
    # Split title into pair of groups ordered by numerator/denominator
    separate_wider_delim(
        pair.name,
        delim='_vs_',
        names=c('SampleID.Numerator', 'SampleID.Denominator')
    ) %>% 
        # {.} -> tmp; tmp
        # tmp %>% 
    extract_all_sample_pair_metadata(
        SampleID.cols=c('SampleID.Numerator', 'SampleID.Denominator'),
        SampleID.fields=c('Edit', 'Celltype', 'Genotype'),
        suffixes=c('Numerator', 'Denominator')
    )
}

load_all_TADCompare_results <- function(
    nom.threshold,
    # fdr.threshold,
    gw.fdr.threshold,
    ...){
    # gw.fdr.threshold=1; fdr.threshold=0.1; nom.threshold=0.05
    list_all_TADCompare_results() %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    nest(filepaths=c(filepath, region)) %>% 
    mutate(
        results=
            # pmap(
            future_pmap(
                .l=.,
                # load_TADCompare_results,
                .f=load_and_correct_TADCompare_results,
                nom.threshold=nom.threshold,
                # fdr.threshold=fdr.threshold,
                gw.fdr.threshold=gw.fdr.threshold,
                boundaries_only=TRUE,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    dplyr::rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        'isDifferential'=is.Differential,
        'DifferenceType'=TAD.Difference.Type
    ) %>% 
    select(-c(filepaths))
}

load_correct_count_TADCompare_results <- function(
    filepaths,
    sig.colname='p.adj.gw',
    ...){
    filepaths %>% 
    load_and_correct_TADCompare_results(
        nom.threshold=1,
        # fdr.threshold=1,
        gw.fdr.threshold=1
    ) %>% 
    # for each thresh, make binary col if TAD difference meets threshold
    mutate(
        "sig.lvl.{sig.colname} < 1e-15" := .data[[sig.colname]] <  1e-15,
        "sig.lvl.{sig.colname} < 1e-10" := .data[[sig.colname]] <  1e-10,
        "sig.lvl.{sig.colname} < 1e-05" := .data[[sig.colname]] <  1e-5,
        "sig.lvl.{sig.colname} < 0.001" := .data[[sig.colname]] <  1e-3,
        "sig.lvl.{sig.colname} < 0.05 " := .data[[sig.colname]] <  0.05,
        "sig.lvl.{sig.colname} < 0.1  " := .data[[sig.colname]] <  0.10,
        "sig.lvl.N.S."                  := .data[[sig.colname]] >= 0.10
        # "sig.lvl.NA"                    := is.na(.data[[sig.colname]])
    ) %>% 
    pivot_longer(
        starts_with('sig.lvl.'),
        names_to='sig.lvl',
        names_prefix='sig.lvl.',
        values_to='meet.sig.lvl'
    ) %>% 
    # Inclusively count how many TAD differences meet each thrshold across categories
    # This produces inclusive counts  for each significance threshold i.e. 
    # the number of TAD differences < 0.1 also includes all differences <= 0.01
    filter(meet.sig.lvl) %>% 
    count(
        isTADBoundary,
        is.Differential,
        TAD.Difference.Type,
        Enriched.Condition,
        region,
        sig.lvl
    )
}

load_correct_count_all_TADCompare_results <- function(){
    list_all_TADCompare_results() %>% 
    # Load all results + correct pvalues genome wide per Sample.Group
    nest(filepaths=c(filepath, region)) %>% 
    mutate(
        results=
            future_pmap(
                .l=.,
                .f=load_correct_count_TADCompare_results,
                .progress=TRUE
            )
    ) %>%
    unnest(results) %>% 
    dplyr::rename(
        'chr'=region,
        'isBoundary'=isTADBoundary, 
        'isDifferential'=is.Differential,
        'DifferenceType'=TAD.Difference.Type
    ) %>% 
    select(-c(filepaths))
}

post_process_TADCompare_results <- function(results.df){
    results.df %>%
    # filter(TAD.method != 'cooltools') %>% 
    filter(isBoundary) %>% 
    mutate(
        # isBoundary=ifelse(isBoundary, 'TAD', 'Not TAD'),
        across(
            c(
                SampleID.Numerator,
                SampleID.Denominator,
                Enriched.Condition
            ),
            ~ str_remove(.x, '.Merged.Merged')
        ),
    ) %>% 
    # mutate(log.p.adj.gw=-log10(p.adj.gw)) %>% 
    select(
        -c(
            isDifferential,
            isBoundary,
            z.thresh,
            window.size,
            gap.thresh
        )
    ) %>% 
    relocate(
        c(
            resolution,
            TAD.method,
            TAD.params,
            SampleID.Numerator, SampleID.Denominator, 
            chr,
            DifferenceType,
            Enriched.Condition
        )
    )
}

