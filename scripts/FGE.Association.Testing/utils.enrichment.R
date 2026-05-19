################################################################################
# Dependencis
################################################################################
library(stringi)
library(furrr)
# library(idr2d)
# library(plyranges)

################################################################################
# Handle raw FGE signal data
################################################################################
# JASPAR 2022 CTCF data
load_CTCF_sites <- function(force.redo=FALSE){
    # https://dozmorovlab.github.io/CTCF/
    # annotation object
    check_cached_results(
        results_file=CTCF_SITE_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                AnnotationHub() %>% 
                # subset specific set of CTCF sites predicted in humans by JASPAR 2022
                subset(
                    species == 'Homo sapiens' &
                    genome == 'hg38' &
                    dataprovider == 'JASPAR 2022' &
                    preparerclass == "CTCF"
                ) %>%
                {.[['AH104727']]} %>%
                # make tidy tibble
                as_tibble() %>% 
                dplyr::rename(
                    'motif'=name,
                    'length'=width,
                    'chr'=seqnames
                ) %>% 
                select(
                    -c(
                        # db,
                        # db.set,
                        # bioset,
                        # genomic.feature,
                        sequence
                    )
                )
            }
    )
}
# ENCODE cCRE data
load_ENCODE_cCREs <- function(){
    # Graphical + explicit definitions of of cCRE types
    # https://screen-v4.wenglab.org/about
    # cCRE.descriptions <- 
    #     tribble(
    #         ~cCRE.type,   ~cCRE.description,
    #         'CA',         'Chromatin Accessible Only',
    #         'CA-CTCF',    'Chromatin Accessible with CTCF',
    #         'CA-H3K4me3', 'Chromatin Accessible with H4K4me3',
    #         'CA-TF',      'Chromatin Accessible with TF',
    #         'PLS',        'Promoter-like',
    #         'TF',         'TF Only',
    #         'dELS',       'Distal enhancer-like'
    #         'pELS',       'Proximal enhancer-like'
    #     )
    ENCODE_CCRE_SITES_FILE %>% 
    read_tsv(
        show_col_types=FALSE,
        col_names=
            c(
                'chr',
                'start', 'end',
                # 'cCREID.P1', 'cCREID.P2',
                'cCREID', 'cCREID.2',
                'cCRE.type'
            )
    ) %>% 
    select(-c(cCREID.2))
}
# TF binding sites for various TFs
load_TF_binding_sites <- function(force.redo=FALSE){
    # stop('Not implemented')
    relevant.TF.symbols <- 
        c(
            'TFAP2A',
            'SMARCA4',
            'TP53BP1',
            'SMARCA2',
            'SOX2',
            'KDM6A',
            'SOX3',
            'FOSL2',
            # 'CTCF',
            'EZH2',
            'MBD4',
            'RAD21',
            'EP300',
            'ZEB1',
            'TAF1',
            'MECP2'
        )
    check_cached_results(
        results_file=FILTERED_TF_SITES_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                RAW_TF_SITES_FILE %>%
                read_tsv(show_col_types=FALSE) %>%
                dplyr::rename(
                    'chr'='#CHROM',
                    'start'=START,
                    'end'=END,
                    'ID'=id,
                    'TF.ClassID'=tfClassId,
                    'TF.Symbol'=tfTitle,
                    'UniprotID'=uniprotId
                ) %>%
                filter(grepl('wildtype|WT|Genotype: control', treatment.set)) %>% 
                filter(!grepl('cancer|vehicle', treatment.set)) %>% 
                filter(grepl('neuron|neural|brain|NSC|iN', cell.set)) %>% 
                filter(!grepl('kidney|leukemia|lung|astoid|astoma|carcinoma', cell.set)) %>% 
                filter(TF.Symbol %in% relevant.TF.symbols) # %>% {.} -> TFs.df
            }
    )
    # TFs.df %>% filter(!grepl('kidney|leukemia|lung|astoid|astoma|carcinoma', cell.set)) %>% count(cell.set) %>% arrange(desc(n))
    # TFs.df %>% filter(grepl('neuron|neural|brain|NSC|iN', cell.set)) %>% count(cell.set) %>% arrange(desc(n))
    # TFs.df %>% filter(grepl('neuron|neural|brain|NSC|iN', cell.set)) %>% distinct(cell.set) %>% arrange(cell.set) %>% print(n=Inf)
    # TFs.df %>% count(treatment.set) %>% arrange(desc(n))
    # TFs.df %>% count(peak-caller.set) %>% arrange(desc(n))
    # TFs.df %>% count(TF.Symbol) %>% arrange(desc(n))
    # TFs.df %>% 
    #     # filter(grepl('control|wildtype|none|untreated', treatment.set)) %>% 
    #     filter(grepl('wildtype|WT|Genotype: control', treatment.set)) %>% 
    #     filter(!grepl('cancer|vehicle', treatment.set)) %>% 
    #     filter(grepl('neuron|neural|brain|NSC|iN', cell.set)) %>% 
    #     filter(!grepl('kidney|leukemia|lung|astoid|astoma|carcinoma', cell.set)) %>% 
    #     # count(TF.Symbol) %>% arrange(desc(n))
    #     filter(TF.Symbol %in% relevant.TF.symbols) %>%
    #     distinct(cell.set) %>% print(n=Inf)
    #     distinct(treatment.set) %>% print(n=Inf)
    # # All FDR DEGs in at least 1 comparison in RNA-Seq + TFBS in WT neural cell
    # load_all_DESeq2_results() %>% 
    #     filter(symbol %in% relevant.TF.symbols) %>%
    #     filter(pvalue < 0.05) %>% count(symbol) %>% print(n=Inf)
    #     filter(padj < 0.1) %>% count(symbol) %>% print(n=Inf)
    #     count(padj < 0.1, pvalue < 0.05, comparison)
}

################################################################################
# Compute binwise summaries for all FGE signals 
################################################################################
map_FGEs_to_bins <- function(
    sites.df,
    bins.df){
    # bins.filepath){
    # paste0('row.index=1; ', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '))
    # map site to bins they are inside
    sites.df %>% 
    right_join(
        # bins.filepath %>% read_tsv(),
        bins.df,
        suffix=c('.FGE', ''),
        by=
            join_by(
                # chr == chrom,
                chr,
                within(x$start, x$end, y$start, y$end)
            )
    )
}

compute_binwise_motif_signals <- function(overlaps.df){
    overlaps.df %>% 
    mutate(log10_qvalue=-log10(qvalue)) %>% 
    group_by(chr, start, end) %>% 
    summarize(
        FGE.overlaps.n=n(),
        across(
            .cols=
                c(
                    score,
                    # pvalue,
                    # qvalue,
                    log10_qvalue
                ),
            .fn=
                list(
                    'min'=min,
                    'q25'=partial(stats::quantile, probs=0.25, na.rm=TRUE),
                    'mean'=mean,
                    'median'=median,
                    'var'=var,
                    'q75'=partial(stats::quantile, probs=0.75, na.rm=TRUE),
                    'max'=max,
                    'total'=sum
                ),
            .names="FGE.{.col}.{.fn}"
        )
    ) %>%
    ungroup()
}

compute_binwise_cCRE_signals <- function(overlaps.df){
    overlaps.df %>% 
    count(
        chr, start, end,
        name='FGE.overlaps.n'
    )
}

compute_binwise_FGE_signals <- function(
    FGE.type, 
    sites.df,
    bins.df,
    # bins.filepath,
    ...){
    map_FGEs_to_bins(
        sites.df,
        # bins.filepath
        bins.df
    ) %>% 
    {
        if (FGE.type == 'cCRE') {
            compute_binwise_cCRE_signals(overlaps.df=.)
        } else if (FGE.type == 'CTCF') {
            compute_binwise_motif_signals(overlaps.df=.)
        } else if (FGE.type == 'TFBS') {
            compute_binwise_motif_signals(overlaps.df=.)
        } else {
            stop(glue('Invalid FGE.type: {FGE.type}'))
        }
        # case_when(
        #     FGE.type == 'cCRE' ~ list(compute_binwise_cCRE_signals(overlaps.df=.)),
        #     FGE.type == 'CTCF' ~ list(compute_binwise_motif_signals(overlaps.df=.)),
        #     FGE.type == 'TFBS' ~ list(compute_binwise_motif_signals(overlaps.df=.)),
        #     .unmatched='error'
        # )
    } # %>% {.[[1]]}
}

compute_all_binwise_signals <- function(
    all.raw.FGE.data.df,
    genomic.bins.df,
    silence=FALSE,
    force.redo=FALSE){
    all.raw.FGE.data.df %>% 
    cross_join(genomic.bins.df) %>%
    mutate(
        results_file=
            file.path(
                FGE_SIGNAL_DIR,
                glue('FGE.type_{FGE.type}'),
                glue('FGE.subtype_{FGE.subtype}'),
                glue("resolution_{resolution}"),
                glue('binwise.signal.summary.stats.tsv')
            )
    ) %>% 
    # only try to generate results that dont exist already
    {
        if (!force.redo) {
            filter(., !file.exists(results_file))
        } else {
            .
        }
    } %>% 
    # pmap(
    future_pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=compute_binwise_FGE_signals,
        force_redo=force.redo,
        silence=silence,
        return_data=FALSE,
        .progress=TRUE
    )
}

list_all_binwise_signal_files <- function(){
    FGE_SIGNAL_DIR %>% 
    parse_results_filelist(suffix='binwise.signal.summary.stats.tsv') %>%
    select(-c(MatrixID))
}

################################################################################
# Hi-C Feature ~ Functional Genomic Element (FGE) Association  Testing
################################################################################
pivot_FGE_binwise_metrics <- function(
    results.df,
    statlist=NULL){
    results.df %>% 
    # Pivot CTCF summary stats to tidy-format
    pivot_longer(
         starts_with('FGE.'),
         names_prefix='FGE.',
         names_to='sumstat',
         values_to='value'
    ) %>% 
    separate_wider_delim(
        sumstat,
        delim=fixed('.'),
        names=c('metric', 'stat'),
        cols_remove=FALSE
    ) %>%
    {
        if (!is.null(statlist)) {
            filter(., stat %in% statlist)
        } else {
            .
        }
    }
}

assign_nearest_feature_to_bins <- function(bins.df, features.df){
    bins.df %>% 
    # match bins to closest downstream feature 
    left_join(
        features.df,
        relationship='many-to-many',
        by=join_by(chr, closest(start >= feature.end))
    ) %>% 
    # match bins to closest upstream feature
    left_join(
        features.df,
        relationship='many-to-many',
        suffix=c('.ds', '.us'),
        by=join_by(chr, closest(end <= feature.start))
    ) %>% 
    # check if any bins match feature exactly
    left_join(
        features.df %>%
        rename_with(.cols=-c(chr), .fn=~ str_replace(.x, '$', '.match')) %>% 
        add_column(isExactMatch=TRUE),
        relationship='many-to-many',
        by=join_by(chr, start == feature.start.match, end == feature.end.match)
    ) %>% 
    mutate(isExactMatch=ifelse(!is.na(isExactMatch), isExactMatch, FALSE)) %>% 
    # pick the closer TAD boundary for each bin
    mutate(
        # for all bins within each TAD, get distance to nearest feature (start or end)
        dist.to.feature.ds=ifelse(is.na(feature.start.ds), Inf, abs(feature.start.ds - start)),
        dist.to.feature.us=ifelse(is.na(feature.end.us),   Inf, abs(feature.end.us - end)),
        dist.to.nearest=
            case_when(
                isExactMatch                             ~  0,
                dist.to.feature.ds <= dist.to.feature.us ~ -dist.to.feature.ds,
                dist.to.feature.ds >  dist.to.feature.us ~  dist.to.feature.us
            ),
        nearest.boundary=
            case_when(
                isExactMatch                             ~ 'match',
                dist.to.feature.ds <= dist.to.feature.us ~ 'ds',
                dist.to.feature.ds >  dist.to.feature.us ~ 'us'
            ),
        # Create new column with only the data for the closest feature (upstream or downstream)
        across(
            ends_with('.match'), 
            ~ case_when(
                # if ds bin is closest get *.ds column data
                nearest.boundary == 'match' ~ .x,
                nearest.boundary == 'ds'    ~ get(str_replace(cur_column(), '.match', '.ds')), 
                nearest.boundary == 'us'    ~ get(str_replace(cur_column(), '.match', '.us'))
            ),
            .names="{str_remove(.col, '.match$')}"
        )
    ) %>% 
    # Re-categorize bins "close" to HiC Features as HiC Features for testing
    mutate(
        isHiCFeature=
            case_when(
                abs(dist.to.nearest / resolution) <= HiCFeatureRadius.bins ~ 'Feature',
                abs(dist.to.nearest / resolution) >  HiCFeatureRadius.bins ~ 'Not Feature' 
            ) %>% 
            # force isHiCFeature == TRUE  to be x group in t.test
            # so "greater" is testing if TAD > non-TAD in terms of CTCF sites overlaps
            # i.e. is the mean number of CTCF sites in/near TAD boundaries > than
            # all bins within that TAD but not near either boundary
            factor(levels=c('Feature', 'Not Feature'))
    ) %>% 
    select(
        -c(
            nearest.boundary,
            ends_with(c('.match', '.ds', '.us'))
        )
    )
}

compute_features_fisher_tests <- function(
    overlaps.df,
    resolution,
    n.FGE.min.thresh,
    ...){
    # Count overlap of bins with enough CTCF sites & bins close enough to HiC Features
    contingency.table <- 
        overlaps.df %>% 
        group_by(isHiCFeature) %>% 
        summarize(
            n.bins.wo.FGE=sum(n.overlaps <  n.FGE.min.thresh),
            n.bins.w.FGE=sum(n.overlaps  >= n.FGE.min.thresh)
        ) %>%
        select(n.bins.w.FGE, n.bins.wo.FGE) %>%
        as.matrix()
    # Calculate enrichment stat for test
    bins.overlap      <- contingency.table[1,1]
    bins.w.feature    <- bins.overlap + contingency.table[1,2]
    bins.w.FGE       <- bins.overlap + contingency.table[2,1]
    bins.total        <- sum(contingency.table)
    # Calculate enrichment stat for test
    expected          <- bins.w.FGE * (bins.w.feature / bins.total)
    variance_term_1   <- (bins.total - bins.w.FGE) / (bins.total - 1)
    variance_term_2   <- (bins.total - bins.w.feature) / (bins.total - 1)
    std_dev           <- sqrt(expected * variance_term_1 * variance_term_2)
    enrichment.zscore <- (bins.overlap - expected) / std_dev
    # Calculate fisher pvalue 
    fisher.test.row <- 
        contingency.table %>% 
        fisher.test(alternative='greater') %>% 
        tidy()
    # tidy data into single row tibble
    list(
        contingency.A=contingency.table[1,1],
        contingency.B=contingency.table[1,2],
        contingency.C=contingency.table[2,1],
        contingency.D=contingency.table[2,2],
        enrichment.zscore=enrichment.zscore,
        metric='overlaps',
        stat='n',
        test='fisher'
    ) %>% 
    as_tibble() %>% 
    bind_cols(fisher.test.row)
}

compute_features_t_tests <- function(overlaps.df) {
    overlaps.df %>% 
    # pivot metrics so can do 1 test per metric+stat (i.e. per row)
    pivot_FGE_binwise_metrics() %>% 
    group_by(metric, stat) %>% 
    summarize(
        # welch's t-test of the mean FGE metric is > near  TAD boundaries or not
        # results.t.test=
            t.test(
                value ~ isHiCFeature,
                alternative='greater' # only care if TADs are enriched for FGEs, not depleted
            ) %>%
            tidy()
    ) %>%
    ungroup() %>% 
    add_column(test='t.test')
}

compute_features_corr_tests <- function(overlaps.df) {
    overlaps.df %>% 
    pivot_FGE_binwise_metrics() %>% 
    group_by(metric, stat) %>% 
    # calcualte test results for each FGEs metric + stat combo
    summarize(
        # are TADs more enriched for FGEs signal closer or farther from boundaries
        results.pearson=
            cor.test(
                x=dist.to.nearest,
                y=value,
                method='pearson',
                exact=FALSE,
                alternative='greater' # want closer to TAD ~ more CTCF signal, so +ve signal only
            ) %>% 
            list(),
        results.spearman=
            cor.test(
                x=dist.to.nearest,
                y=value,
                method='spearman',
                exact=FALSE,
                alternative='greater' # want closer to TAD ~ more CTCF signal, so +ve signal only
        ) %>% 
        list(),
        results.kendall=
            cor.test(
                x=dist.to.nearest,
                y=value,
                method='kendall',
                exact=FALSE,
                alternative='greater' # want closer to TAD ~ more CTCF signal, so +ve signal only
        ) %>% 
        list()
    ) %>%
    ungroup() %>% 
    pivot_longer(
        starts_with('results.'),
        names_prefix='results.',
        names_to='test',
        values_to='test.results'
    ) %>% 
    rowwise() %>% 
    mutate(test.results=tidy(test.results)) %>% 
    unnest(test.results)
}

calculate_enrichment_tests <- function(
    bins.df,
    features.df,
    test.type,
    ...){
    # paste0('row.idx=1; ', paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.idx]]', sep='', collapse='; '))
    # row.idx=1; annotation=tmp$annotation[[row.idx]]; feature.type=tmp$feature.type[[row.idx]]; motif=tmp$motif[[row.idx]]; method=tmp$method[[row.idx]]; Sample.Group=tmp$Sample.Group[[row.idx]]; resolution=tmp$resolution[[row.idx]]; bins.df=tmp$bins.df[[row.idx]]; features.df=tmp$features.df[[row.idx]]; HiCFeatureRadius.bins=tmp$HiCFeatureRadius.bins[[row.idx]]; n.CTCF.min.thresh=tmp$n.CTCF.min.thresh[[row.idx]]
    # Re-classify bins as being "HiC Features" if they are close enoughy i.e. 
    # within HiCFeatureRadius.bins bins of the actual feature
    assign_nearest_feature_to_bins(
        bins.df,
        features.df
    ) %>%
    {
        # calculate fisher's exact test pvalue of whether bins that are at/near TAD boundaries
        # are more likely to have > n.CTCF.min.thresh CTCF sites overlapping them than bins
        # that are at least HiCFeatureRadius.bins bins away from a TAD boundary
        if (test.type == 'fisher') {
            compute_features_fisher_tests(
                overlaps.df=.,
                ...
            ) %>%
            select(-c(method, alternative))
        } else if (test.type == 't.test') {
            # Directly test if summary stats over CTCF site qvalues/scores are different closer to Features
            compute_features_t_tests(overlaps.df=.) %>% 
            select(-c(estimate1, estimate2, parameter, statistic, method, alternative))
        # Test if distance to feature is correlated with CTCF stats
        } else if (test.type == 'corr.test') {
            compute_features_corr_tests(overlaps.df=.)
        } else {
            stop(glue('invalid test type: {test.type}'))
        }
    } %>% 
    dplyr::rename_with(.cols=-c('test'), .fn=~str_replace(., '^', 'test.'))
}

calculate_all_feature_enrichments <- function(
    overlaps.df,
    p.corr.group.cols=c(),
    ...){
    # overlaps.df=TAD.CTCF.overlaps.df;
    # overlaps.df=TAD.CTCF.overlaps.df %>% cross_join(expand_grid(HiCFeatureRadius.bins=c(0, 1, 2, 3), n.CTCF.min.thresh=c(1, 5, 10, 20, 30, 40) ));
    overlaps.df %>% 
            # {.} -> tmp; tmp
    # calculate enrichment test results across all conditions + params + hyper-params
    mutate(
        test.results=
            pmap(
            # future_pmap(
                .l=.,
                .f=calculate_enrichment_tests,
                ...,
                .progress=TRUE
            )
    )  %>% 
    unnest(test.results) %>% 
    select(-c(bins.df, features.df)) %>% 
    # correct tests across resolutions + TAD calling methods + tests
    group_by(
        across(
            all_of(
                intersect(
                    colnames(.), 
                    c(
                        'test',
                        'test.metric',
                        'test.stat',
                        p.corr.group.cols
                    )
                )
            )
        )
    ) %>% 
    mutate(p.adj=p.adjust(test.p.value, method='BH')) %>% 
    ungroup() %>%
    # calculate log pvalues for plotting
    mutate(
        log.p.value=-log10(test.p.value),
        log.p.adj=-log10(p.adj)
    )
}

post_process_test_results <- function(
    test.df,
    sig.expr=NULL){
    test.df %>% 
    {
        if (!is.null(sig.expr)) {
            mutate(., is.Enriched=ifelse(!! rlang::parse_expr(sig.expr), sig.expr, 'N.S.'))
        } else {
            .
        }
    } %>% 
    {
        if ('motif' %in% colnames(.)) {
            mutate(., motif=str_remove(motif, '_CORE_vertebrates'))
        } else {
            .
        }
    } %>% 
    {
        if ('HiCFeatureRadius.bins' %in% colnames(.)) {
            mutate(
                .,
                HiCFeatureRadius.bins.fct=
                    glue('<= {HiCFeatureRadius.bins} away') %>% 
                    fct_reorder(HiCFeatureRadius.bins)
            )
        } else {
            .
        }
    {
    } %>% 
        if ('n.FGE.min.thresh' %in% colnames(.)) {
            mutate(
                .,
                n.FGE.min.thresh.fct=
                    glue('>= {n.FGE.min.thresh} FGEs') %>% 
                    fct_reorder(n.FGE.min.thresh)
            )
        } else {
            .
        }
    } %>% 
    dplyr::rename(
        'metric'=test.metric,
        'stat'=test.stat
    )
}

