################################################################################
# Dependencies
################################################################################

################################################################################
# Compute statistical testing results for HiF ~ DEG associations
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

compute_fisher_tests_within <- function(
    HiF.DEG.associations.df,
    DEG.p.adj.thresh=0.1,
    resolution,
    ...){
    ### NEED TO SEPARATEBY DEG COMAPRISON IN FUNCTIONAS OUTSIDE OF THIS
    # Count overlap of genes that are DEGs and associated with a HiF
    contingency.table <- 
        HiF.DEG.associations.df %>% 
        mutate(isDEG=ifelse(padj < DEG.p.adj.thresh, 'DEG', 'N.S')) %>% 
        group_by(isDEG) %>% 
        {
            # if associations are defined as the gene/gene-linked locus being 'within' the HiF
            if (association.strategy == 'within') {
                summarize(
                    .,
                    n.genes.within.a.HiF=sum(!is.na(FeatureID)),
                    n.genes.outside.all.HiF=sum(is.na(FeatureID))
                ) %>%
                select(n.genes.within.a.HiF, n.genes.outside.all.HiF)
            # # if associations are defined as the gene/gene-linked locus being 'nearby' the HiF
            # } else if (association.strategy == 'nearby') {
            #     mutate(
            #         .,
            #         is.HiF.Near.Associated.Gene=
            #             (HiF.to.Associated.Gene.Distance / resolution) < nearby.dist.thresh.bins
            #     ) %>% 
            #     summarize(
            #         n.genes.near.a.HiF=sum(is.HiF.Near.Associated.Gene),
            #         n.genes.far.from.all.HiFs=sum(!is.HiF.Near.Associated.Gene)
            #     ) %>%
            #     select(n.genes.near.a.HiF, n.genes.far.from.all.HiFs)
            } else {
                stop(glue('Invalid association.strategy: {association.strategy}'))
            }
        } %>% 
        as.matrix()
    # Calculate enrichment stat for fisher test over genes
    DEGs.associated.with.a.HiF      <- contingency.table[1,1]
    all.DEGs                        <- DEGs.associated.with.a.HiF + contingency.table[1,2]
    all.genes.associated.with.a.HiF <- DEGs.associated.with.a.HiF + contingency.table[2,1]
    all.genes                       <- sum(contingency.table)
    # Calculate enrichment stat for test
    expected          <- all.genes.associated.with.a.HiF * (all.DEGs / all.genes)
    variance_term_2   <- (all.genes - all.DEGs                       ) / (all.genes - 1)
    variance_term_1   <- (all.genes - all.genes.associated.with.a.HiF) / (all.genes - 1)
    std_dev           <- sqrt(expected * variance_term_1 * variance_term_2)
    enrichment.zscore <- (DEGs.associated.with.a.HiF - expected) / std_dev
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

compute_t_tests <- function(
    gene.HiF.mappings,
    padj.thresh=0.1,
    alternative='greater', # only care if TADs are enriched for FGEs, not depleted
    ...){
    # gene.HiF.mappings %>% head(1) %>% t()
    gene.HiF.mappings %>% 
    mutate(isDEG=padj.DESeq2 < padj.thresh) %>% 
    mutate(is.gene.associated.with.HiF=!is.na(FeatureID)) %>% 
    # group_by(comparison.DESeq2) %>% 
    summarize(
        # welch's t-test of the mean FGE metric is > near  TAD boundaries or not
        results_ks.test=
            ks.test(
                value ~ is.gene.associated.with.HiF,
                alternative=alternative
            ) %>%
            tidy()
        results_t.test=
            t.test(
                value ~ is.gene.associated.with.HiF,
                alternative=alternative
            ) %>%
            tidy()
    ) %>%
    ungroup() %>% 
    add_column(test='t.test')
}

compute_corr_tests <- function(overlaps.df) {
    overlaps.df %>% 
    pivot_FGE_binwise_metrics() %>% 
    # calcualte test results for each FGEs metric + stat combo
    group_by(metric, stat) %>% 
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

calculate_association_test_results <- function(
    HiF.DEG.associations.df,
    test.type,
    ...){
    # paste(c('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()'), collapse='; ')
    # row.index=1; filepath=tmp$filepath[[row.index]]; association.type=tmp$association.type[[row.index]]; association.subtype=tmp$association.subtype[[row.index]]; association.source=tmp$association.source[[row.index]]; association.strategy=tmp$association.strategy[[row.index]]; resolution=tmp$resolution[[row.index]]; HiF.type=tmp$HiF.type[[row.index]]; HiF.scope=tmp$HiF.scope[[row.index]]; SampleID=tmp$SampleID[[row.index]]; HiF.DEG.associations.df=tmp$HiF.DEG.associations.df[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% t()
    # Re-classify bins as being "HiC Features" if they are close enoughy i.e. 
    # within HiCFeatureRadius.bins bins of the actual feature
    HiF.DEG.associations.df %>% 
    group_by(comparison) %>% 
    {
        # calculate fisher's exact test pvalue of whether bins that are at/near TAD boundaries
        # are more likely to have > n.CTCF.min.thresh CTCF sites overlapping them than bins
        # that are at least HiCFeatureRadius.bins bins away from a TAD boundary
        if (test.type == 'fisher') {
            compute_fisher_tests_within(
                HiF.DEG.associations.df=.,
                ...
            ) %>%
            select(-c(method, alternative))
        # } else if (test.type == 't.test') {
        #     # Directly test if summary stats over CTCF site qvalues/scores are different closer to Features
        #     compute_features_t_tests(HiF.DEG.associations.df=.) %>% 
        #     select(-c(estimate1, estimate2, parameter, statistic, method, alternative))
        # # Test if distance to feature is correlated with CTCF stats
        # } else if (test.type == 'corr.test') {
        #     compute_features_corr_tests(HiF.DEG.associations.df=.)
        } else {
            stop(glue('invalid test type: {test.type}'))
        }
    } %>% 
    dplyr::rename_with(.cols=-c('test'), .fn=~str_replace(., '^', 'test.'))
}

calculate_all_association_test_results <- function(
    all.HiF.DEG.associations.df,
    p.corr.group.cols=c(),
    ...){
    all.HiF.DEG.associations.df %>% 
            # {.} -> tmp; tmp
    # calculate enrichment test results across all conditions + params + hyper-params
    mutate(
        test.results=
            pmap(
            # future_pmap(
                .l=.,
                .f=calculate_association_test_results,
                ...,
                .progress=TRUE
            )
    )  %>% 
    unnest(test.results) %>% 
    select(-c(ends_with('.df'))) %>% 
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

################################################################################
# Stratify input associations by HiF + DEG specific features for testing more specific hypotheses
################################################################################

        # filter(association.status != 'HiF Only') %>% 
        # count(association.status, padj.DESeq2 < 0.01) %>% 
        # pivot_wider(names_from=`padj.DESeq2 < 0.01`, names_prefix='is.DEG.', values_from='n') %>%
        # dplyr::arrange(-dplyr::row_number()) %>% 
        # select(starts_with('is.DEG.')) %>%
        # relocate(is.DEG.TRUE, is.DEG.FALSE) %>% 
        # as.matrix() %>%
        # fisher.test(alternative='greater') %>%
        # tidy()
