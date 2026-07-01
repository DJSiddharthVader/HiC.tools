################################################################################
# Dependencies
################################################################################
library(here)
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.HiFs.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.association.R'))
    library(tidyverse)
    library(magrittr)
})
# Indirect gene associations
all.indirect.gene.links.df <- 
    enhancers.df %>% 
    inner_join(
        load_all_DESeq2_results(
            rename.for.comparison=TRUE,
            as.granges=FALSE
        ) %>% 
        select(-c(chr, start, end)),
        relationship='many-to-many',
        by=join_by(Target.Gene.Symbol)
    ) %>% 
    mutate(across(c(start, end), as.integer)) %>% 
    nest(
        associations.df=
            -c(
                comparison.DESeq2,
                association.source,
                association.subtype,
                association.type
            )
    ) %>%
    mutate(
        associations.df=
            pmap(
                .l=list(associations.df),
                .f=
                    function(associations.df) {
                        associations.df %>% 
                        dplyr::rename('seqnames'=chr) %>% 
                        as_granges()
                    }
            )
    )
# Direct gene associations
deg.results.df <- 
    load_all_DESeq2_results(
        rename.for.comparison=TRUE,
        as.granges=TRUE
    )
all.direct.gene.links.df <- 
    tribble(
        ~association.type, ~association.subtype, ~association.source,
        'Direct',          'Direct',             'Direct',           
    ) %>%
    cross_join(
        deg.results.df %>%
        as_tibble() %>% 
        nest(associations.df=-c(comparison.DESeq2)) %>%
        mutate(associations.df=pmap(.l=list(associations.df), .f=as_granges))
    )

all.mappings.df <- 
    bind_rows(
        all.direct.gene.links.df,
        all.indirect.gene.links.df
    ) %>% 
    filter(association.subtype %in% c('Direct', 'enhancer.linked.promoter')) %>% 
    inner_join(
        all.diff.TADs.df %>%
        mutate(comparison.DESeq2=glue('{Sample.Group.Numerator}_vs_{Sample.Group.Denominator}')),
        by='comparison.DESeq2'
    ) %>%
    mutate(
        plot.df=
            pmap(
                .l=.,
                .f=
                    function(associations.df, diff.TADs.df, ...){
                        # join_overlap_inner_within(
                        join_overlap_left_within(
                            # all.mappings.df$associations.df[[41]],
                            # all.mappings.df$diff.TADs.df[[41]],
                            associations.df,
                            diff.TADs.df,
                            suffix=c('.association', '.TAD')
                        ) %>%
                        as_tibble() %>%
                        mutate(
                            is.TAD.differential=
                                ifelse(is.na(is.TAD.differential), FALSE, is.TAD.differential)
                        ) %>% 
                        mutate(abs.logFC.DESeq2=abs(log2FoldChange.DESeq2)) %>% 
                        mutate(is.DEG=padj.DESeq2 < 0.1) %>%
                        dplyr::rename('chr'=seqnames) %>% 
                        select(-c(strand))
                    }
            )
    ) %>% 
    select(-c(comparison.DESeq2, associations.df, diff.TADs.df, Sample.Group.Numerator, Sample.Group.Denominator)) %>% 
    filter(resolution == 50000) %>% 
    select(-c(association.source, TADCompare.params, TAD.method, association.type)) %>% 
    unnest(plot.df) %T>%
    write_tsv(file.path(RESULTS_DIR, './deg.tad.enhancer.mappings.tsv'))

plot.df <- 
    read_tsv('./deg.tad.enhancer.mappings.tsv') %>%
    # dplyr::rename('Edit.Numerator'=Sample.Group.Edit) %>% 
    # dplyr::rename('logFC'=abs.logFC.DESeq2) %>% 
    dplyr::rename('logFC'=log2FoldChange.DESeq2) %>% 
    mutate(
        Edit.Numerator=
            glue('{Sample.Group.Edit} DEL vs WT') %>%
            factor(levels=c('CTCF DEL vs WT','RAD21 DEL vs WT','WAPL DEL vs WT','NIPBL DEL vs WT'))
    ) %>% 
    mutate(
        is.TAD.differential=ifelse(is.TAD.differential, 'Differential TAD', 'Non-Differential TAD'),
        is.DEG=ifelse(is.DEG, 'Gene p adj < 0.1', 'N.S.')
    )

make_base_plot2 <- function(
    data.df,
    outlier.size=0.1,
    strip.size=9,
    x.size=7,
    y.title.size=7){
    ggplot(
        data.df %>% filter(is.DEG == 'Gene p adj < 0.1'),
        aes(
            x=is.TAD.differential,
            y=logFC
        )
    ) +
    geom_boxplot(
        # fill='blue',
        # color='black',
        # outliers=TRUE
        outliers=FALSE
    ) +
    # Make boxplot of counts
    # stat_compare_means(
    #     comparisons=
    #         list(
    #             c('Differential TAD', 'Non-Differential TAD')
    #         ),
    #     size=2,
    #     aes(label=paste0('t-test p-value=', after_stat(p.format))),
    #     # label='p.signif',
    #     label.x.npc='left',
    #     label.y.npc='top',
    # ) +
    ylab('Gene Expression log2(FC)') +
    # labs(fill='DEG Status') +
    geom_hline(
        yintercept=0,
        linewidth=0.5,
        color='black',
        linetype='dashed'
    ) +
    facet_nested(
        association.subtype ~ Edit.Numerator,
        # independent='y',
        scales='free'
    ) +
    # theme_classic() +
    theme(
        legend.position='top',
        legend.title=element_blank(),
        strip.background=element_rect(fill='grey70'),
        strip.text=element_text(size=strip.size),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=x.size),
        axis.title.y=element_text(size=y.title.size),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        axis.line=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks.length=unit(3, "pt"),
        axis.title=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=10
            ),
        axis.text=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=8
            )
    )
}
plot.df %>% 
    # make_boxplot_batch(mode='neither') %>%
    # make_base_plot() %>%
    make_base_plot2() %>%
    ggsave(
        file.path('./direct.and.enhancer.mediated.TAD.DEG.association.pdf'), .,
        height=5, width=7, unit='in'
    )

# fisher testing
fisher.df <- 
    plot.df %>%
    nest(
        data.df=
            -c(
                resolution,
                association.subtype,
                Sample.Group.Edit
            )
    ) %>%
    mutate(
        fisher.results=
            pmap(
                .l=.,
                .f=
                    function(data.df, ...){
                        # data.df %>%
                        contingency.table <- 
                            # fisher.df$data.df[[8]] %>% 
                            data.df %>% 
                            count(is.DEG, is.TAD.differential) %>% 
                            pivot_wider(names_from=is.TAD.differential, values_from=n) %>%
                            select(-c(is.DEG)) %>% 
                            as.matrix()
                        # Calculate enrichment stat for fisher test over genes
                        if (all(dim(contingency.table) == c(2,2))) {
                            DEGs.in.dTAD      <- contingency.table[1,1]
                            all.DEGs          <- contingency.table[1,1] + contingency.table[1,2]
                            all.genes.in.dTAD <- contingency.table[1,1] + contingency.table[2,1]
                            all.genes         <- sum(contingency.table)
                            # Calculate enrichment stat for test
                            expected          <- all.genes.in.dTAD * (all.DEGs / all.genes)
                            variance_term_2   <- (all.genes - all.DEGs           ) / (all.genes - 1)
                            variance_term_1   <- (all.genes - all.genes.in.dTAD) / (all.genes - 1)
                            std_dev           <- sqrt(expected * variance_term_1 * variance_term_2)
                            enrichment.zscore <- (DEGs.in.dTAD - expected) / std_dev
                            # Calculate fisher pvalue 
                            fisher.test.row <- 
                                if (!is.na(min(contingency.table))) {
                                    contingency.table %>% 
                                    fisher.test(alternative='greater') %>% 
                                    tidy()
                                } else {
                                    tibble_row()
                                }
                            # tidy data into single row tibble
                            list(
                                DEGs.in.dTAD=contingency.table[1,1],
                                all.DEGs=contingency.table[1,1] + contingency.table[1,2],
                                all.genes.in.dTAD=contingency.table[1,1] + contingency.table[2,1],
                                genes.not.in.dTAD=contingency.table[2,2],
                                all.genes=sum(contingency.table),
                                enrichment.zscore=enrichment.zscore,
                                test='fisher'
                            ) %>% 
                            as_tibble() %>% 
                            bind_cols(fisher.test.row)
                        } else {
                            tibble_row()
                        }
                }
            )
    ) %>%
    unnest(fisher.results) %>% 
    mutate(
        Edit.Numerator=
            glue('{Sample.Group.Edit} DEL vs WT') %>%
            factor(levels=c('CTCF DEL vs WT','RAD21 DEL vs WT','WAPL DEL vs WT','NIPBL DEL vs WT'))
    )
# fisher.df %>% head(2) %>% t()

make_fisher_plot <- function(
    data.df,
    outlier.size=0.1,
    strip.size=9,
    x.size=7,
    y.title.size=7){
    ggplot(
        data.df,
        aes(
            x=Edit.Numerator,
            y=-log10(p.value),
            fill=association.subtype
        )
    ) +
    geom_col(position='dodge', color='black') +
    # ylab("Fishers Exact Test log(p-value)") +
    ylab("Fishers Exact Test -log10(p-value)") +
    geom_hline(
        yintercept=-log10(0.05),
        linewidth=0.5,
        color='black',
        linetype='dashed'
    ) +
    # facet_nested(
    #     association.subtype ~ Edit.Numerator,
    #     independent='y',
    #     scales='free'
    # ) +
    # theme_classic() +
    theme(
        legend.position='top',
        legend.title=element_blank(),
        strip.background=element_rect(fill='grey70'),
        strip.text=element_text(size=strip.size),
        axis.title.x=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=x.size),
        axis.title.y=element_text(size=y.title.size),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        axis.line=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks.length=unit(3, "pt"),
        axis.title=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=10
            ),
        axis.text=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=8
            )
    )
}

# fisher.df %>% t()
fisher.df %>% 
    make_fisher_plot() %>%
    ggsave(
        file.path('./fisher.direct.and.enhancer.mediated.TAD.DEG.association.pdf'), .,
        height=5, width=7, unit='in'
    )

