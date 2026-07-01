library(here)
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(BiocParallel)
    library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
    library(magrittr)
    library(tidyverse)
})
# load DIR data
differential.contacts.df <- 
    check_cached_results(
        results_file=ALL_MULTIHICCOMPARE_RESULTS_FILE,
        # force_redo=TRUE,
        results_fnc=load_all_multiHiCCompare_results,
        sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC,
        sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS,
        # resolutions=c(100, 50, 25, 10) * 1e3,
        resolutions=c(100, 50, 25) * 1e3,
        # resolutions=c(100, 50) * 1e3,
        gw.fdr.threshold=0.1,
        fdr.threshold=0.1,
        nom.threshold=0.05
    ) %>% 
    standardize_data_cols() %>%
    rename_with(~ str_replace(.x, 'Sample.Group', 'SampleID')) %>% 
    MAKE_SAMPLE_GROUP_PAIR_TYPE_COLUMNS() %>% 
    rename_with(~ str_replace(.x, 'SampleID', 'Sample.Group'))
plot.df <-
    differential.contacts.df %>% 
    select(
        resolution, merged,
        # Edit.Numerator, Edit, Genotype, Comparison.Group, Sample.Group,
        Edit.Numerator, Comparison, Comparison.Group,
        chr,
        logFC, log.p.adj.gw
    ) %>% 
    filter(resolution %in% c('100Kb', '50Kb')) %>% 
    filter(Comparison.Group  %in% c('Edit.MT vs  All.WT', 'Edit.MT vs Edit.WT')) %>% 
    mutate(merged=ifelse(merged == 'Merged', 'Merged Matrices', 'Individual Replicates')) %>% 
    mutate(
        Comparison.Group=
            str_replace_all(Comparison.Group, 'MT', 'DEL') %>% 
            str_replace_all('Edit', 'GOG') %>% 
            str_replace_all(fixed('.'), ' ')
    ) %>% 
    mutate(
        Comparison=
            factor(
                str_replace_all(Comparison, fixed('.'), ' '),
                str_replace_all(levels(Comparison), fixed('.'), ' ')
            )
    ) %>% 
    # mutate(numerator=str_split_i(Comparison, ' vs ', 1) %>% unlist())
# plot.df %>% pull(numerator) %>% unique()
    mutate(
        numerator=str_split_i(Comparison, ' vs ', 1) %>% unlist(),
        numerator=
            factor(
                paste(str_replace_all(numerator, fixed('.'), ' '), 'vs WT'),
                levels=c('CTCF BIALLELIC vs WT', 'CTCF DEL vs WT', 'RAD21 DEL vs WT', 'WAPL DEL vs WT', 'NIPBL DEL vs WT')
            )
    )
# make boxplot of logFC of DIRs
make_base_plot <- function(
    data.df,
    outlier.size=0.1,
    strip.size=9,
    x.size=7,
    y.title.size=7){
    my_comparisons <- 
        list(
            c('CTCF BIALLELIC vs WT', 'CTCF DEL vs WT'),
            c('CTCF BIALLELIC vs WT', 'RAD21 DEL vs WT'),
            # c('RAD21 DEL', 'WAPL DEL'),
            c('WAPL DEL vs WT', 'NIPBL DEL vs WT')
        )
    # my_comparisons <- 
    #     data.df[['numerator']] %>%
    #     unique() %>%
    #     as.character() %>% 
    #     combn(m=2, simplify=FALSE)
    ggplot(
        data.df,
        aes(
            x=numerator,
            y=logFC,
            color=Edit.Numerator,
        )
    ) +
    geom_boxplot(
        aes(fill=Edit.Numerator),
        color='black',
        outliers=FALSE
    ) +
    geom_jitter(
        data=
            data.df %>%
            group_by(Comparison.Group, merged, Comparison, resolution) %>% 
            filter(logFC < quantile(logFC, 0.25) - 1.5 * IQR(logFC) | logFC > quantile(logFC, 0.75) + 1.5 * IQR(logFC)) %>% 
            ungroup(), 
        show.legend=FALSE,
        size=outlier.size
    ) +
    # Make boxplot of counts
        stat_compare_means(
            comparisons=my_comparisons,
            size=2,
            aes(label=paste0('t-test p-value=', after_stat(p.format))),
            # label='p.signif',
            label.x.npc='left',
            label.y.npc='top',
        ) +
    ylab('log(FC) Contact Frequency') +
    labs(fill='GOG Model') +
    geom_hline(
        yintercept=0,
        linewidth=0.5,
        color='black',
        linetype='dashed'
    ) +
    # theme_classic() +
    theme(
        legend.position='top',
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
make_boxplot_batch <- function(data.df, mode, ...){
    data.df %>% 
    mutate(across(c(Comparison.Group, resolution, merged,), function(x) {as.factor(as.character(x))})) %>% 
    {
        if (mode == 'all.wt'){ 
            filter(., !grepl('Merged', merged)) %>% 
            make_base_plot(...) +
            facet_nested(
                resolution ~ Comparison.Group,
                independent='y',
                scales='free'
            )
        } else if (mode == 'merged'){ 
            filter(., !grepl('All', Comparison)) %>% 
            make_base_plot(...) +
            facet_nested(
                resolution ~ merged,
                independent='y',
                scales='free'
            )
        } else if (mode == 'neither'){ 
            filter(., !grepl('Merged', merged)) %>% 
            filter(!grepl('All', Comparison)) %>% 
            make_base_plot(...) +
            facet_nested(
                resolution ~ .,
                independent='y',
                scales='free'
            )
        } else {
            make_base_plot(., ...) +
            facet_nested(
                resolution ~ merged + Comparison.Group,
                independent='y',
                scales='free'
            )
        }
    }
}
# plot.df %>% 
#     make_boxplot_batch(mode='all.wt') %>%
#     ggsave(
#         file.path(SCRIPT_DIR, 'all-cohesin.grant.2026.DIR.logFC.boxplot.pdf'), ., 
#         height=5, width=7, unit='in'
#     )
plot.df %>% 
    make_boxplot_batch(mode='neither') %>%
    ggsave(
        file.path(SCRIPT_DIR, 'neither-cohesin.grant.2026.DIR.logFC.boxplot.pdf'), ., 
        height=5, width=5, unit='in'
    )
# plot.df %>% 
#     make_boxplot_batch(mode='both') %>%
#     ggsave(
#         file.path(SCRIPT_DIR, 'both-cohesin.grant.2026.DIR.logFC.boxplot.pdf'), ., 
#         height=5, width=7, unit='in'
#     )
