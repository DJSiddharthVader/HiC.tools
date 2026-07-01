######################################################################
# Depdendencies
######################################################################
library(here)
here::i_am('scripts/weiner.replication/make.replication.figures.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'utils.annotations.R'))
    source(file.path(SCRIPT_DIR, 'weiner.replication/utils.weiner.replication.R'))
    library(magrittr)
    library(tidyverse)
})

######################################################################
# Comparison parameters
######################################################################
options(dplyr.summarize.inform=FALSE, scipen=999)
# Comparison hyper-params
resolutions <- c(100, 50, 10) * 1e3
hyper.params.df <- 
    expand_grid(
        # test.method=c('t.test', 'wilcox.test', 'anova', 'kruskal.test'),
        test.method=c('t.test', 'wilcox.test'),
        normalization=c('raw', 'balanced'),
        resolution=resolutions,
        comparison.group=c('regions', 'Genotype')
    ) # %>% 
# list of regions to compare
region.comparisons <- 
    tribble(
        ~region.P2,    ~region.P1,     ~region.background,
        "16p11.2 CNV", "16p Telomere", "16p",
        "16p11.2 CNV", "16p Telomere", "16q",
        "16p11.2 CNV", "16p Telomere", "chr16",
        "16p11.2",     "16p Telomere", "16p",
        "16p11.2",     "16p Telomere", "16q",
        "16p11.2",     "16p Telomere", "chr16",
        # "16p",         "16p",          "16q",
        # "16q",         "16p",          "16p"
    )

######################################################################
# Set up comparison inputs
######################################################################
# list all comparisons across Samples, region-pairs and hyper-params
all.comparisons.df <- 
    list_all_region_comparisons(
        huper.params.df,
        region.comparisons
    ) %>% 
    mutate(
        resolution.label=scale_numbers(resolution, force_chr=TRUE),
        region.comparison=glue('{region.P2} - {region.P1} [{region.background}]'),
        results_file=
            file.path(
                WEINER_REPLICATION_PLOTS_DIR,
                glue('comparison_{comparison.group}'),
                glue('normalization_{normalization}'),
                glue('isMerged_{isMerged}'),
                glue('resolution_{resolution}'),
                glue('{region.P1}-{region.P2}-{region.background}-region.vs.distanct.matched.ctrl.pdf')
            ) %>%
            str_replace_all(' ', '.')
    ) %>% 
    filter(test.method == 't.test') %>% 
    filter(normalization == 'raw') %>% 
    filter(isMerged) %>% 
    filter(resolution %in% c(100000)) %>% 
    separate_wider_delim(
        SampleID,
        delim='.',
        names=c(NA, 'Celltype', 'Genotype', NA, NA), 
        cols_remove=FALSE
    )
    # {.}; all.comparisons.df

######################################################################
# Compare regions differences
######################################################################
txt.size=10
all.comparisons.df %>%
    filter(comparison.group == 'regions') %>% 
    filter(Genotype == 'WT') %>% 
    select(-c(Celltype, Genotype)) %>% 
    nest(
        plot.df=
            -c(
                # stat.var, facet.row, 
                # SampleID,
                region.comparison,
                results_file,
                resolution, resolution.label,
                normalization,
                test.method,
                isMerged
            )
    ) %>% 
    pmap(
        .l=.,
        .f=plot_contacts_regions_boxplot,
        stat.var='group',
        y.var='IF',
        fill.var='group',
        facet.row='Genotype',
        facet.col='Celltype',
        scales='free_y',
        width=6,
        height=4,
        label.x.npc=c(0.25),
        label.y.npc=c(0.48),
        n_size=3,
        exclude_background=FALSE,
        ylim=TRUE,
        axis.text.x=element_text(size=txt.size, angle=35, hjust=1),
        axis.text.y=element_text(size=txt.size),
        strip.text=element_text(size=txt.size),
        legend.text=element_text(size=txt.size),
        legend.key.size=unit(10, 'pt'),
        axis.title.x=element_blank(),
        .progress=TRUE
    )

######################################################################
# Compare genotype differences
######################################################################
all.comparisons.df %>%
    filter(comparison.group == 'Genotype') %>% 
    select(-c(Celltype, Genotype)) %>% 
    nest(
        plot.df=
            -c(
                # SampleID,
                region.comparison,
                results_file,
                resolution, resolution.label,
                normalization,
                test.method,
                isMerged
            )
    ) %>% 
    pmap(
        .l=.,
        .f=plot_contacts_regions_boxplot,
        stat.var='Genotype',
        y.var='IF',
        fill.var='Genotype',
        facet.col='Celltype',
        facet.row='group',
        scales='free_y',
        width=6,
        height=4,
        n_size=3,
        exclude_background=TRUE,
        ylim=FALSE,
        axis.text.x=element_text(size=txt.size),
        axis.text.y=element_text(size=txt.size),
        strip.text=element_text(size=txt.size),
        legend.text=element_text(size=txt.size),
        legend.key.size=unit(10, 'pt'),
        .progress=TRUE
    )

