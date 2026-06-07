################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/Delta.Expression.Association.Testing/link.HiFs.to.Genes.R')
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

################################################################################
# Params to specify
################################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('resolutions', 'threads', 'force'),
        has.positional=FALSE
    )
if (parsed.args$threads > 1){
    # plan(multisession, workers=N_WORKERS_FOR_PARLLELIZATION)
    options(future.globals.maxSize=1.23 * 1024**3)
    plan(multisession, workers=parsed.args$threads)
} else {
    plan(sequential)
}
# Params neede to define if a HiF is associated with a gene for a given association strategy or not
# Some features are association specific
# parsed.args$resolutions <- c(50, 25) * 1e3
parsed.args$resolutions <- c(50) * 1e3
association.params.df <- 
    HiF.association.hyper.params.df %>% 
    cross_join(tibble(resolution=parsed.args$resolutions))
# load DEG results for all conditions
deg.results.df <- 
    load_all_DESeq2_results() %>% 
    as_iranges()

################################################################################
# Combine all HiFs + testing params into neat parsable table
################################################################################
# combine everything, keeps things tidy by nesting miscelaneous feature data into a single column
# will be unnested and saved as individual columns
all.HiFs.df <- 
    c(
        # 'compartment.switches',
        # 'compartment.region',
        # 'loop.nesting',
        # 'loop.anchor',
        # 'loop',
        'TAD.Boundary',
        'TAD'
    ) %>% 
    combine_all_per_condition_HiFs(association.params.df=association.params.df)
    
################################################################################
# Generate "Direct" HiF ~ Gene associations
################################################################################
all.HiFs.df %>% 
    add_column(association.type='Direct') %>% 
    add_column(association.subtype='Direct') %>% 
    add_column(association.source='Direct') %>% 
    # create nested directory of files with HiF ~ Gene mappings for easy parsing/loading
# Direct gene associations
all.direct.gene.links.df <- 
    load_gene_annotations() %>%
    standardize_data_cols() %>% 
    mutate(
        seqnames=Target.Gene.chr,
        start=Target.Gene.start,
        end=Target.Gene.end
        results_file=
            file.path(
                HIF_GENE_ASSOCIATION_MAPPING_DIR,
                glue('association.type_{association.type}'),
                glue('association.subtype_{association.subtype}'),
                glue('association.source_{association.source}'),
                glue('association.strategy_{association.strategy}'),
                glue('HiF.type_{HiF.type}'),
                glue('resolution_{resolution}'),
                glue('{SampleID}-HiF.Gene.Associations.tsv')
            )
    ) %>% 
    as_granges() %>% 
    tribble(
        ~association.type, ~association.subtype, ~association.source, ~associations.df,
        'Direct',          'Direct',             'Direct',            .
        {.} -> tmp; tmp
    # map every gene within/narby
    pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=map_HiFs_to_genes_directly,
        deg.results.df=deg.results.df,
        .progress=TRUE
    )
# DEG results for all genes
deg.results.df <- 
    prep_DESeq2_results_for_associations(force.redo=FALSE)
    # prep_DESeq2_results_for_associations(force.redo=TRUE)

################################################################################
# Generate "Indirect" HiF ~ Gene associations i.e. those defined by functional elements
################################################################################
# load clean functional element annotations data 
all.indirect.associations.df <- 
    ALL_CLEAN_GENE_LOCUS_ASSOCIATIONS_FILE %>% 
    nest(
        associations.df=
            -c(
                association.source,
                association.subtype,
                association.type
            )
    ) %>%
    mutate(associations.df=pmap(.l=list(associations.df), .f=as_iranges)) %>% 
# join HiFs to functional loci
all.HiFs.df %>% 
    inner_join(
        all.indirect.associations.df,
    ) %>% 
    # create nested directory of files with HiF ~ Gene mappings for easy parsing/loading
    mutate(
        results_file=
            file.path(
                HIF_GENE_ASSOCIATION_MAPPING_DIR,
                glue('association.type_{association.type}'),
                glue('association.subtype_{association.subtype}'),
                glue('association.source_{association.source}'),
                glue('association.strategy_{association.strategy}'),
                glue('HiF.type_{HiF.type}'),
                glue('resolution_{resolution}'),
                glue('{SampleID}-HiF.Gene.Associations.tsv')
            )
    ) %>% 
        {.} -> tmp; tmp
    # map every gene within/narby
    pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=map_HiFs_to_genes_indirectly,
        deg.results.df=deg.results.df,
        all.indirect.associations.df=all.indirect.associations.df,
        .progress=TRUE
    )

