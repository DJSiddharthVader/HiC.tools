################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/Delta.Expression.Association.Testing/link.diff.HiFs.to.Genes.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.HiFs.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.association.R'))
    # source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    # source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    # source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    # source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
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
# load DEG results for all conditions
deg.results.df <- 
    load_all_DESeq2_results() %>% 
    as_iranges()

################################################################################
# Generate "Direct" differential HiF ~ Gene associations
################################################################################
# repeat same but for differential features i.e. diff.HiFs are from a pair of conditions being compared
# instead of HiFs, which are from individual conditions
all.diff.HiFs.dfs <- 
    list(
        # 'compartment.switches',
        # 'compartment.region',
        # 'loop.nesting',
        # 'loop.anchor',
        # 'loop',
        # 'DIR.anchor',
        'DIR',
        'TAD.Boundary',
        'TAD'
    ) %>% 
    combine_all_between_condition_HiFs(association.params.df=association.params.df)

################################################################################
# Generate "Direct" HiF ~ Gene associations
################################################################################
# gene.positions.df <- load_gene_annotations() %>% as_iranges()
all.HiFs.df %>% 
    add_column(association.type='Direct') %>% 
    add_column(association.subtype='Direct') %>% 
    add_column(association.source='Direct') %>% 
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
                glue('normalization_{normalization}'),
                glue('{SampleID}-HiF.Gene.Associations.tsv')
            )
    ) %>% 
        {.} -> tmp; tmp
    # map every gene within/narby
    pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=map_HiFs_to_genes_directly,
        deg.results.df=deg.results.df,
        .progress=TRUE
    )

################################################################################
# Generate "Indirect" HiF ~ Gene associations i.e. those defined by functional elements
################################################################################

