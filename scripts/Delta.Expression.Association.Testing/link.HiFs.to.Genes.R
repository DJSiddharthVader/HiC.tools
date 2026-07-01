######################################################################
# Dependencies
######################################################################
library(here)
# here::i_am('scripts/Delta.Expression.Association.Testing/link.HiFs.to.Genes.R')
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
# CLI args
parsed.args <- 
    handle_CLI_args(
        args=c('resolutions', 'threads', 'force'),
        has.positional=FALSE
    )
# parallelization
if (parsed.args$threads > 1){
    # options(future.globals.maxSize=1.23 * 1024**3)
    plan(multisession, workers=parsed.args$threads)
} else {
    plan(sequential)
}

######################################################################
# Combine all HiFs + testing params into neat parsable table
######################################################################
# combine all HiFs across HiF.types keeps things tidy by nesting feature-specific data into a single column
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
    combine_all_per_condition_HiFs()
# load functional element annotations i.e. functinoal loci linked each linked to a specific gene defined by 
# a specific functional mechanisms (association type/subtype) from a specific dataset (annotation source)
# Direct gene associations
all.direct.gene.links.df <- 
    load_gene_annotations() %>%
    standardize_data_cols() %>% 
    mutate(
        seqnames=Target.Gene.chr,
        start=Target.Gene.start,
        end=Target.Gene.end
    ) %>% 
    as_granges() %>% 
    tribble(
        ~association.type, ~association.subtype, ~association.source, ~associations.df,
        'Direct',          'Direct',             'Direct',            .
    )
# Indirect gene associations
all.indirect.gene.links.df <- 
    ALL_CLEAN_GENE_LOCUS_ASSOCIATIONS_FILE %>% 
    readRDS() %>% 
    filter(association.type %in% c('ABC.enhancer'))
# Parameters to guide mapping of HiFs to "associated" gene-linked functional loci
association.params.df <- 
    HiF.association.strategies.df %>% 
    cross_join(tibble(resolution=c(100, 50, 25) * 1e3))
    # cross_join(tibble(resolution=parsed.args$resolutions))
# DEG results for all genes
deg.results.df <- 
    prep_DESeq2_results_for_associations(force.redo=FALSE)
    # prep_DESeq2_results_for_associations(force.redo=TRUE)

######################################################################
# Map HiFs ~ Genes using Gene positions + gene-associated functional loci
######################################################################
# Now join all the input data together via matching relevant params 
# so now each row represents a specific set of HiF ~ Gene mappings (associations) to 
# save to an output file and use for downstream statistical testing
all.HiFs.df %>% 
    # Define association strategy for each type of Hi-C features
    inner_join(
        association.params.df,
        relationship='many-to-many',
        by=
            join_by(
                HiF.type,
                resolution
            )
    ) %>% 
    # Compare all sets of HiFs against all sets of gene-associated functional loci 
    cross_join(
        bind_rows(
            all.direct.gene.links.df,
            all.indirect.gene.links.df
        )
    ) %>% 
    # create output filepath using association metadata
    mutate(
        results_file=
            pmap_chr(
                .l=.,
                .f=make_per_condition_mapping_results_filepath,
                .progress=FALSE
            )
    ) %>% 
    # save mappings of HiFs to all associated genes + genes with 0 associated HiFs
    pmap(
        .l=.,
        .f=check_cached_results,
        results_fnc=associate_gene_links_to_HiFs,
        force_redo=TRUE,
        return_data=FALSE,
        deg.results.df=deg.results.df,
        .progress=TRUE
    )

