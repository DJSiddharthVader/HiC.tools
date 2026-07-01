######################################################################
# Dependencies
######################################################################
library(here)
# here::i_am('scripts/Delta.Expression.Association.Testing/coallate.all.gene.associated.functioal.loci.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'FGE.Association.Testing/utils.enrichment.R'))
    source(file.path(SCRIPT_DIR, 'Delta.Expression.Association.Testing/utils.association.R'))
    library(tidyverse)
    library(magrittr)
})

######################################################################
# Load Association data i.e. functional associations between loci and Genes
######################################################################
# ABC-based enhancers-promoter pairs
enhancers.df <- 
    bind_rows(
        # load_internal_ABC_enhancers() %>% 
        #     add_column(association.source='Internal'),
        load_nasser_ABC_enhancers() %>%
            select(-c(CellType)) %>% 
            add_column(association.source='Nasser_2021')
    ) %>% 
    add_column(association.type='ABC.enhancer')
# eQTLs linking genes to loci with specific variants
eQTLs.df <- 
    check_cached_results(
        results_file=ALL_RELEVANT_EQTLS_FILE,
        force_redo=parsed.args$force.redo,
        # force_redo=TRUE,
        results_fnc=download_and_clean_all_eQTL_files,
        p.adj.thresh=0.1,
        force.redo=FALSE
    ) %>% 
    unite(
        'association.subtype',
        sep='#',
        remove=FALSE,
        c(
            tissue_label,
            condition_label,
            quant_method
        )
    ) %>% 
    add_column(association.source='EBI') %>% 
    add_column(association.type='eQTL')
# # TFs linking genes to the locations of TFs genes that target them
# TFs.df <- 
#     tibble()
# TF binding sites, linking genes to bindings sites of TFs that target them
TFBS.df <- 
    load_TF_binding_sites() %>% 
    # remove any TF with < 2K site genome-wide
    add_count(TF.Symbol) %>% filter(n > 2000) %>% select(-c(n)) %>% 
    select(
        -c(
            ID,
            antibody.set,
            cell.set,
            exp.set,
            treatment.set,
            `peak-caller.set`,
            `peak-caller.count`,
            `peak-caller.list`,
            # peak.count,
            # summit,
            TF.ClassID
        )
    ) %>% 
    dplyr::rename('Target.Gene.Symbol'=TF.Symbol) %>% 
    mutate('association.subtype'=Target.Gene.Symbol) %>% 
    add_column(association.source='GTRD') %>% 
    add_column(association.type='TFBS')

######################################################################
# Combine all functional elements and add gene metadata for downstream analysis
######################################################################
# gene.linked.functional.annotations=list(TFBS.df, eQTLs.df, enhancers.df)
# gene.annotations.df=load_gene_annotations()
# Now combine all annotations together
# for each gene-associated functional locus, join the gene position + metadata + EnsemblID
all.functional.loci.gene.associations <- 
    check_cached_results(
        results_file=ALL_CLEAN_GENE_LOCUS_ASSOCIATIONS_FILE,
        # force_redo=parsed.args$force.redo,
        force_redo=TRUE,
        results_fnc=nest_and_combine_all_functional_loci_dataset,
        gene.linked.functional.annotations=
            list(
                TFBS.df,
                eQTLs.df,
                # TFs.df,
                enhancers.df
            ),
        gene.annotations.df=load_gene_annotations(),
    )

