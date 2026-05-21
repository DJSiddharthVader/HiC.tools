################################################################################
# Dependencies
################################################################################

################################################################################
# eQTL data
################################################################################
get_remote_filepaths_for_eQTLs_of_interest <- function(){
    EQTL_REMOTE_FILEPATHS_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    filter(sample_group %in% c()) %>% 
    filter(tissue_label %in% c()) %>% 
    filter(condition_label %in% c()) %>% 
    filter(quant_method %in% c()) %>% 
    select(
        # study_id,
        # dataset_id,
        study_label,
        sample_group,
        # tissue_id,
        tissue_label,
        condition_label,
        sample_size,
        # quant_method
        ftp_cs_path  # path to credible set of eQTLs
    ) %>%
    dplyr::rename('eQTL.filpath'=ftp_cs_path)
}

download_and_parse_eQTL_file <- function(eQTL.filpath){
    tmp.filepath <- tempfile()
    EQTL_REMOTE_FILEPATHS_LINK %>% download.file(tmp.filepath)
    tmp.filepath %>%
    {.}
}

download_and_parse_all_eQTL_files <- function(force_redo=FALSE){
    get_remote_filepaths_for_eQTLs_of_interest() %>%
    mutate(results_file=file.path(CACHED_EQTLS_DIR, basename(eQTL.filpath))) %>% 
    pmap(
        .f=check_cached_results,
        results_fnc=download_and_parse_eQTL_file,
        force_redo=force_redo,
        .progress=TRUE
    ) %>%
    bind_rows()
}

################################################################################
# ABC data
################################################################################
load_internal_ABC_enhancers <- function(){
    INTERNAL_ABC_SCORES_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    # dplyr::rename(
    #     'association.subtype'=,
    #     'EnhancerID'=,
    #     'Target.Gene.TSS.Enhancer.distance'=,
    #     'Target.Gene.TSS'=,
    #     'Target.Gene.Symbol'=
    # ) %>% 
    add_column(
        association.source='Internal',
        association.type='ABC.enhancer'
    )
}

load_nasser_ABC_enhancers <- function(){
    # every row is the position of an enhancer and what gene it is linked to and how
    NASSER_ABC_SCORES_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    # Elements with an ABC score > 0.015 are typically considered "significant" connections. 
    # according to Nasser et al. 2021
    filter(ABC.Score > 0.15) %>% 
    mutate(
        association.subtype=
            case_when(
                isSelfPromoter ~ 'self.promoter',
                !isSelfPromoter ~ glue('enhancer.linked.{class}')
            )
    ) %>% 
    dplyr::rename(
        'EnhancerID'=name,
        'Target.Gene.TSS.Enhancer.distance'=distance,
        'Target.Gene.TSS'=TargetGeneTSS,
        'Target.Gene.Symbol'=TargetGene
    )
}

################################################################################
# Utils
################################################################################
load_gene_annotations <- function(gene.types=NULL){
    # gene types
    GENE_ANNOTATIONS_FILE %>%
    read_tsv(
        show_col_types=FALSE,
        col_names=
            c(
                'chr',
                'start',
                'end',
                'strand',
                'EnsemblID',
                'Gene.Symbol',
                'Gene.Type'
            )
    ) %>%
    {
        if (!is.null(gene.types)) {
            filter(., Gene.Type %in% gene.types)
        } else {
            .
        }
    }
    # 19913 protein_coding
    # 10214 processed_pseudogene
    #  7484 lincRNA
    #  5497 antisense
    #  2662 unprocessed_pseudogene
    #  2221 misc_RNA
    #  1909 snRNA
    #  1879 miRNA
    #  1067 TEC
    #   943 snoRNA
    #   898 sense_intronic
    #   853 transcribed_unprocessed_pseudogene
    #   555 processed_transcript
    #   549 rRNA
    #   472 transcribed_processed_pseudogene
    #   188 IG_V_pseudogene
    #   183 sense_overlapping
    #   144 IG_V_gene
    #   123 transcribed_unitary_pseudogene
    #   106 TR_V_gene
    #    95 unitary_pseudogene
    #    79 TR_J_gene
    #    49 scaRNA
    #    47 bidirectional_promoter_lncRNA
    #    38 polymorphic_pseudogene
    #    37 IG_D_gene
    #    33 TR_V_pseudogene
    #    32 3prime_overlapping_ncRNA
    #    22 pseudogene
    #    22 Mt_tRNA
    #    18 IG_J_gene
    #    14 IG_C_gene
    #     9 IG_C_pseudogene
    #     8 ribozyme
    #     6 TR_C_gene
    #     5 sRNA
    #     4 TR_J_pseudogene
    #     4 TR_D_gene
    #     3 non_coding
    #     3 IG_J_pseudogene
    #     2 translated_processed_pseudogene
    #     2 Mt_rRNA
    #     1 vaultRNA
    #     1 scRNA
    #     1 macro_lncRNA
    #     1 IG_pseudogene
}

                        )
                )
    ) %>% 
        by=
            join_by(
                chr,
                Target.Gene.Symbol == Gene.Symbol
            )
    )
}

################################################################################
# Directly link genes to HiFs by overlap/proximity
################################################################################
map_HiFs_to_genes_directly <- function(
    HiFs.df,
    gene.positionsd.df,
    min.overlap.frac,
    recip.min.overlap.frac,
    ...){
    HiFs.df %>% 
    {.}
}

################################################################################
# Indirectly link genes to HiFs by overlap/proximity with gene-associated functinoal elements
################################################################################

