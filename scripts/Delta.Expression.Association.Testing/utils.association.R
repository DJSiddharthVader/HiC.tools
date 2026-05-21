################################################################################
# Dependencies
################################################################################
library(plyranges)

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
# Load Gene info + Differential results
################################################################################
load_gene_annotations <- function(
    gene.types=
        c(
            'lincRNA',
            'snRNA',
            'miRNA',
            'snoRNA',
            # 'transcribed_unprocessed_pseudogene',
            # 'processed_transcript',
            # 'transcribed_processed_pseudogene',
            'protein_coding'
        )
    ){
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
    } %>% 
    standardize_data_cols(skip.resolution=TRUE)
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

load_all_DESeq2_results <- function(force.redo=FALSE){
    check_cached_results(
        results_file=DESEQ2_RESULTS_FILE,
        force_redo=force.redo,
        results_fnc=
            function(suffix='-DESeq2.tsv') {
                # List all results files
                DESEQ2_DATA_DIR %>% 
                list.files(
                    pattern=suffix,
                    full.names=TRUE,
                    recursive=TRUE
                ) %>% 
                tibble(filepath=.) %>%
                mutate(info=str_remove(basename(filepath), suffix)) %>% 
                filter(!grepl('iPSC', info)) %>% 
                # Tidy pairwise metadata
                separate_wider_delim(
                    info,
                    # delim='-',
                    delim='_vs_',
                    names=c('Sample.Group.Numerator', 'Sample.Group.Denominator'),
                ) %>%
                # load DEG results
                mutate(
                    results=
                        pmap(
                            list(filepath),
                            read_tsv,
                            show_col_types=FALSE
                        )
                ) %>%
                unnest(results) %>%
                # rename_with(.f=~ str_replace(.x, '^', 'DESeq2.')) %>% 
                # # run TRADEtools to define transcriptome-wide effects
                # clean up columns
                dplyr::rename('EnsemblID'=ensemblid) %>% 
                mutate(gene.length=end - start) %>% 
                select(
                    -c(
                        filepath,
                        Sample.Group.Numerator, Sample.Group.Denominator,
                        row_index,
                        # strand,
                        geneid,
                        stat
                    )
                )
            }
    ) %>% 
    standardize_data_cols()
}

################################################################################
# Directly link genes to HiFs by overlap/proximity
################################################################################
map_HiF_to_genes_within <- function(
    HiFs.df,
    deg.results.df,
    min.Gene.overlap.frac,
    ...){
    deg.results.df %>% 
    # join_overlap(
    #     HiFs.df,
    left_join(
        HiFs.df,
        suffix=c('.Gene', '.HiF'),
        by=
            join_by(
                chr,
                within(start.Gene, end.Gene, start.HiF, end.HiF)
            )
    ) %>% 
    filter(Gene.overlap.frac > min.Gene.overlap.frac)
}

map_HiF_to_genes_nearby <- function(
    HiFs.df,
    deg.results.df,
    nearby.threshold,
    ...){
    # deg.results.df %>% 
    stop('Not Implemented')
}

map_HiFs_to_genes_directly <- function(
    HiFs.df,
    deg.results.df,
    association.strategy,
    ...){
    # paste0('row.index=1; ', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '), '; tmp %>% head(row.index) %>% tail(1) %>% t()')
    # row.index=161; resolution=tmp$resolution[[row.index]]; normalization=tmp$normalization[[row.index]]; SampleID=tmp$SampleID[[row.index]]; HiF.type=tmp$HiF.type[[row.index]]; HiFs.df=tmp$HiFs.df[[row.index]]; fuzzy.matching.threshold.bins=tmp$fuzzy.matching.threshold.bins[[row.index]]; association.strategy=tmp$association.strategy[[row.index]]; frac.reciprocal.matching.overlap=tmp$frac.reciprocal.matching.overlap[[row.index]]; association.type=tmp$association.type[[row.index]]; association.subtype=tmp$association.subtype[[row.index]]; association.source=tmp$association.source[[row.index]]; results_file=tmp$results_file[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% t()
    if (association.strategy == 'nearby') {
        map_HiF_to_genes_nearby(
            HiFs.df,
            deg.results.df,
            ...
        )
    } else if (association.strategy == 'within') {
        map_HiF_to_genes_within(
            HiFs.df,
            deg.results.df,
            ...
        )
    } else {
        stop(glue('Invalid association.strategy: {association.strategy}'))
    }
}

################################################################################
# Indirectly link genes to HiFs by overlap/proximity with gene-associated functinoal elements
################################################################################

