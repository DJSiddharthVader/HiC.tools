######################################################################
# Dependencies
######################################################################
library(plyranges)

######################################################################
# Fixed association analaysis parameter sets
######################################################################
HiF.association.strategies.df <- 
    tribble(
        ~HiF.type,            ~association.strategy,
        'DIR.anchors',        'nearby',
        'DIR',                'within',
        'compartment.switch', 'nearby',
        'compartment.region', 'within',
        'loop.nesting',       'within',
        'loop.anchor',        'within',
        'loop',               'within',
        'TAD.Boundary',       'nearby',
        'TAD',                'within'
    )

######################################################################
# Load Gene info + Differential results
######################################################################
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
                'Symbol',
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
    dplyr::rename_with(~ str_replace(.x, '^', 'Target.Gene.')) %>% 
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
                mutate(comparison=str_remove(basename(filepath), suffix)) %>% 
                # filter(!grepl('iPSC', filepath)) %>% 
                # Tidy pairwise metadata
                separate_wider_delim(
                    comparison,
                    delim='_vs_',
                    names=c('Sample.Group.Numerator', 'Sample.Group.Denominator'),
                    cols_remove=FALSE
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
                dplyr::rename('EnsemblID'=ensemblid) %>% 
                dplyr::rename('Symbol'=symbol) %>% 
                dplyr::rename('Gene.Type'=type) %>% 
                select(
                    -c(
                        filepath,
                        comparison,
                        # Sample.Group.Numerator, Sample.Group.Denominator,
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

prep_DESeq2_results_for_associations <- function(force.redo=FALSE){
    load_all_DESeq2_results(force.redo=force.redo) %>% 
    rename_with(
        .,
        .fn=~str_replace(.x, '$', '.DESeq2'),
        .cols=-c(starts_with('Sample.Group.'), Gene.Type, Symbol, EnsemblID, chr, start, end)
    ) %>% 
    rename_with(
        .fn=~str_replace(.x, '^', 'Target.Gene.'),
        # .cols=!ends_with('.DESeq2')
        .cols=c(Symbol, EnsemblID, chr, start, end)
    ) %>% 
    mutate(comparison=glue('{Sample.Group.Numerator} vs {Sample.Group.Denominator}'))
}

######################################################################
# Combine all indirect gene-linked functional loci annotations together 
######################################################################
# eQTL data
get_remote_filepaths_for_eQTLs_of_interest <- function(){
    EQTL_REMOTE_FILEPATHS_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    filter(
        tissue_label %in% 
            c(
                'brain (DLPFC)',
                'brain (cerebellum)',
                'brain (putamen)',
                'brain (substantia nigra)',
                'brain (amygdala)',
                'brain (anterior cingulate cortex)',
                'brain (caudate)',
                'brain (cortex)',
                'brain (hippocampus)',
                'brain (hypothalamus)',
                'brain (nucleus accumbens)',
                'brain (spinal cord)',
                'neocortex',
                'microglia',
                'neuron,',
                'neural progenitor',
                'sensory neuron',
                'serotonergic neuron',
                'tibial nerve',
                'dopaminergic neuron',
                'astrocyte',
                'dendritic cell',
                'plasmacytoid dendritic cell',
                # controls?
                'neutrophil',
                'lung'
            )
    ) %>% 
    filter(str_detect(condition_label, 'naive')) %>% 
    filter(
        quant_method %in% 
            c(
                'ge',
                'tx',
                'txrev'
            )
    )
}

download_and_clean_all_eQTL_files <- function(
    p.adj.thresh=0.1,
    force.redo=FALSE){
    get_remote_filepaths_for_eQTLs_of_interest() %>%
    dplyr::rename('url'=ftp_cs_path) %>% 
    mutate(
        eQTL.results=
            pmap(
                .l=.,
                .f=
                    function(url, p.adj.thresh, force.redo, ...){ 
                        destfile <- file.path(CACHED_EQTLS_DIR, basename(url))
                        if (!file.exists(destfile) & !force.redo) {
                            download.file(
                                url=url,
                                destfile=destfile
                            )
                        } 
                        destfile %>%
                        gzfile() %>% 
                        read_tsv(show_col_types=FALSE) %>%
                        mutate(p.adj=p.adjust(pvalue, method='BH')) %>% 
                        filter(p.adj < p.adj.thresh)

                    },
                p.adj.thresh=p.adj.thresh,
                force.redo=force.redo,
                .progress=TRUE
            )
    ) %>% 
    unnest(eQTL.results) %>% 
    # parse eQTL locus locations into separate columns 
    separate_wider_delim(
        region,
        delim=':',
        names=c('chr', 'region')
    ) %>% 
    separate_wider_delim(
        region,
        delim=fixed('-'),
        names=c('start', 'end', 'tmp'),
        too_few='align_start'
    ) %>% 
    # idk wierd case where start coords are negative?
    filter(is.na(tmp)) %>% select(-c(tmp)) %>% 
    # only keep most significant eQTL at each locus
    add_count(
        gene_id,
        chr, start, end,
        name='n.eQTLs.at.locus'
    ) %>% 
    group_by(
        gene_id,
        chr, start, end
    ) %>% 
    slice_min(p.adj, n=1) %>% 
    ungroup() %>% 
    dplyr::rename(
        'PiP'=pip,
        'csID'=cs_id,
        'cs.size'=cs_size,
        'VariantID'=variant,
        'Target.EnsemblID'=molecular_trait_id,
        'Target.Gene.EnsemblID'=gene_id
    ) %>% 
    select(
        -c(
            beta, se, z, cs_min_r2,
            ends_with('_id'),
            ends_with('_path'), 
            url
        )
    )
}
# ABC data
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
    mutate(is.enhancer.sig=ABC.Score > 0.15) %>% 
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
# combine everything
nest_and_combine_all_functional_loci_dataset <- function(
    gene.linked.functional.annotations,
    gene.annotations.df){
    # gene.linked.functional.annotations=list(TFBS.df, eQTLs.df, enhancers.df); gene.annotations.df=load_gene_annotations()
    pmap(
        .l=list(df=gene.linked.functional.annotations),
        .f=
            function(df, gene.annotations.df, ...){
                df %>% 
                {
                    if (length(setdiff(c('Target.Gene.EnsemblID', 'Target.Gene.Symbol'), colnames(.))) == 0) {
                        left_join(
                            .,
                            gene.annotations.df,
                            relationship='many-to-many',
                            by=
                                join_by(
                                    Target.Gene.Symbol,
                                    Target.Gene.EnsemblID
                                    # Target.Gene.Symbol == Symbol,
                                    # Target.Gene.EnsemblID == EnsemblID
                                )
                        )
                    } else if (length(setdiff(c('Target.Gene.EnsemblID'), colnames(.))) == 0) {
                        left_join(
                            .,
                            gene.annotations.df,
                            relationship='many-to-many',
                            by=join_by(Target.Gene.EnsemblID)
                            # by=join_by(Target.Gene.EnsemblID == EnsemblID)
                                    
                        )
                    } else if (length(setdiff(c('Target.Gene.Symbol'), colnames(.))) == 0) {
                        left_join(
                            .,
                            gene.annotations.df,
                            relationship='many-to-many',
                            by=join_by(Target.Gene.Symbol)
                            # by=join_by(Target.Gene.Symbol == Symbol,)
                        )
                    } else {
                        .
                    }
                }
            },
        gene.annotations.df=gene.annotations.df,
        .progress=TRUE
    ) %>% 
    bind_rows() %>% 
    # nest specific associations by type
    nest(
        associations.df=
            -c(
                association.source,
                association.subtype,
                association.type
            )
    ) %>%
    # convert associations info to granges
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
}

######################################################################
# Map HiFs to associated Genes directly or via linked functional loci 
######################################################################
map_HiF_to_genes_with_associations_within <- function(
    HiFs.df,
    associations.df,
    deg.results.df,
    ...){
    # paste('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()', sep='; ')
    # row.index=250 ; resolution=tmp$resolution[[row.index]]; SampleID=tmp$SampleID[[row.index]]; HiF.type=tmp$HiF.type[[row.index]]; HiFs.df=tmp$HiFs.df[[row.index]]; association.strategy=tmp$association.strategy[[row.index]]; association.type=tmp$association.type[[row.index]]; association.subtype=tmp$association.subtype[[row.index]]; association.source=tmp$association.source[[row.index]]; associations.df=tmp$associations.df[[row.index]]; results_file=tmp$results_file[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% t()
    # map all associations to each HiF they are inside, include all associations outside all HiFs
    associations.df %>% 
    mutate(association.status='Gene-link within HiF') %>%
    # keep all gene-associations even if they overlap 0 HiFs 
    join_overlap_left(HiFs.df) %>% 
    as_tibble() %>% 
    mutate(
        association.status=
            ifelse(
                is.na(association.status),
                'Gene-link outside HiFs',
                association.status
            )
    ) %>% 
    left_join(
        deg.results.df,
        by=select(deg.results.df, starts_with('Target.Gene.')) %>% colnames()
    ) %>% 
    # count(is.na(FeatureID), is.na(association.status), is.na(padj.DESeq2), as.character(comparison))
    dplyr::rename('chr'=seqnames) %>% 
    select(-c(width, strand))
}

map_HiF_to_nearest_associated_gene <- function(
    HiFs.df,
    associations.df,
    ...){
    # paste('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '), 'tmp %>% head(row.index) %>% tail(1) %>% t()', sep='; ')
    # row.index=1; resolution=tmp$resolution[[row.index]]; SampleID=tmp$SampleID[[row.index]]; HiF.type=tmp$HiF.type[[row.index]]; HiFs.df=tmp$HiFs.df[[row.index]]; association.strategy=tmp$association.strategy[[row.index]]; association.subtype=tmp$association.subtype[[row.index]]; association.source=tmp$association.source[[row.index]]; association.type=tmp$association.type[[row.index]]; comparison.DESeq2=tmp$comparison.DESeq2[[row.index]]; associations.df=tmp$associations.df[[row.index]]; results_file=tmp$results_file[[row.index]]; tmp %>% head(row.index) %>% tail(1) %>% t()
    # for each gene-link, get the nearest HiF
    # can perform fisher testing with arbitrary thresholdign later
    join_nearest(
        associations.df,
        HiFs.df,
        distance=TRUE
    ) %>% 
    as_tibble() %>% 
    dplyr::rename('Gene.link.HiF.distance'=distance) %>% 
    # Set coords to be nearest HiF coords, not gene coords
    dplyr::rename('chr'=seqnames) %>% 
    select(-c(width, strand))
}

associate_gene_links_to_HiFs <- function(
    HiFs.df,
    associations.df,
    association.strategy,
    deg.results.df,
    ...){
    if (association.strategy == 'nearby') {
        map_HiF_to_nearest_associated_gene(
            HiFs.df,
            associations.df,
            deg.results.df,
            ...
        )
    } else if (association.strategy == 'within') {
        map_HiF_to_genes_with_associations_within(
            HiFs.df,
            associations.df,
            deg.results.df,
            ...
        )
    } else {
        stop(glue('Invalid association.strategy: {association.strategy}'))
    }
}

make_per_condition_mapping_results_filepath <- function(
    association.type,
    association.subtype,
    association.source,
    association.strategy,
    resolution,
    HiF.type,
    SampleID,
    # comparison.DESeq2,
    ...){
    file.path(
        HIF_GENE_ASSOCIATION_MAPPING_DIR,
        # association metadata
        glue('association.type_{association.type}'),
        glue('association.subtype_{association.subtype}'),
        glue('association.source_{association.source}'),
        glue('association.strategy_{association.strategy}'),
        # HiF metadata
        glue('resolution_{resolution}'),
        glue('HiF.type_{HiF.type}'),
        glue('HiF.scope_per.condition'),
        glue('{SampleID}-HiF.Gene.Associations.tsv')
        # glue('HiF.SampleID_{SampleID}'),
        # HiFs for this SampleID, DEG results from all comparisons
        # glue('{comparison.DESeq2}-HiF.Gene.Associations.tsv')
    )
}

list_all_per_condition_gene_association_HiF_mapping_results <- function(){
    HIF_GENE_ASSOCIATION_MAPPING_DIR %>%
    parse_results_filelist(
        suffix='-HiF.Gene.Associations.tsv',
        filename.column.name='SampleID'
    ) %>%
    mutate(
        gene.HiF.mappings.df=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                show_col_types=FALSE,
                progress=FALSE
            )
    )
}

make_between_condition_mapping_results_filepath <- function(
    association.type,
    association.subtype,
    association.source,
    association.strategy,
    resolution,
    HiF.type,
    Numerator,
    Denominator,
    ...){
    file.path(
        HIF_GENE_ASSOCIATION_MAPPING_DIR,
        # association metadata
        glue('association.type_{association.type}'),
        glue('association.subtype_{association.subtype}'),
        glue('association.source_{association.source}'),
        glue('association.strategy_{association.strategy}'),
        # HiF metadata
        glue('resolution_{resolution}'),
        glue('HiF.type_{HiF.type}'),
        glue('HiF.scope_between.conditions'),
        # per pair of conditions that Hif is differential between
        glue('{Numerator}-{Denominator}-differential.HiF.Gene.Associations.tsv')
    )
}

list_all_bewteen_condition_gene_association_HiF_mapping_results <- function(){
    HIF_GENE_ASSOCIATION_MAPPING_DIR %>%
    parse_results_filelist(
        suffix='-differential.HiF.Gene.Associations.tsv',
        filename.column.name='Comparison'
    ) %>%
    separate_wider_delim(
        Comparison,
        delim='-',
        names=c('SampleID.Numerator', 'SampleID.Denominator')
    ) %>%
    mutate(
        gene.HiF.mappings.df=
            pmap(
                .l=list(filepath),
                .f=read_tsv,
                show_col_types=FALSE,
                progress=FALSE
            )
    )
}

