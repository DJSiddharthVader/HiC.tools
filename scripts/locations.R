################################################################################
# Global Stuff
################################################################################
SCRIPT_DIR            <- file.path(BASE_DIR, 'scripts')
SAMPLE_METADATA_FILE  <- file.path(BASE_DIR, SAMPLE_METADATA_FILENAME) 
RESULTS_DIR           <- file.path(BASE_DIR, 'results')
MCOOL_DIR             <- file.path(RESULTS_DIR, 'coolers_library')

################################################################################
# Reference files
################################################################################
REF_DIR                <- file.path(BASE_DIR, 'reference.files')
GENE_ANNOTATIONS_FILE  <- file.path(REF_DIR, 'homo_sapiens_GRCh38_92.gene_annotations.strand.info.txt') 
GENOME_BINS_FILES_DIR  <- file.path(REF_DIR, 'genome.bins')
GENOME_TRACK_FILES_DIR <- file.path(REF_DIR, 'genome.tracks')
GENOME_NAME            <- 'hg38'
GRCH38_DIR             <- file.path(BASE_DIR, 'GRCh38.Reference')
CHROMSIZES_FILE        <- file.path(GRCH38_DIR, 'GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes')

################################################################################
# Sample QC Results
################################################################################
SAMPLE_QC_DIR                   <- file.path(RESULTS_DIR, 'sample.QC')
MARGINAL_COVERAGE_DIR           <- file.path(SAMPLE_QC_DIR, 'coverage')
DISTANCE_EXPECTED_CONTACTS_DIR  <- file.path(SAMPLE_QC_DIR, 'expected.coverage')
RESOLTION_COVERAGE_SUMAMRY_FILE <- file.path(SAMPLE_QC_DIR, 'resolution.coverage.summaries.tsv')
MIN_SAMPLE_RESOLUTION_FILE      <- file.path(SAMPLE_QC_DIR, 'minimum.viable.resolutions.tsv')

################################################################################
# HiCRep results
################################################################################
HICREP_DIR          <- file.path(RESULTS_DIR, 'hicrep')
HICREP_RESULTS_DIR  <- file.path(HICREP_DIR, 'results')
HICREP_RESULTS_FILE <- file.path(HICREP_DIR, 'all.HiCRep.Scores.tsv')

################################################################################
# TAD Annotations
################################################################################
TAD_DIR                          <- file.path(RESULTS_DIR, 'TADs')
TAD_RESULTS_DIR                  <- file.path(TAD_DIR, 'results_TADs')
# hiTAD results
HITAD_TAD_RESULTS_DIR            <- file.path(TAD_RESULTS_DIR, 'TAD.method_hiTAD')
HITAD_TAD_RESULTS_FILE           <- file.path(TAD_DIR, 'all.hiTAD.TADs.tsv')
HITAD_SCORE_RESULTS_FILE         <- file.path(TAD_DIR, 'all.hiTAD.ADI.Scores.tsv')
# cooltools results
COOLTOOLS_TAD_RESULTS_DIR        <- file.path(TAD_RESULTS_DIR, 'TAD.method_cooltools')
COOLTOOLS_TAD_RESULTS_FILE       <- file.path(TAD_DIR, 'all.cooltools.TADs.tsv')
COOLTOOLS_SCORE_RESULTS_FILE     <- file.path(TAD_DIR, 'all.cooltools.Insulation.Scores.tsv')
# ConsensusTAD results
CONSENSUSTAD_TAD_RESULTS_DIR     <- file.path(TAD_RESULTS_DIR, 'TAD.method_ConsensusTAD')
CONSENSUSTAD_TAD_RESULTS_FILE    <- file.path(TAD_DIR, 'all.ConsensusTAD.TADs.tsv')
CONSENSUSTAD_SCORE_RESULTS_FILE  <- file.path(TAD_DIR, 'all.ConsensusTAD.Consensus.Scores.tsv')
# Joint results
ALL_TAD_RESULTS_FILE             <- file.path(TAD_DIR, 'all.all.TADs.tsv')
ALL_TAD_SCORES_FILE              <- file.path(TAD_DIR, 'all.all.TAD.scores.tsv')
ALL_TAD_BOUNDARIES_FILE          <- file.path(TAD_DIR, 'all.all.TAD.Boundaries.tsv')
ALL_MOC_FILE                     <- file.path(TAD_DIR, 'all.all.TAD.MoCs.tsv')
# TAD_BED_FILES_DIR                <- file.path(TAD_DIR, 'BED.files')
# TADCompare results
TADCOMPARE_DIR                   <- file.path(TAD_DIR, 'results_TADCompare')
TADCOMPARE_RESULTS_FILE          <- file.path(TAD_DIR, 'all.TADCompare.results.tsv')
TADCOMPARE_COUNTS_RESULTS_FILE   <- file.path(TAD_DIR, 'all.TADCompare.counts.results.tsv')

################################################################################
# Loop Annotations
################################################################################
LOOPS_DIR                               <- file.path(RESULTS_DIR, 'loops')
LOOP_RESULTS_DIR                        <- file.path(LOOPS_DIR, 'results_loops')
ALL_COOLTOOLS_LOOPS_RESULTS_FILE        <- file.path(LOOPS_DIR, 'all.cooltools.loops.tsv')
# FILTERED_COOLTOOLS_LOOPS_RESULTS_FILE   <- file.path(LOOPS_DIR, 'filtered.cooltools.loops.tsv')
# LOOP_BED_FILES_DIR                      <- file.path(LOOPS_DIR, 'BED.files')
# IDR2D Analysis
LOOPS_IDR2D_DIR                         <- file.path(LOOPS_DIR, 'results_IDR2D')
ALL_IDR2D_RESULTS_FILE                  <- file.path(LOOPS_DIR, 'all.cooltools.IDR2D.results.tsv')
# Nesting analysis
ALL_LOOP_VALENCY_RESULTS_FILE           <- file.path(LOOPS_DIR, 'all.cooltools.valency.results.tsv')
ALL_LOOP_NESTING_RESULTS_DIR            <- file.path(LOOPS_DIR, 'results_nesting')
ALL_LOOP_NESTING_RESULTS_FILE           <- file.path(LOOPS_DIR, 'all.cooltools.nesting.results.tsv')

################################################################################
# Compartment Annotations
################################################################################
COMPARTMENTS_DIR              <- file.path(RESULTS_DIR, 'compartments')
COMPARTMENTS_RESULTS_DIR      <- file.path(COMPARTMENTS_DIR, 'results_compartments')
ALL_COMPARTMENTS_RESULTS_FILE <- file.path(COMPARTMENTS_DIR, 'all.cooltools.compartments.tsv')
COMPARTMENT_SADDLE_FILES_DIR  <- file.path(COMPARTMENTS_DIR, 'results_saddle')
COMPARTMENT_BED_FILES_DIR     <- file.path(COMPARTMENTS_DIR, 'BED_FILES')

################################################################################
# Differential Contact results from multiHiCCompare
################################################################################
DIFFERENTIAL_CONTACTS_DIR             <- file.path(RESULTS_DIR, 'DifferentialContacts')
ALL_MULTIHICCOMPARE_RESULTS_FILE      <- file.path(DIFFERENTIAL_CONTACTS_DIR, 'all.multiHiCCompare.results.tsv')
FILTERED_MULTIHICCOMPARE_RESULTS_FILE <- file.path(DIFFERENTIAL_CONTACTS_DIR, 'filtered.multiHiCCompare.results.tsv')
MULTIHICCOMPARE_SIG_RESULTS_FILE      <- file.path(DIFFERENTIAL_CONTACTS_DIR, 'all.multiHiCCompare.n.results.tsv')

################################################################################
# Functional Genomic Element (FGE) Data + Results
################################################################################
FGE_RAW_DIR                  <- file.path(REF_DIR, 'raw.FGE.data')
RAW_TF_SITES_FILE            <- file.path(FGE_RAW_DIR, 'Homo_sapiens_meta_clusters.interval')
FILTERED_TF_SITES_FILE       <- file.path(FGE_RAW_DIR, 'filtered.Homo_sapiens_meta_clusters.interval.tsv')
ENCODE_CCRE_SITES_FILE       <- file.path(FGE_RAW_DIR, 'GRCh38-cCREs.bed')
CTCF_SITE_FILE               <- file.path(FGE_RAW_DIR, 'CTCF.annotations.tsv')
GENE_DESERT_REGIONS_FILE     <- file.path(FGE_RAW_DIR, 'hg38_gene_deserts.tsv')
GENE_CONSTRAINTS_FILE        <- file.path(FGE_RAW_DIR, 'gene.constraints.CNVR.tsv')
GENOME_GTF_FILE              <- file.path(FGE_RAW_DIR, 'gencode.v38.annotation.gtf.gz')
# Processed FGE annotation data to calculate bin-wise singal results
FGE_DIR                      <- file.path(RESULTS_DIR, 'FGE.Association.Testing')
# pre-computed bin-wise FGE signal results dir
FGE_SIGNAL_DIR               <- file.path(FGE_DIR, 'binwise.FGE.signals')  
# association testng results 
FGE_TEST_RESULTS_DIR         <- file.path(FGE_DIR, 'testing.results')  
FGE_FISHER_TEST_RESULTS_FILE <- file.path(FGE_DIR, 'all.fisher.test.results.tsv')
FGE_TTEST_TEST_RESULTS_FILE  <- file.path(FGE_DIR, 'all.ttest.test.results.tsv')
FGE_CORR_TEST_RESULTS_FILE   <- file.path(FGE_DIR, 'all.corr.test.results.tsv')

################################################################################
# Delta Expression Association Testing
################################################################################
# DEG results
DESEQ2_DATA_DIR             <- file.path(RESULTS_DIR, 'RNASeq', 'DESeq2')
DESEQ2_RESULTS_FILE         <- file.path(RESULTS_DIR, 'RNASeq', 'all.DESeq2.results.tsv')
# Genomic data to link genes to Hi-C features via non-coding annotations
GENE_LOCUS_ASSOCIATIONS_DIR <- file.path(REF_DIR, 'raw.gene.locus.association.data')
CACHED_EQTLS_DIR            <- file.path(GENE_LOCUS_ASSOCIATIONS_DIR, 'raw.eQTL.data')
NASSER_ABC_SCORES_FILE      <- file.path(GENE_LOCUS_ASSOCIATIONS_DIR, 'H1_Derived_Neuronal_Progenitor_Cultured_Cells-Roadmap_abc_scores.tsv')
INTERNAL_ABC_SCORES_FILE    <- file.path(GENE_LOCUS_ASSOCIATIONS_DIR, 'gnomAD_SV_v3.vs.ABC_H1_Derived_Neuronal_Progenitor_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz')
EQTL_REMOTE_FILEPATHS_FILE  <- file.path(GENE_LOCUS_ASSOCIATIONS_DIR, 'tabix_ftp_paths.tsv')
# Association testing results
EXPRESSION_ASSOCIATION_DIR              <- file.path(RESULTS_DIR, 'Delta.Expression.Association.Testing')
ALL_CLEAN_GENE_LOCUS_ASSOCIATIONS_FILE  <- file.path(EXPRESSION_ASSOCIATION_DIR, 'all.cleaned .Gene.Locus.Associations.rds')
HIF_GENE_ASSOCIATION_MAPPING_DIR        <- file.path(EXPRESSION_ASSOCIATION_DIR, 'HiF.Gene.associations')
EXPRESSION_ASSOCIATION_TEST_RESULTS_DIR <- file.path(EXPRESSION_ASSOCIATION_DIR, 'association.test.results')

