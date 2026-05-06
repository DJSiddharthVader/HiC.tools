###################################################
# Dependencies
###################################################
library(here)
here::i_am('scripts/functional.enrichment/generate.cCRE.enrichment.cmds.R')
BASE_DIR <- here()
# BASE_DIR <- '/data/talkowski/Samples/cohesin_project/HiC'
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(SCRIPT_DIR, 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'functional.enrichment/utils.enrichment.R'))
    source(file.path(SCRIPT_DIR, 'TADs/utils.TADs.R'))
    source(file.path(SCRIPT_DIR, 'loops/utils.loops.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    source(file.path(SCRIPT_DIR, 'DifferentialContacts/utils.multiHiCCompare.R'))
    library(tidyverse)
    library(magrittr)
})
options(scipen=999)
options(show_col_types=FALSE)
options(future.globals.maxSize=2.5 * (1024 ** 3))
# plan(sequential)
plan(multisession, workers=N_WORKERS_FOR_PARLLELIZATION)

###################################################
# CLI Args
###################################################
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force'),
        has.positional=FALSE
    )
FORCE_REDO <- parsed.args$force.redo

###################################################
# Generate HiC Features BED files
###################################################
generate_all_TAD_bed_files(force_redo=FORCE_REDO)
generate_all_Loop_bed_files(force_redo=FORCE_REDO)
# generate_all_Compartment_bed_files(force_redo=FORCE_REDO)
generate_all_DIR_bed_files(force_redo=FORCE_REDO)

###################################################
# Generate Annotation BED files 
###################################################
# First separate ENCODE cCREs into separate bed files per type
load_encode_ccres() %>% 
    split_cCRE_annotations()
# Also generate CTCF Site bed files from JASPAR 2022 Annotations
load_CTCF_sites() %>%
    split_CTCF_sites()
# calculate bin-wise overlap of CTCF sites + su at multiple resolutions
GENOME_BINS_FILES_DIR %>% 
    parse_results_filelist(suffix='.tsv') %>% 
    select(-c(MatrixID)) %>% 
    cross_join(list_all_annotation_bed_files()) %>% 
    mutate(
        results_file=
            file.path(
                BINWISE_FUNCTIONAL_ENRICHMENT_DIR,
                glue('annotation_{annotation}'),
                glue('resolution_{resolution}'),
                glue('{annotation}-overlaps.per.bin.tsv')
            )
    ) %>% 
    # Generate `bedtools intersect` commands to calculate how many cCRES overlap each bin at each resolution
    make_bedtools_cmds(
        force_redo=FORCE_REDO,
        bed_cmds_file=file.path(FUNCTIONAL_ENRICHMENT_DIR, 'all.binwise.intersect.cmds.txt')
    )

###################################################
# Generate bedtools intersect cmds for cCRE enrichment
###################################################
# First list all dirs of bed files for all annotated feautres
HIC_FEATURES_BED_DIRS <- 
    c(
        TAD_BED_FILES_DIR,
        LOOP_BED_FILES_DIR,
        # COMPARTMENT_BED_FILES_DIR,
        MULTIHICCOMPARE_BED_FILES_DIR
    )
# Generate bedtools intersect cmds to calculate annotation overlaps across HiC features
HIC_FEATURES_BED_DIRS %>% 
    tibble(bed_dir=.) %>% 
    # List all bed files in all feature directories
    mutate(
        bed_files=
            pmap(
                .l=.,
                .f=list_bed_files,
                suffix='.bed'
                # suffix='-.*.bed(pe)?'
            )
    ) %>%
    unnest(bed_files) %>% 
    # Annotation bed files to count overlaps with for each features
    cross_join(list_all_annotation_bed_files()) %>% 
    mutate(
        results_file=
            file.path(
                FUNCTIONAL_ENRICHMENT_DIR,
                basename(bed_dir),
                param_dir,
                glue('annotation_{annotation}'),
                glue('{SampleID}-overlaps.tsv')
            )
    ) %>% 
        # {.} -> feature.bed.files.df; bed_cmds_file=ALL_ENRICHMENT_CMDS_FILE; feature.bed.files.df
    # Generate `bedtools intersect` commands to calculate how many cCRES overlap with each HiC feature
    make_bedtools_cmds(
        force_redo=FORCE_REDO,
        bed_cmds_file=file.path(FUNCTIONAL_ENRICHMENT_DIR, 'all.featurewise.intersect.cmds.txt')
    )
# message(glue('Generated file with all overlap commands: '))
# Now you can run the all the generated commands in ALL_TAD_RESULTS_FILE 
# using `gnu parallel` to generate all the ovelrap results files (1 per line)

