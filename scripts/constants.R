# Dependencies
library(tibble)
library(dplyr)
library(magrittr)
library(glue)
options(dplyr.summarise.inform = FALSE)
# location of scripts 
if (grepl('/home/', BASE_DIR)) {
    SCRIPT_DIR <- file.path(BASE_DIR, '../remote.16p/scripts')
} else {
    SCRIPT_DIR <- here('scripts')
}
# Sample metadata filename
SAMPLE_METADATA_FILENAME <- 'HiC.16p.sample_metadata.tsv'
# Factor levels for variaous metadata categories
ALL_SAMPLE_GROUPS <- 
    c(
        '16p.NSC.DUP',
        '16p.NSC.DEL',
        '16p.NSC.WT',
        '16p.iN.DUP',
        '16p.iN.DEL',
        '16p.iN.WT'
    )
ALL_SAMPLE_GROUP_COMPARISONS <- 
    tribble(
        ~Sample.Group.Numerator, ~Sample.Group.Denominator,
        # '16p.iN.DUP',       '16p.iN.DEL', 
        # '16p.NSC.DUP',      '16p.NSC.DEL',
        # '16p.NSC.DUP',      '16p.iN.DUP',
        # '16p.NSC.DEL',      '16p.iN.DEL',
        '16p.NSC.WT',       '16p.iN.WT',
        '16p.iN.DUP',       '16p.iN.WT',  
        '16p.iN.DEL',       '16p.iN.WT',  
        '16p.NSC.DUP',      '16p.NSC.WT',
        '16p.NSC.DEL',      '16p.NSC.WT'
    )
# functions which determine which sample group will be numerator/denominator in comparisons
SAMPLE_GROUP_PRIORITY_FNC <- function(Sample.Group){
    # FC is determined by the edger::exactTest() function called in multiHiCCompare
    # https://github.com/dozmorovlab/multiHiCcompare/blob/dcfe4aaa8eaef45e203f3d7f806232bb613d2c9b/R/glm.R#L69
    # According to the docs for exactTest()
    # """Note that the first group listed in the pair is the baseline for the comparison—so if the pair is c("A","B") then the comparison is B - A, so genes with positive log-fold change are up-regulated in group B compared with group A (and vice versa for genes with negative log-fold change)."""
    # So with the factor level that comes FIRST is used as the baseline i.e. DENOMINATOR
    # https://www.quantargo.com/help/r/latest/packages/edgeR/NEWS/exactTest
    # so for 16p.NSC.DEL vs 16p.NSC.WT we want to force 16p.NSC.DEL be numerator 
    # therefore we make it the SECOND factor level i.e. have  a larger priority number 
    # i.e. this works when the priority of 16p.NSC.DEL > 16p.NSC.WT
    # so if we create a factor based on this priority the factor hass the following levels
    # Levels: 16p.NSC.WT 16p.NSC.DEL
    # which results in the call exactTest(groups=c(16p.NSC.WT, 16p.NSC.DEL))
    # which is DEL as the FC numerator
    case_when(
        str_detect(Sample.Group, 'NSC.DUP') ~ 1,  # always numerator in FCs
        str_detect(Sample.Group,  'iN.DUP') ~ 2,
        str_detect(Sample.Group, 'NSC.DEL') ~ 3,
        str_detect(Sample.Group,  'iN.DEL') ~ 4,
        str_detect(Sample.Group, 'NSC.WT')  ~ 5,
        str_detect(Sample.Group,  'iN.WT')  ~ 6,
        TRUE                                ~ Inf
    )
}
# have consistent colors for genotypes across figures
GENOTYPE_COLORS <- 
    c(
        'WT'='#b3b3ff',
        'DEL'='#ff0000',
        'DUP'='#0000ff'
    )
CELLTYPE_COLORS <- 
    c(
        'NSC'='#b3b3ff',
        'iN'='#ff0000'
    )
CHROMOSOMES <- 
    c(
        paste0('chr', 1:22), 
        'chrX',
        'chrY'
    )
# Genomic region boundaries to refernce during analysis plotting
# 16p deleteion and telomere are to replicate Fig 4B in this paper:
# https://www.nature.com/articles/s41588-022-01203-y
GENOMIC_REGIONS <- 
    tribble(
        ~region.group, ~region,                    ~region.chr, ~region.start, ~region.end,
        "16p",         "16p11.2 CNV",              "chr16",          29488679,    30188679,
        "16p",         "16p Telomere",             "chr16",                 0,     5149999,
        "16p",         "16p11.2",                  "chr16",          24288679,    30188679,
        "16p",         "16p",                      "chr16",                 0,    36800000,
        "16p",         "16q",                      "chr16",          36800001,    90338345,
        "16p",         "chr16",                    "chr16",                 0,    90338345,
        "10q",         "chr10q CNV",               "chr10",          79791504,    87313680,
        "RGDs",        "1q21",                     "chr1",          146081967,   148779515,
        "RGDs",        "7q11.23_WBS",              "chr7",           73174795,    74833380,
        "RGDs",        "8p23.1",                   "chr8",             523442,     5539949,
        "RGDs",        "15q_AS_PWS_large",         "chr15",          22420444,    28730223,
        "RGDs",        "15q13.3",                  "chr15",          30377807,    32441361,
        "RGDs",        "16p13.11",                 "chr16",          14816373,    16353400,
        "RGDs",        "17q12",                    "chr17",          36222499,    38194356,
        "RGDs",        "17q21.31",                 "chr17",            452200,     1272274,
        "RGDs",        "22q11.21_DGS_VCFS_common", "chr22",          18518837,    21562827,
	) %>%	
    mutate(
        region.UCSC=glue("{region.chr}:{region.start}-{region.end}"),
        region.dist=region.end - region.start
    )
