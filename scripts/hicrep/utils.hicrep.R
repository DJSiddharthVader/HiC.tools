######################################################################
# Dependencies
######################################################################
# library(tidyverse)
# library(magrittr)
library(furrr)
library(glue)

######################################################################
# Load resutls
######################################################################
get_smoothing_param <- function(resolution){
    case_when(
        resolution == 1000000 ~  0,
        resolution == 500000  ~  1,
        resolution == 100000  ~  3,
        resolution == 50000   ~  4,
        resolution == 40000   ~  5,
        resolution == 25000   ~ 10,
        resolution == 10000   ~ 20,
        resolution == 5000    ~ 40
    )
}

load_hicrep_results <- function(
    filepath,
    ...){
    read_tsv(
        filepath,
        skip=2,
        progress=FALSE,
        show_col_types=FALSE,
        col_names=c('SCC')
    )
}

list_all_hicrep_results_files <- function(){
    # List all results files that have been generated
    HICREP_RESULTS_DIR %>% 
    parse_results_filelist(
        suffix='-hicrep.txt',
        filename.column.name='file.pair',
        param_delim='_'
    ) %>%
    # Separate IDs of 2 matrices being compared for each results file
    separate_wider_delim(
        file.pair,
        delim='-',
        names=
            c(
                'SampleID.P1',
                'SampleID.P2'
            )
    ) %>% 
    # standardize sample order
    mutate(
        SampleID.Numerator=
            case_when(
                SampleID.P1 >  SampleID.P2 ~ SampleID.P1,
                SampleID.P1 <  SampleID.P2 ~ SampleID.P2,
                SampleID.P1 == SampleID.P2 ~ SampleID.P1,
            ),
        SampleID.Denominator=
            case_when(
                SampleID.Numerator == SampleID.P1 ~ SampleID.P2,
                SampleID.Numerator == SampleID.P2 ~ SampleID.P1
            )
    ) %>% 
    select(-c(SampleID.P1, SampleID.P2)) %>% 
    # Now extract sample metadata from IDs and 
    extract_all_sample_pair_metadata(
        SampleID.cols=c('SampleID.Numerator', 'SampleID.Denominator'),
        SampleID.fields=c('Edit', 'Celltype', 'Genotype', 'CloneID', 'TechRepID'),
        suffixes=c('Numerator', 'Denominator')
    )
}

load_all_hicrep_results <- function(){
    # Load all files generated from ./scripts/run.hicrep.sh
    list_all_hicrep_results_files() %>% 
    # Load hicrep scores for each comparison
    mutate(
        hicrep.results=
            # future_pmap(
            pmap(
                .l=.,
                .f=load_hicrep_results,
                .progress=TRUE
            )
    ) %>%
    unnest(hicrep.results) %>% 
    select(-c(filepath))
}

post_proces_hicrep_results <- function(results.df){
    results.df %>% 
    mutate(
        isMerged=
            case_when(
                CloneID.Numerator == 'Merged'     & CloneID.Denominator == 'Merged'     ~ 'Merged vs Merged',
                CloneID.Numerator == 'Individual' & CloneID.Denominator == 'Merged'     ~ 'Individual vs Merged',
                CloneID.Numerator == 'Merged'     & CloneID.Denominator == 'Individual' ~ 'Individual vs Merged',
                CloneID.Numerator == 'Individual' & CloneID.Denominator == 'Individual' ~ 'Individual vs Individual',
                TRUE ~ NA
            ) %>%
            factor(
                levels=
                    c(
                        'Merged vs Merged',
                        'Individual vs Merged',
                        'Individual vs Individual'
                    )
            )
    ) %>% 
    filter(isMerged != 'Individual vs Merged') %>%
    mutate(
        window.size=scale_numbers(window.size, force_chr=TRUE),
        is.downsampled.fct=
            ifelse(is.downsampled, 'Downsampled', 'Original') %>%
            factor(levels=c('Downsampled', 'Original')),
        within.Edit=    ifelse(Edit.Numerator     == Edit.Denominator,     'Same Edit',     'Different Edit'),
        within.Celltype=ifelse(Celltype.Numerator == Celltype.Denominator, 'Same Celltype', 'Different Celltype'),
        within.Genotype=ifelse(Genotype.Numerator == Genotype.Denominator, 'Same Genotype', 'Different Genotype'),
        Genotype.Group=
            case_when(
                Genotype.Numerator == 'WT' & Genotype.Denominator == 'WT' ~ 'WT vs WT',
                Genotype.Numerator == 'WT' & Genotype.Denominator != 'WT' ~ ' * vs WT',
                Genotype.Numerator != 'WT' & Genotype.Denominator == 'WT' ~ ' * vs WT',
                TRUE                                                      ~ 'Other'
            ) %>%
            factor(levels=c(' * vs WT', 'WT vs WT', 'Other'))
    )
}

######################################################################
# plotting
######################################################################
make_heatmap_plotdf <- function(
    hicrep.results,
    ...){
    hicrep.results %>% 
    group_by(resolution, Sample.ID) %>% 
    summarize(`Mean HiCRep Score`=mean(hicrep.score)) %>%
    ungroup() %>% 
    # Get sample metadata for both samples per pair
    separate_wider_delim(
        Sample.ID,
        delim=' vs ',
        names=c('A.Sample.ID', 'B.Sample.ID')
    ) %>%
    separate_wider_delim(
        c(A.Sample.ID, B.Sample.ID),
        delim=fixed('.'),
        names_sep='.',
        names=c('Edit', 'Genotype', 'SampleNumber', 'Celltype'),
        cols_remove=FALSE
    ) %>%
    rename(
        'A.Sample.ID'=A.Sample.ID.A.Sample.ID,
        'B.Sample.ID'=B.Sample.ID.B.Sample.ID
    ) %>% 
    rename_with(
        .fn=~ str_replace(.x, '(A|B)\\.Sample\\.ID\\.(.*)$', '\\1.\\2'),
        .cols=matches('^(A|B).Sample.ID.')
    )
}

