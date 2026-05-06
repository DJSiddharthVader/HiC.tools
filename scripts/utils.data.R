###############################################################################################################
# Dependencies
###############################################################################################################
# library(tidyverse)
# library(magrittr)
library(tictoc)
library(glue)
library(optparse)
library(future)
# library(HiCExperiment)
# library(plyranges)
# library(hictkR)

###############################################################################################################
# Pairsing and Caching
###############################################################################################################
install_packages <- function(){
    install.packages(
        c(
            'ggplot2',
            'ggpubr',
            'ggh4x',
            'ggridges',
            'GGally',
            'scales',
            'ggpointdensity',
            'viridis',
            'ggstatsplot',
            'tidyverse',
            'stringi',
            'glue',
            'here',
            'magrittr',
            'dplyr',
            'optparse',
            'cluster',
            'furrr',
            'future',
            'tictoc',
            'devtools',
            'PRIMME',
            'Matrix',
            'HiCcompare',
            'BiocManager'
        )
    )
    BiocManager::install(c('TADCompare'))
    BiocManager::install("multiHiCcompare")
    pak::pak("jasonwong-lab/gghic")
    pak::pak("paulsengroup/hictkR")
    }

check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    show_col_types=FALSE,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        message("No results file, will just return data")
        results <- results_fnc(...)
        return_data <- TRUE
    } else {
        # Set read/write functions based on filetype
        output_filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        message(output_filetype)
        if (output_filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output_filetype %in% c('.txt', '.tsv')) {
            load_fnc <- partial(read_tsv, show_col_types=show_col_types)
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {output_filetype}'))
        }
        if (file.exists(results_file) & !force_redo) {
            message(glue('{results_file} exists, not recomputing results'))
            if (return_data) {
                message('Loading cached results...')
                results <- load_fnc(results_file)
            }
        # or force recompute+cache the data and return it
        } else {
            if (file.exists(results_file) & force_redo) {
                message(glue('{results_file} exists, recomputing results anyways'))
            } else {
                message(glue('No cached results, generating: {results_file}'))
            }
            dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
            results <- results_fnc(...) %T>% save_fnc(results_file)
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
    # Or dont save the results, just return the results
    toc()
    if (return_data) {
        return(results)
    } else { 
        return(invisible(NULL))
    }
}

parse_results_filelist <- function(
    input_dir,
    suffix,
    filename.column.name='MatrixID',
    pattern=NA,
    param_delim='_',
    parse_filepath_to_columns=TRUE,
    ...){
    # input_dir=file.path(TAD_DIR, 'method_cooltools'); suffix='-TAD.tsv'; filename.column.name='MatrixID'; pattern=NA; param_delim='_';
    # !!NOTICE!!
    # This will break if any parameter_dir name has a param_delim character in the name or value, not as the delimiter
    # This shouldnt  break if the filename has a single param_delim character in it 
    # i.e. min_resolution-0.45 is fine when param_delim='-'
    #      min_resolution_0.45 is will break with any value for param_delim
    # input_dir=HICREP_DIR; suffix='-hicrep.txt'; filename.column.name='file.pair'; param_delim='_'
    suffix_pattern <- glue('*{suffix}$')
    # List all results files that exist
    input_dir %>% 
    list.files(
        pattern=suffix_pattern,
        recursive=TRUE,
        full.names=FALSE
    ) %>% 
    tibble(fileinfo=.) %>%
    mutate(filepath=file.path(input_dir, fileinfo)) %>% 
    {
        if (is.na(pattern)) {
            .
        } else {
            filter(., grepl(pattern, filepath))
        }
    } %>% 
    # Extract param info from directory names into structured columns
    {
        if (parse_filepath_to_columns) {
            separate_longer_delim(
                .,
                fileinfo,
                delim='/'
            ) %>%
            mutate(
                fileinfo=
                    ifelse(
                        grepl(suffix_pattern, fileinfo),
                        paste(
                            filename.column.name,
                            fileinfo,
                            sep=param_delim
                        ) %>% 
                        str_remove(suffix),
                        fileinfo
                    )
            ) %>% 
            separate_wider_delim(
                fileinfo,
                delim='_',
                too_many="merge",  # in case filenames have delim inside
                names=
                    c(
                        'Parameter',
                        'Value'
                    )
            ) %>%
            pivot_wider(
                names_from=Parameter,
                values_from=Value
            )
        } else {
            .
        }
    } %>% 
    # Fix column types based on content
    readr::type_convert()
}

handle_CLI_args <- function(
    args=c('threads', 'force', 'resolutions'),
    has.positional=FALSE){
    # parse listed arguments
    parsed.args <- 
        OptionParser() %>%
        {
            if ('threads' %in% args){
                add_option(
                    .,
                    c('-t', '--threads'),
                    type='integer',
                    default=length(availableWorkers()),
                    dest='threads'
                )
            } else {
                .
            }
        } %>% 
        {
            if ('force' %in% args){
                add_option(
                    .,
                    c('-f', '--force'),
                    action='store_true',
                    default=FALSE,
                    dest='force.redo'
                )
            } else {
                .
            }
        } %>% 
        {
            if ('resolutions' %in% args){
                add_option(
                    .,
                    c('-r', '--resolutions'),
                    type='character',
                    # default=paste(c(5, 10, 25, 50, 100) * 1e3, collapse=','),
                    default=paste(c(10, 25, 50, 100) * 1e3, collapse=','),
                    dest='resolutions'
                )
            } else {
                .
            }
        } %>% 
        parse_args(positional_arguments=TRUE)
    # parse list of resolutions if passed
    # if ('resolutions' %in% args){
    if (is.character(parsed.args$options$resolutions)) {
        parsed.args$options$resolutions <- 
            parsed.args$options$resolutions %>%
            str_split(',') %>%
            lapply(as.integer) %>%
            unlist()
    }
    # return positional args if supplied
    if (has.positional){
        parsed.args
    } else {
        parsed.args$options
    }
}

###############################################################################################################
# Format stuff
###############################################################################################################
glue_sample_ID <- function(
    Edit,
    Celltype,
    Genotype,
    CloneID,
    TechRepID,
    prefix='',
    ...) {
    glue('{prefix}{Edit}.{prefix}{Celltype}.{prefix}{Genotype}.{prefix}{CloneID}.{prefix}{TechRepID}')
}

get_info_from_SampleIDs <- function(
    df,
    SampleID.col='SampleID',
    SampleID.delim='.',
    SampleID.fields=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID'
        ),
    field.suffix='',
    include_merged_col=FALSE,
    keep_id=TRUE,
    ...){
    SampleID.fields <- 
        if (field.suffix == ''){ 
            ifelse(
                is.na(SampleID.fields),
                NA,
                SampleID.fields
            )
        } else {
            ifelse(
                is.na(SampleID.fields),
                NA,
                glue('{SampleID.fields}{SampleID.delim}{field.suffix}')
            )
        }
    df %>%
    {
        if (include_merged_col){
            mutate(
                .,
                isMerged=
                    ifelse(
                        grepl('Merged', !!sym(SampleID.col)),
                        'Merged',
                        'Individual'
                    ) %>% 
                    factor()
            )
        } else {
            .
        }
    } %>% 
    # Split SampleID into separate metadata columns specified as input
    separate_wider_delim(
        all_of(SampleID.col),
        delim=fixed(SampleID.delim),
        names=SampleID.fields,
        cols_remove=!keep_id
    )
}

get_info_from_MatrixIDs <- function(
    df,
    MatrixID.col='MatrixID',
    MatrixID.fields=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID',
            NA,
            'ReadFilter',
            NA
        ),
    MatrixID.delim=fixed('.'),
    SampleID.col='SampleID',
    SampleID.fields=
        c(
            'Edit',
            'Celltype',
            'Genotype',
            'CloneID',
            'TechRepID'
        ),
    SampleID.delim=fixed('.'),
    field.suffix='',
    include_merged_col=FALSE,
    keep_id=FALSE,
    ...){
    # MatrixID.col='MatrixID'; MatrixID.fields=c('Edit', 'Celltype', 'Genotype', 'CloneID', 'TechRepID', NA, 'ReadFilter', NA); MatrixID.delim=fixed('.'); SampleID.col='SampleID'; SampleID.fields=c('Edit', 'Celltype', 'Genotype'); SampleID.delim=fixed('.'); field.suffix=''; include_merged_col=FALSE; keep_id=FALSE;
    MatrixID.fields <- 
        if (field.suffix == ''){ 
            ifelse(
                is.na(MatrixID.fields),
                NA,
                MatrixID.fields
            )
        } else {
            ifelse(
                is.na(MatrixID.fields),
                NA,
                glue('{MatrixID.fields}{MatrixID.delim}{field.suffix}')
            )
        }
    SampleID_glue_pattern <- 
        paste0(
            '{', SampleID.fields, '}', 
            collapse=SampleID.delim
        )
    # Set up some column names/variables
    col_names <- c(MatrixID.fields[!is.na(MatrixID.fields)], 'isMerged')
    # extract + format metadata
    df %>%
    {
        if (include_merged_col){
            mutate(
                .,
                isMerged=
                    ifelse(
                        grepl('Merged', !!sym(MatrixID.col)),
                        'Merged',
                        'Individual'
                    ) %>% 
                    factor()
            )
        } else {
            .
        }
    } %>% 
    separate_wider_delim(
        all_of(MatrixID.col),
        delim=MatrixID.delim,
        names=MatrixID.fields,
        too_many='merge',
        cols_remove=!keep_id
    ) %>% 
    # create SampleID
    {
        if (SampleID.col != '') {
            mutate(., !!SampleID.col := glue(SampleID_glue_pattern))
        } else {
            .
        }
    }
}

scale_numbers <- function(
    numbers,
    accuracy=2,
    force_numeric=FALSE,
    force_chr=FALSE){
    if ((is.character(numbers) & force_numeric) | (is.factor(numbers) & force_numeric)) {
        numbers %>%
        as.character() %>% 
        tibble(resolution.str=.) %>%
        mutate(
            suffix=str_extract(resolution.str, '[KMG]b'),
            magnitude=
                case_when(
                    suffix == 'Kb' ~ 1e3,
                    suffix == 'Mb' ~ 1e6,
                    suffix == 'Gb' ~ 1e9,
                    TRUE ~ 1
                ),
            resolution=
                resolution.str %>% 
                str_remove('[KMG]b') %>%
                as.integer() %>%
                multiply_by(magnitude)
        ) %>% 
        pull(resolution) #%>% format(scientific=FALSE) %>% as.numeric()
    } else if (is.numeric(numbers) & force_chr) {
        numbers %>%
        tibble(resolution=.) %>%
        # mutate(resolution=ifelse(resolution == 0, 1, resolution)) %>% 
        mutate(
            magnitude=resolution %>% {log10(. + 1)} %>% floor() %>% {. %/% 3} %>% {. * 3},
            suffix=
                case_when(
                    magnitude == 0 ~ 'bp',
                    magnitude == 3 ~ 'Kb',
                    magnitude == 6 ~ 'Mb',
                    magnitude == 9 ~ 'Gb',
                    TRUE ~ ''
                ),
            resolution.digits=signif(resolution / 10**magnitude, digits=accuracy),
            resolution.str=glue('{resolution.digits}{suffix}')
        ) %>% 
        mutate(resolution.str=fct_reorder(resolution.str, resolution)) %>% 
        pull(resolution.str)
    } else {
        numbers
    }
}

rename_chrs <- function(
    chrs, 
    to_numeric=FALSE,
    to_label=FALSE){
    if (is.factor(chrs)) {
        chrs
    } else if (is.numeric(chrs) & to_numeric){
        chrs
    } else if (is.character(chrs) & to_numeric){
        case_when(
            chrs == 'X'                  ~ 23,
            chrs == 'Y'                  ~ 24,
            chrs %in% as.character(1:22) ~ chrs,
            grepl('chr', chrs)           ~ str_remove(chrs, 'chr'),
            TRUE                         ~ NA
        ) %>%
        as.integer()
    } else if (is.numeric(chrs) & to_label) {
        case_when(
            chrs == 23                   ~ 'X',
            chrs == 24                   ~ 'Y',
            chrs  >  0 & chrs < 23       ~ as.character(chrs),
            TRUE                         ~ NA
        ) %>% 
        as.character() %>% 
        paste0('chr', .) %>%
        factor(levels=CHROMOSOMES)
    } else if (is.character(chrs) & to_label){
        case_when(
            # chrs == 'Genome.Wide'        ~ 'Genome.Wide',
            chrs == 'X'                  ~ 'chrX',
            chrs == 'Y'                  ~ 'chrY',
            chrs %in% as.character(1:22) ~ paste0('chr', chrs),
            grepl('chr', chrs)           ~ chrs,
            TRUE                         ~ NA
        ) %>%
        as.character() %>% 
        factor(levels=c(CHROMOSOMES))
        # factor(levels=c(CHROMOSOMES, 'Genome.Wide'))
    } else {
        stop('Invalid input to rename_chrs()')
    }
}

standardize_data_cols <- function(
    results.df,
    skip.isGenome=FALSE,
    skip.isMerged=FALSE,
    skip.resolution=FALSE,
    skip.chr=FALSE,
    to_numeric=FALSE,
    to_label=TRUE,
    ...){
    results.df %>% 
    {
        if ('chr' %in% colnames(.) & !skip.chr) {
            if ('Genome.Wide' %in% .$chr & !skip.isGenome) {
                mutate(
                    .,
                    chr=
                        chr %>%
                        rename_chrs(to_numeric=to_numeric, to_label=to_label) %>% 
                        factor(levels=c(CHROMOSOMES, 'Genome.Wide'))
                ) %>% 
                mutate(
                    isGenome=
                        ifelse(chr == 'Genome.Wide', 'Genome.Wide', 'Per.Chr') %>% 
                        factor(levels=c('Genome.Wide', 'Per.Chr'))
                )
            } else {
                mutate(
                    .,
                    chr=
                        chr %>% 
                        rename_chrs(to_numeric=to_numeric, to_label=to_label) %>% 
                        factor(levels=CHROMOSOMES)
                )
            } 
        } else {
            .
        }
    } %>% 
    {
        if ('isGenome' %in% colnames(.) & !skip.isGenome) {
            mutate(., isGenome=factor(isGenome, levels=c('Per.Chr', 'Genome.Wide')))
        } else {
            .
        }
    } %>% 
    {
        if ('isMerged' %in% colnames(.) & !skip.isMerged) {
            if (is.logical(results.df$isMerged)) {
                mutate(
                    .,
                    isMerged=
                        ifelse(isMerged, 'Merged', 'Individual') %>%
                        factor(levels=c('Merged', 'Individual'))
                )
            } else {
                mutate(., isMerged=factor(isMerged, levels=c('Merged', 'Individual')))
            }
        } else {
            .
        }
    } %>% 
    {
        if ('resolution' %in% colnames(.) & !skip.resolution) {
            mutate(
                .,
                resolution=
                    resolution %>%
                    scale_numbers(force_numeric=TRUE) %>%
                    scale_numbers(force_chr=TRUE),
            )
        } else {
            .
        }
    }
}

###############################################################################################################
# Parsing Sample Identifiers <-> Sample Metadata
###############################################################################################################
###############################################################################################################
# Load Specific Data
###############################################################################################################
load_sample_metadata <- function(filter=TRUE){
    SAMPLE_METADATA_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    mutate(isMerged=CloneID == 'Merged') %>% 
    {
        if(filter) {
            filter(., Included)
        } else {
            .
        }
    }
}

load_chr_sizes <- function(){
    CHROMOSOME_SIZES_FILE %>% 
    read_tsv(
        show_col_types=FALSE,
        col_names=c('chr', 'chr.size.bp')
    )
}

###############################################################################################################
# Load mcool files
###############################################################################################################
list_included_samples <- function(){
    load_sample_metadata() %>%
    filter(Included) %>%
    pull(SampleID)
}

load_mcool_file <- function(
    filepath,
    resolution=100000,
    normalization="NONE",
    range1="",
    range2="",
    cis=TRUE,
    type='df',
    include_ends=FALSE,
    ...){
    if (normalization %in% c('weight', 'balanced', 'ICE')){
        normalization <- "weight"
    } else {
        normalization <- "NONE"
    }
    if (range2 == "") {
        range2 <- range1
    }
    filepath %>% 
    File(resolution=resolution) %>% 
    fetch(
        range1=range1,
        range2=range2,
        normalization=normalization,
        type=type,
        join=TRUE,
        query_type='UCSC'
    ) %>% 
    # format column names
    {
        if (type == 'df') {
            if (cis) {
                as_tibble(.) %>%
                filter(chrom1 == chrom2) %>% 
                {
                    if (include_ends){
                        dplyr::rename(
                            .,
                            'chr.A'=chrom1,
                            'start.A'=start1,
                            'end.A'=end1,
                            'chr.B'=chrom2,
                            'start.B'=start2,
                            'end.B'=end2,
                            'IF'=count
                        )
                    } else {
                        select(
                            .,
                            -c(
                                chrom2,
                                end1,
                                end2
                            )
                        ) %>% 
                        dplyr::rename(
                            'chr'=chrom1,
                            'range1'=start1,
                            'range2'=start2,
                            'IF'=count
                        )
                    }
                }
            } else {
                as_tibble(.) %>%
                dplyr::rename(
                    'chr1'=chrom1,
                    'chr2'=chrom2,
                    'range1'=start1,
                    'range2'=start2,
                    'IF'=count
                )
            }
        } else {
            as.matrix(.)
        }
    }
}

list_mcool_files <- function(
    pattern='.hg38.mapq_30.1000.mcool',
    only_use_included_samples=TRUE, 
    ...){
    # List all cooler files 
    COOLERS_DIR %>% 
    list.files(
        pattern=pattern,
        recursive=TRUE,
        full.names=TRUE
    ) %>% 
    tibble(filepath=.) %>% 
    # parse sample metadata
    mutate(MatrixID=str_remove(basename(filepath), '.mcool')) %>% 
    get_info_from_MatrixIDs(include_merged_col=TRUE) %>% 
    {
        if (only_use_included_samples){
            filter(., SampleID %in% list_included_samples() | grepl('.Merged.Merged', SampleID))
        } else {
            .
        }
    }
}

load_mcool_files <- function(
    pattern='.hg38.mapq_30.1000.mcool',
    resolutions=NULL,
    normalizations=NULL,
    regions.df=NULL,
    range1s=NULL,
    range2s=NULL,
    progress=TRUE,
    return_metadata_only=FALSE,
    keep_metadata_columns=FALSE,
    ...){
    # Define all genomic regions to load contacts 
    regions.df <- 
        if (is.null(regions.df)) {
            # Get the whole genome
            if ((is.null(range1s)) & (is.null(range2s))) {
                tibble()
                # tibble(
                #     range1=CHROMOSOMES,
                #     range2=CHROMOSOMES
                # )
            # get all contacts within all regions in range1s
            } else if (is.null(range2s)) {
                tibble(
                    range1=range1s,
                    range2=range1s
                )
            # get all contacts between all pairs of regions only (not intra-region contacts)
            } else {
                expand_grid(
                    range1=range1s,
                    range2=range2s
                )
            }
        # Just load the specified regions (1 region per row: chr, start, end)
        } else {
            regions.df
        }
    # List all regions for all samples
    list_mcool_files(
        pattern=pattern,
        resolutions=resolutions,
        normalizations=normalizations,
    ) %>% 
    join_all_rows(regions.df) %>% 
    # Load contacts if specified or just return sample metadata + filepaths + regions
    {
        if (return_metadata_only) {
            .
        } else {
            mutate(
                .,
                contacts=
                    purrr::pmap(
                        .l=.,
                        .f=load_mcool_file,
                        .progress=progress
                    )
            ) %>% 
            select(-c(filepath, range1, range2)) %>% 
            unnest(contacts)
        }
    }
}

###############################################################################################################
# Group samples by condition
###############################################################################################################
set_up_sample_groups <- function(
    sample.groups,
    # resolutions=c(),
    use_merged=FALSE){
    # sample.groups=ALL_SAMPLE_GROUPS; resolutions=parsed.args$resolutions; use_merged=FALSE
    list_mcool_files() %>%
    # get_min_resolution_per_matrix() %>% 
    distinct() %>% 
    {
        if (use_merged) {
            filter(., isMerged)
        } else {
            filter(., !isMerged)
        }
    } %>% 
    # if any group spans multiple edits
    nest(samples.df=-c(isMerged)) %>% 
    cross_join(
        sample.groups %>%
        tibble(Sample.Group=.) %>% 
        mutate(Sample.Group.Pattern=str_replace_all(Sample.Group, 'All', '.*')),
    ) %>%
    # subset relevant samples for each comparison
    rowwise() %>% 
    mutate(
        samples.df=
            samples.df %>%
            mutate(
                Sample.Group=
                    case_when(
                        str_detect(SampleID, Sample.Group.Pattern) ~ Sample.Group,
                        TRUE ~ NA
                    )
            ) %>%
            filter(!is.na(Sample.Group)) %>% 
            list()
    ) %>%
    # # minimum and max resoltion of all individual matrices per comparison
    # mutate(
    #     resolution.min=min(samples.df$resolution),
    #     resolution.max=max(samples.df$resolution)
    # ) %>% 
    ungroup() %>%
    # mutate(resolution=list(unique(c(resolutions, resolution.min, resolution.max)))) %>%
    # unnest(resolution) %>% 
    # mutate(
    #     resolution.type=
    #         case_when(
    #             resolution == resolution.max ~ 'max',
    #             resolution == resolution.min ~ 'min',
    #             TRUE                         ~ NA
    #         )
    # ) %>% 
    # select(
    #     -c(
    #         resolution.min,
    #         resolution.max,
    #         ends_with('.Pattern')
    #     )
    # )
    select(-c(ends_with('.Pattern')))
}

set_up_sample_comparisons <- function(
    comparison.groups,
    merging='individual'){
    # comparison.groups=ALL_SAMPLE_GROUP_COMPARISONS; resolutions=parsed.args$resolutions; merging='individual'
    # get info + filepaths for all contact matrices
    list_mcool_files() %>%
    # get_min_resolution_per_matrix() %>% 
    # distinct() %>% 
    {
        if (merging == 'individual' ) {
            filter(., !isMerged)
        } else if (merging == 'merged' ) {
            filter(., isMerged)
        } else {
            .
        }
    } %>%
    # select(-c(isMerged)) %>% 
    mutate(
        across(
            starts_with('Sample.Group.'),
            ~ str_replace_all(.x, 'All', '.*'),
            .names='{.col}.Pattern'
        )
    ) %>% 
    # nest(samples.df=everything()) %>%
    nest(samples.df=-c(isMerged)) %>%
    # Now group samples by condition,
    cross_join(
        comparison.groups %>% 
        mutate(
            across(
                starts_with('Sample.Group.'),
                ~ str_replace_all(.x, 'All', '.*'),
                .names='{.col}.Pattern'
            )
        )
    ) %>% 
    # subset relevant samples for each comparison
    rowwise() %>% 
    mutate(
        samples.df=
            samples.df %>%
            mutate(
                Sample.Group=
                    case_when(
                        str_detect(SampleID, Sample.Group.Numerator.Pattern) ~ Sample.Group.Numerator,
                        str_detect(SampleID, Sample.Group.Denominator.Pattern) ~ Sample.Group.Denominator,
                        TRUE ~ NA
                    )
            ) %>%
            filter(!is.na(Sample.Group)) %>% 
            list()
    ) %>%
    # # minimum and max resoltion of all individual matrices per comparison
    # mutate(
    #     resolution.min=min(samples.df$resolution),
    #     resolution.max=max(samples.df$resolution)
    # ) %>% 
    # # List every comparison + resolution that is either a min or max for 1 comparison
    # ungroup() %>%
    # mutate(resolution=list(unique(c(resolutions, resolution.min, resolution.max)))) %>%
    # unnest(resolution) %>% 
    # mutate(
    #     resolution.type=
    #         case_when(
    #             resolution == resolution.max ~ 'max',
    #             resolution == resolution.min ~ 'min',
    #             TRUE                         ~ NA
    #         )
    # ) %>% 
    # select(-c(resolution.min, resolution.max)) %>%
    ungroup() %>% 
    select(-c(ends_with('.Pattern')))
}

set_foldchange_direction_as_factor <- function(
    results.df,
    sample_group_priority_fnc,
    group1_colname='Sample.Group.Left', 
    group2_colname='Sample.Group.Right',
    ...){
    # results.df=hicrep.df; sample_group_priority_fnc=SAMPLE_GROUP_PRIORITY_FNC; group1_colname='SampleID.P1'; group2_colname='SampleID.P2'
    # Use this to make sure that when testing is done by edger::exactTest(), which relies only 
    # on factor levels to determine fold-change direction, that the explicitly stated numerator
    group1_priority <- glue('{group1_colname}.Priority') 
    group2_priority <- glue('{group2_colname}.Priority') 
    # set factor levels for all possible sampleID/group values
    # will have the approriate factor level detected by edger
    results.df %>% 
    # Determine sample group priority i.e. 
    # which sample group is numerator in fold-change values (lower priority value) and 
    # which sample group is denominator in fold change values (higher priority value)
    # use priority function to determine which sample group should represent the numerator
    # in the fold change 
    mutate(
        across(
            c(group1_colname, group2_colname),
            ~ sample_group_priority_fnc(.x),
            .names='{.col}.Priority'
        )
    ) %>%
    # Create explicit+consistent Numerator column, i.e. if there is a log(FC) > 0 then
    # the Numerator is enriched for the signal being compared
    # The Denominator is set to be the opposite group
    mutate(
        Sample.Group.Numerator=
            case_when(
                .data[[group1_priority]] <  .data[[group2_priority]] ~ .data[[group1_colname]],
                .data[[group1_priority]] >  .data[[group2_priority]] ~ .data[[group2_colname]],
                .data[[group1_priority]] == .data[[group2_priority]] ~ sort(c(.data[[group1_colname]], .data[[group2_colname]]))[[1]],
                # TRUE ~ paste(list(.data[[group1_priority]], .data[[group2_priority]]), collapse='-')
                TRUE                                                 ~ NA
            ),
        Sample.Group.Denominator=
            case_when(
                Sample.Group.Numerator == .data[[group1_colname]] ~ .data[[group2_colname]],
                Sample.Group.Numerator == .data[[group2_colname]] ~ .data[[group1_colname]],
                # TRUE ~ paste(sort(c(.data[[group1_colname]], .data[[group2_colname]]))[[1]], collapse='-')
                TRUE                                              ~ NA
            )
    ) %>% 
    # clean up unneded columns since we now have explicity numerator/denominator labels
    select(
        -c(
            ends_with('.Priority'),
            group1_colname,
            group2_colname
        )
    )
}

###############################################################################################################
# Handle pairs of samples
###############################################################################################################
join_all_rows <- function(
    df1,
    df2=NULL,
    cols_to_match=c(),
    ...){
    if (is.null(df2)){ 
        df2 <- df1
    }
    if (length(cols_to_match) == 0) {
        cross_join(
            df1,
            df2,
            ...
        ) 
    } else {
        full_join(
            df1,
            df2,
            relationship='many-to-many',
            by=cols_to_match,
            ...
        )
    }
}

get_all_row_combinations <- function(
    df1,
    cols_to_match=c(),
    suffixes=c('.P1', '.P2'),
    keep_self=FALSE,
    ...){
    idx.cols <- paste0('row.idx', suffixes, sep='')
    # Get all combinations of rows with matching attributes (cols_to_pair)
    join_all_rows(
        df1 %>% mutate(row.idx=row_number()),
        df2=NULL,
        cols_to_match=cols_to_match,
        suffix=suffixes,
        ...
    ) %>% 
    {
        if (keep_self){
            .
        } else {
            filter(., !!sym(idx.cols[1]) != !!sym(idx.cols[2]))
        }
    } %>% 
    # remove redundant combinations i.e. A~B vs B~A -> just keep A~B
    rowwise() %>% 
    # mutate(pair.idx=paste(sort(idx.cols), collapse='-')) %>%
    mutate(
        pair.idx=
            sort(c(!!sym(idx.cols[1]), !!sym(idx.cols[2]))) %>% 
            paste(collapse='-')
    ) %>%
    ungroup() %>% 
    distinct(pair.idx, .keep_all=TRUE) %>%
    select(-c(matches('\\.idx')))
}

extract_sample_pair_metadata <- function(
    SampleIDs.df,
    suffixes,
    ...){
    # Separate IDs of 2 matrices being compared for each results file
    # Extract sample metadata from IDs
    SampleIDs.df %>% 
    # Split SampleID into separate metadata columns specified as input
    get_info_from_SampleIDs(
        SampleID.col=colnames(SampleIDs.df)[[1]],
        field.suffix=suffixes[[1]],
        ...
    ) %>% 
    # Split SampleID into separate metadata columns specified as input
    get_info_from_SampleIDs(
        SampleID.col=colnames(SampleIDs.df)[[2]],
        field.suffix=suffixes[[2]],
        ...
    )
}

extract_all_sample_pair_metadata <- function(
    data.df,
    SampleID.fields,
    SampleID.cols=c('SampleID.P1', 'SampleID.P2'),
    suffixes=c('P1', 'P2'),
    ...){
    data.df %>% 
    mutate(
        tidy.metadata=
            extract_sample_pair_metadata(
                SampleIDs.df=
                    select(
                        .data=., 
                        all_of(c(SampleID.cols[[1]], SampleID.cols[[2]]))
                    ),
                suffixes=suffixes,
                SampleID.fields=SampleID.fields,
                keep_id=FALSE,
                ...
            )
    ) %>%
    unnest(tidy.metadata)
}

tidy_pair_metadata <- function(
    sampleID.pairs.df,
    SampleID.fields,
    suffixes,
    delim='.',
    keep_separate_metadata_fields=FALSE,
    ...){
    # sampleID.pairs.df=tmp %>% select(all_of(c('Sample.Group.Numerator', 'Sample.Group.Denominator'))); SampleID.fields=c('Edit', 'Celltype', 'Genotype'); suffixes=c('Numerator', 'Denominator'); delim='.'; keep_separate_metadata_fields=TRUE
    sampleID.pairs.df %>% 
    # necessary to make sure pairs stay unique when pivoting
    mutate(row.index=row_number()) %>% 
    # Separate IDs of 2 matrices being compared for each results file
    extract_sample_pair_metadata(
        keep_id=FALSE,
        SampleID.fields=SampleID.fields,
        delim=delim,
        suffixes=c('P1', 'P2'),
        # ...
    ) %>% 
    # Now pivot to longer possible format i.e. 2 * nrow(sampleID.pairs.df) * length(SampleID.fields) rows
    pivot_longer(
        ends_with(c('P1', 'P2')),
        names_to='metadata.field',
        values_to='value'
    ) %>%
    separate_wider_delim(
        metadata.field,
        delim=delim,
        names=c('metadata.field', 'pair.index')
    ) %>% 
    # Now each there are nrow(sampleID.pairs.df) * length(SampleID.fields) rows
    pivot_wider(
        names_from='pair.index',
        values_from='value'
    ) %>% 
    # sort pair values so that output column is consistent e.g.
    # you will only ever see "DEL vs WT" and not "WT vs DEL", regardless of P1/P2 ordering
    # rowwise() %>% 
    mutate(
        metadata.value=
            case_when(
                P1 >= P2 ~ glue('{P1} vs {P2}'),
                P1 <  P2 ~ glue('{P2} vs {P1}'),
                TRUE     ~ NA
            )
            # sort(c(P1, P2)) %>% 
            # paste(collapse=' vs ')
    ) %>%
    # ungroup() %>% 
    # return to the same number of rows as the input, with separate columns for each field for P1/P2
    pivot_wider(
        names_from=metadata.field,
        # names_glue="{metadata.field}{delim}{.value}",
        names_glue=glue("{metadata.field}[delim]{.value}", .open='[', .close=']'),
        values_from=c(metadata.value, P1, P2)
    ) %>%
    # Clean up columns/column names
    dplyr::rename_with(~ str_remove(.x, '.metadata.value')) %>%
    {
        if (keep_separate_metadata_fields){
            dplyr::rename_with(., ~ str_replace(.x, 'P1$', suffixes[[1]])) %>%
            dplyr::rename_with(   ~ str_replace(.x, 'P2$', suffixes[[2]]))
        } else {
            select(., -ends_with(c('P1', 'P2')))
        }
    } %>% 
    select(-c(row.index))
}

enumerate_pairwise_comparisons <- function(
    data.df,
    sample.group.comparisons=NULL,
    pair_grouping_cols=c(),
    SampleID.fields=NULL,
    sampleID_col='SampleID',
    delim='.',
    suffixes=c('P1', 'P2'),
    ...){
    data.df %>% 
    # get all possible pairs of rows of data.df 
    join_all_rows(
        cols_to_match=pair_grouping_cols,
        suffix=paste0(delim, suffixes)
    ) %>% 
    # only keep explicitly specified comparisons
    {
        if (is.null(sample.group.comparisons)) {
            .
        } else {
            inner_join(
                .,
                sample.group.comparisons,
                by=colnames(sample.group.comparisons)
            )
        }
    } %>% 
    # format sample pair metadata
    {
        if(!is.null(SampleID.fields)) {
            extract_all_sample_pair_metadata(
                data.df=.,
                SampleID.fields=SampleID.fields,
                SampleID.cols=paste0(sampleID_col, delim, suffixes),
                SampleID.delim=delim,
                suffix=suffixes
            )
            # mutate(
            #     .,
            #     tidy.metadata=
            #         tidy_pair_metadata(
            #             sampleID.pairs.df=
            #                 select(
            #                     .data=., 
            #                     all_of(c(sampleID_col.P1, sampleID_col.P2))
            #                 ),
            #             SampleID.fields=SampleID.fields,
            #             suffixes=suffixes,
            #             delim=delim,
            #             ...
            #         )
            # ) %>%
            # unnest(tidy.metadata)
        } else {
            .
        }
    }
}

###################################################
# Boundary-based analyes
###################################################
define_boundary_pairs <- function(
    boundaries.P1,
    boundaries.P2,
    only_keep_overlaps=TRUE,
    ...){
    # boundaries.P1, boundaries.P2 are both tibbles with 
    # the following 3 columns: start, end, length
    # add indices to track all pairs of boundaries
    # only keep overlapping pairs
    inner_join(
        boundaries.P1 %>% mutate(idx=row_number()),
        boundaries.P2 %>% mutate(idx=row_number()),
        suffix=c('.P1', '.P2'),
        by=join_by(overlaps(x$start, x$end, y$start, y$end))
    ) %>% 
    # cross_join(
    #     boundaries.P1 %>% mutate(idx=row_number()),
    #     boundaries.P2 %>% mutate(idx=row_number()),
    #     suffix=c('.P1', '.P2')
    # ) %>%
    # {
    #     if (only_keep_overlaps){
    #         filter(
    #             .,
    #             dplyr::overlaps(start.P1, end.P1, start.P2, end.P2)
    #             # between(start.P1, start.P2, end.P2) |
    #             # between(end.P1,   start.P2, end.P2) |
    #             # between(start.P2, start.P1, end.P1) |
    #             # between(end.P2,   start.P1, end.P1)
    #         )
    #     } else {
    #         .
    #     }
    # } %>%
    # compute similarity metrics for each pair of boundaries 
    # i.e. (start.P1, end.P2) vs (start.P2, end.P2)
    rowwise() %>% 
    mutate( 
        # calcualte difference in boundary sites
        start.diff=start.P2 - start.P1,
        end.diff=end.P2 - end.P1,
        length.diff=length.P2 - length.P1,
        # calculate regions intersection
        rightmost.start=max(start.P1, start.P2),
        leftmost.end=min(end.P1, end.P2),
        intersection=leftmost.end  - rightmost.start,
        # calculate regions union
        leftmost.start=min(start.P1, start.P2),
        rightmost.end=max(end.P1, end.P2),
        union=rightmost.end - leftmost.start,
        # calcualte difference in boundary overlap
        jaccard=intersection / union,
        moc.inner=((intersection**2) / (length.P1 * length.P2))
    ) %>%
    ungroup() %>% 
    # only keep pairs of boundaries that are the most similar of all pairs
    filter(jaccard > 0) %>% 
    # Only keep the most simialar loop (highest jaccard index) for each loop
    group_by(idx.P1) %>% 
    slice_max(jaccard) %>% 
    # break ties with smallest difference in start/end coordinates
    slice_min(start.diff + end.diff, n=1) %>% 
    ungroup() %>% 
    # Repeat for other side of pairs
    group_by(idx.P2) %>% 
    slice_max(jaccard) %>% 
    slice_min(start.diff + end.diff, n=1) %>% 
    ungroup() %>% 
        # {.} -> dbp3; dbp3
    select(
        -c(
            idx.P1, 
            idx.P2,
            rightmost.start,
            leftmost.end,
            intersection,
            leftmost.start,
            rightmost.end,
            union
        )
    )
}

get_boundary_pairs <- function(
    boundaries.P1,
    boundaries.P2,
    only_keep_overlaps=TRUE,
    ...){
    # For all pairs of boundaries between P1 vs P2 only keeps the most similar pair for each
    # so if |P1| = 10 boundary sets (10 rows) and |P2| = 5
    # then for each set in P1 I only keep 1 set from P2, the "most similar" one 
    # i.e. largest jaccard index
    # MoC normalization constant
    nboundaries.P1 <- nrow(boundaries.P1)
    nboundaries.P2 <- nrow(boundaries.P2)
    norm_const <- 1 / (sqrt(nboundaries.P1 * nboundaries.P2) - 1)
    define_boundary_pairs(
        boundaries.P1,
        boundaries.P2,
        only_keep_overlaps=only_keep_overlaps,
    ) %>%
    add_column(MoC.norm=norm_const)
}

calculate_boundary_similarities <- function(
    boundaries.P1,
    boundaries.P2,
    only_keep_overlaps=TRUE,
    ...){
    # paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row_index]]', collapse='; '); row_index=1
    # resolution=tmp$resolution[[row_index]]; method.P1=tmp$method.P1[[row_index]]; Sample.Group.P1=tmp$Sample.Group.P1[[row_index]]; chr=tmp$chr[[row_index]]; boundaries.P1=tmp$boundaries.P1[[row_index]]; method.P2=tmp$method.P2[[row_index]]; Sample.Group.P2=tmp$Sample.Group.P2[[row_index]]; boundaries.P2=tmp$boundaries.P2[[row_index]]; Edit=tmp$Edit[[row_index]]; Genotype=tmp$Genotype[[row_index]]
    # boundaries.P1=boundaries.P1 %>% rename('length'=TAD.length)
    # boundaries.P2=boundaries.P1
    # MoC normalization constant
    nboundaries.P1 <- nrow(boundaries.P1)
    nboundaries.P2 <- nrow(boundaries.P2)
    norm_const <- 1 / (sqrt(nboundaries.P1 * nboundaries.P2) - 1)
    # Compute MoC 
    moc.df <- 
        # For all pairs of boundaries between P1 vs P2 only keeps the most similar pair for each
        # so if |P1| = 10 boundary sets (10 rows) and |P2| = 5
        # then for each set in P1 I only keep 1 set from P2, the "most similar" one 
        # i.e. largest jaccard index
        define_boundary_pairs(
            boundaries.P1,
            boundaries.P2,
            only_keep_overlaps=only_keep_overlaps
        ) %>% 
        summarize(
            n.TADPairs=n(),
            MoC=(sum(moc.inner) - 1) * norm_const
        )
        # pivot_longer(
        #     everything(),
        #     names_to='tmp',
        #     values_to='value'
        # ) %>%
        # separate_wider_delim(
        #     tmp,
        #     delim='.',
        #     names=c('context', 'stat')
        # )
    moc.df
    # # summarize similarity statistics across those "most similar" boundary pair sets
    # stats.df <- 
    #     pairs.df %>% 
    #     dplyr::rename_with(~ str_remove(.x, '.diff$')) %>% 
    #     select(c(start, end, length, jaccard)) %>% 
    #     pivot_longer(
    #          everything(),
    #          names_to='context',
    #          values_to='value'
    #      ) %>% 
    #     group_by(context) %>%
    #     reframe(broom::tidy(summary(value))) %>% 
    #     # reframe(skimr::skim_without_charts(value)) %>% 
    #     pivot_longer(
    #         -c(context),
    #         names_to='stat',
    #         values_to='value'
    #     )
    # # return tidy stats table
    # bind_rows(
    #     stats.df,
    #     moc.df
    # )
}

calculate_all_boundary_similarities <- function(
    boundaries.df,
    pair_grouping_cols,
    sampleID_col,
    SampleID.fields,
    sample.group.comparisons=NULL,
    only_keep_overlaps=TRUE,
    ...){
    # boundaries.df=TADs.df %>% select(resolution, method, Sample.Group, chr, start, end, TAD.length, TAD.start.score, TAD.end.score, starts_with('TAD.inner.')) %>% nest(boundaries=c(start, end, TAD.length, TAD.start.score, TAD.end.score, starts_with('TAD.inner.'))); only_keep_overlaps=TRUE; sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS %>% dplyr::rename('Sample.Group.P1'=Sample.Group.Numerator, 'Sample.Group.P2'=Sample.Group.Denominator,); SampleID.fields=c('Edit', NA, 'Genotype'); pair_grouping_cols=c('resolution', 'chr'); sampleID_col='Sample.Group'; 
    # boundaries.df Must contain the following columns
    # SampleInfo: nested tibble of Sample attributes
    # 3 columns; start, end, length indicating each set of boundaries
    # pair_grouping_cols: list of cols, only compare pairs of boundaries sets that match
    boundaries.df %>% 
    enumerate_pairwise_comparisons(
        sample.group.comparisons=sample.group.comparisons,
        pair_grouping_cols=pair_grouping_cols,
        SampleID.fields=SampleID.fields,
        sampleID_col=sampleID_col,
        suffixes=c('.P1', '.P2'),
    ) %>% 
        # {.} -> tmp;  tmp
        # tmp %>% filter(Sample.Group.P1 == Sample.Group.P2)
    # Finally compute all MoCs for all listed pairs of annotations
    mutate(
        similarities=
            # pmap(
            future_pmap(
                .l=.,
                .f=calculate_boundary_similarities,
                # .f=get_boundary_pairs,
                only_keep_overlaps=only_keep_overlaps,
                .progress=TRUE
            )
    ) %>%
    unnest(similarities) %>%
    select(-c(boundaries.P1, boundaries.P2))
}

