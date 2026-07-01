######################################################################
# Dependencies
######################################################################
# library(tidyverse)
# library(magrittr)
library(tictoc)
library(glue)
library(optparse)
library(future)
library(furrr)
# library(HiCExperiment)
library(plyranges)
# library(hictkR)
# for parsing from .cool/.mcool file names
ALL.METADATA.FIELDS <- 
    c(
        'Edit',
        'Celltype',
        'Genotype',
        'CloneID',
        'TechRepID',
        'RefGenome',
        'ReadFilter',
        'MinResolution'
        
    )

######################################################################
# Pairsing and Caching
######################################################################
check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    show_col_types=FALSE,
    silence=FALSE,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        if (!silence) { message("No results file, will just return data") }
        results <- results_fnc(...)
        return_data <- TRUE
    } else {
        # Set read/write functions based on filetype
        output.filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        if (!silence) { message(output.filetype) }
        if (output.filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output.filetype %in% c('.txt', '.tsv')) {
            load_fnc <- partial(read_tsv, show_col_types=show_col_types)
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {output.filetype}'))
        }
        if (file.exists(results_file) & !force_redo) {
            if (!silence) { message(glue('{results_file} exists, not recomputing results')) }
            if (return_data) {
                if (!silence) { message('Loading cached results...') }
                results <- load_fnc(results_file)
            }
        # or force recompute+cache the data and return it
        } else {
            if (file.exists(results_file) & force_redo) {
                if (!silence) { message(glue('{results_file} exists, recomputing results anyways')) }
            } else {
                if (!silence) { message(glue('No cached results, generating: {results_file}')) }
            }
            dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
            results <- results_fnc(...) %T>% save_fnc(results_file)
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
    # Or dont save the results, just return the results
    if (!silence) { toc() }
    if (return_data) {
        return(results)
    } else { 
        return(invisible(NULL))
    }
}

parse_results_filelist <- function(
    input_dir,
    suffix,
    filename.column.name='filename',
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

save_cmds_to_file <- function(
    cmds.df,
    cmds.output.filepath=NULL,
    force_redo=FALSE){
    cmds.df %>% 
    # Only include cmds generating outputfiles that dont exist
    {
        if (!force_redo) {
            filter(., !file.exists(output.filepath))
        } else {
            .
        }
    } %>% 
    {
         if (!is.null(cmds.output.filepath)) {
            select(., cmd) %>% 
            write_tsv(
                cmds.output.filepath,
                col_names=FALSE
            )
            message(glue('Saved all cmds to generat results to: {cmds.output.filepath}'))
         } else {
             .
         }
    }
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

######################################################################
# Format stuff
######################################################################
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

######################################################################
# Parsing Sample Identifiers <-> Sample Metadata
######################################################################
parse_metadata_from_names <- function(
    df,
    info.format,
    info.colname=NULL,
    include_merged_col=FALSE,
    keep_id=TRUE,
    prefix=NULL,
    suffix=NULL,
    delim='.',
    ...) {
    info.colname <- 
        ifelse(
            is.null(info.colname),
            info.format,
            info.colname
        )
    field.names <- 
        case_when(
            info.format  == 'Sample.Group' ~ list(ALL.METADATA.FIELDS[1:3]),
            info.format  == 'SampleID'     ~ list(ALL.METADATA.FIELDS[1:5]), 
            # info.format  == 'MatrixID'     ~ list(ALL.METADATA.FIELDS[c(1:5, 7)])
            info.format  == 'MatrixID'     ~ list(ALL.METADATA.FIELDS)
        ) %>% 
        unlist() %>% 
        {
            if (!is.null(prefix)) {
                paste(prefix, ., sep=delim)
            } else {
                .
            }
        } %>% 
        {
            if (!is.null(suffix)) {
                paste(., suffix, sep=delim)
            } else {
                .
            }
        }
    # Split SampleID into separate metadata columns specified as input
    df %>% 
    {
        if (include_merged_col & info.format != 'Sample.Group') {
            mutate(
                .,
                isMerged=
                    ifelse(
                        grepl('Merged', !!sym(info.colname)),
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
        all_of(info.colname),
        delim=fixed(delim),
        names=field.names,
        cols_remove=!keep_id
    )
}

build_name_from_metadata <- function(
    df,
    info.format,
    info.colname=NULL,
    delim='.',
    sep='',
    in.prefix=NULL,
    in.suffix=NULL,
    out.prefix='',
    out.suffix='',
    ...) {
    # name.col='SampleID'; delim='.'; sep=''; prefix=''; suffix='';
    info.colname <- 
        ifelse(is.null(info.colname), info.format, info.colname)
    name.str <- 
        paste0(out.prefix, info.colname, out.suffix, collapse=delim)
    value.glue.str <- 
        case_when(
            info.format  == 'Sample.Group' ~ list(ALL.METADATA.FIELDS[1:3]),
            info.format  == 'SampleID'     ~ list(ALL.METADATA.FIELDS[1:5]), 
            info.format  == 'MatrixID'     ~ list(ALL.METADATA.FIELDS),
            # info.format  == 'MatrixID'     ~ list(ALL.METADATA.FIELDS[c(1:5, 7)])
            .unmatched='error'
        ) %>%
        unlist() %>% 
        {
            if (!is.null(in.prefix)) {
                paste(in.prefix, ., sep=delim)
            } else {
                .
            }
        } %>% 
        {
            if (!is.null(in.suffix)) {
                paste(., in.suffix, sep=delim)
            } else {
                .
            }
        } %>% 
        paste0( 
            '{', ., '}', 
            collapse=delim
        )
    df %>%
    mutate( "{name.str}" := glue(value.glue.str))
}

convert_SampleID_to_SampleGroup <- function(
    df,
    info.colname=NULL,
    include.metadata=FALSE,
    keep_id=FALSE,
    include_merged_col=FALSE,
    delim='.',
    sep='',
    prefix='',
    suffix='',
    ...) {
    # info.colname=NULL; include.metadata=FALSE; include_merged_col=FALSE; keep_id=FALSE; delim='.'; sep=''; prefix=''; suffix=''
    df %>% 
    parse_metadata_from_names(
        info.format='SampleID',
        info.colname=info.colname,
        include_merged_col=include_merged_col,
        keep_id=keep_id,
        delim=delim,
        prefix='SampleMetadata',
        suffix=NULL
    ) %>% 
    build_name_from_metadata(
        info.format='Sample.Group',
        delim=delim,
        sep=sep,
        in.prefix='SampleMetadata',
        in.suffix=NULL,
        out.prefix=prefix,
        out.suffix=suffix
    ) %>%
    {
        if (!include.metadata) {
            select(., -starts_with('SampleMetadata.'))
        } else {
            dplyr::rename_with(., ~str_remove(.x, '^SampleMetadata.', ''))
        }
    }
}

convert_MatrixID_to_SampleID_and_SampleGroup <- function(
    df,
    info.colname=NULL,
    include.metadata=FALSE,
    keep_id=FALSE,
    delim='.',
    sep='',
    prefix='',
    suffix='',
    ...) {
    # info.colname=NULL; include.metadata=FALSE; include_merged_col=FALSE; keep_id=FALSE; delim='.'; sep=''; prefix=''; suffix=''
    df %>% 
    parse_metadata_from_names(
        info.format='MatrixID',
        info.colname=info.colname,
        include_merged_col=TRUE,
        keep_id=keep_id,
        delim=delim,
        prefix='SampleMetadata',
        suffix=NULL
    ) %>% 
    build_name_from_metadata(
        info.format='SampleID',
        delim=delim,
        sep=sep,
        in.prefix='SampleMetadata',
        in.suffix=NULL,
        out.prefix=prefix,
        out.suffix=suffix
    ) %>% 
    build_name_from_metadata(
        info.format='Sample.Group',
        delim=delim,
        sep=sep,
        in.prefix='SampleMetadata',
        in.suffix=NULL,
        out.prefix=prefix,
        out.suffix=suffix
    ) %>%
    {
        if (!include.metadata) {
            select(., -starts_with('SampleMetadata.'))
        } else {
            dplyr::rename_with(., ~str_remove(.x, '^SampleMetadata.', ''))
        }
    }
}

######################################################################
# Load Specific Data
######################################################################
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

filter_included_samples <- function(df){
    included.samples <- 
        load_sample_metadata() %>%
        filter(Included) %>%
        pull(SampleID)
    df %>% 
    filter(SampleID %in% included.samples | grepl('.Merged.Merged', SampleID))
}

load_chr_sizes <- function(){
    CHROMOSOME_SIZES_FILE %>% 
    read_tsv(
        show_col_types=FALSE,
        col_names=c('chr', 'chr.size.bp')
    ) %>% 
    standardize_data_cols(skip.resolution=FALSE)
}

calculate_cum_chr_sizes <- function(chr.sizes.df){
    # chr.sizes.df <- load_chr_sizes()
    chr.sizes.df %>% 
    mutate(idx=as.integer(chr) - 1) %>%
    mutate(
        cum.size=
            pmap(
                .l=.,
                .f=
                    function(idx, chr.sizes.df, ...){ 
                        chr.sizes.df %>% 
                        head(idx) %>% 
                        summarize(cum.chr.size.bp=sum(chr.size.bp))
                    },
                chr.sizes.df=chr.sizes.df,
                .progress=TRUE
            )
    ) %>%
    unnest(cum.size) %>%
    select(chr, chr.size.bp, cum.chr.size.bp)
}

calculate_position_frac <- function(
    regions.df,
    cum.chr.sizes.df,
    position.cols=c('start', 'end')){
    total.genome.size <- 
        cum.chr.sizes.df %>%
        pull(cum.chr.size.bp) %>%
        max()
    regions.df %>%
    left_join(cum.chr.sizes.df, by=join_by(chr)) %>% 
    mutate(
        across(
            .cols=all_of(position.cols),
            .fns=
                list(
                    'genome.pos'=~ .x + cum.chr.size.bp,
                    'chr.pct'=~ .x / chr.size.bp,
                    'genome.pct'=~ (.x + cum.chr.size.bp) / total.genome.size
                ),
            .names="{.col}.{.fn}"
        )
    ) %>%
    select(-c(chr.size.bp, cum.chr.size.bp))
}

######################################################################
# Reciprocal Genomic Disorder Regions 
######################################################################
load_RGD_regions <- function(as.granges=TRUE){
    RGD_REGIONS_FILE %>%
    read_tsv(show_col_types=FALSE, comment='# ') %>% 
    select(
        RGD,
        Critical_Region.UCSC
    ) %>%
    separate_wider_delim(
        Critical_Region.UCSC,
        delim=':',
        names=c('chr', 'range')
    ) %>%
    separate_wider_delim(
        range,
        delim='-',
        names=c('start', 'end')
    ) %>%
    filter(!is.na(chr)) %>%
    mutate(across(c(start, end), as.integer)) %>% 
    {
        if (as.granges) {
            dplyr::rename(., 'seqnames'=chr) %>% 
            as_granges()
        } else {
            .
        }
    }
}

is_region_colocalized_with_RGD <- function(
    regions.df,
    strategy='within',
    ...){
    rgd.regions.df <- load_RGD_regions(as.granges=TRUE)
    regions.df %>% 
    dplyr::rename('seqnames'=chr) %>% 
    as_granges() %>% 
    {
        if (strategy == 'within') {
            join_overlap_left_within(
                .,
                rgd.regions.df
            )
        } else if (strategy == 'overlaps') {
            join_overlap_left(
                .,
                rgd.regions.df
            )
        } else {
            stop(glue('Invalid strategy for deciding co-localization: {strategy}'))
        }
    } %>%
    as_tibble() %>% 
    dplyr::rename(
        'RGD.status'=RGD,
        'chr'=seqnames
    ) %>% 
    mutate(RGD.binary=ifelse(is.na(RGD.status), NA, 'Within RGD')) %>% 
    select(-c(width, strand))
}

######################################################################
# Calculate MoC comparing two sets of regions (TAD or compartments)
######################################################################
calculate_MoC <- function(
    regions.df.P1,
    regions.df.P2,
    resolution,
    ...){
    # paste('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '), sep='; ')
    # regions.df.P1, regions.df.P2 are both tibbles with 
    # the following 3 columns: chr/seqnames start, end
    # MoC nor to granges objects for finding overlapping regions
    granges.P1 <- 
        regions.df.P1 %>%
        mutate(length=end - start) %>% 
        # mutate(idx=row_number()) %>% 
        { if ('chr' %in% colnames(.)) { dplyr::rename(., 'seqnames'=chr) } else { . } } %>% 
        as_granges()
    granges.P2 <- 
        regions.df.P2 %>%
        mutate(length=end - start) %>% 
        # mutate(idx=row_number()) %>% 
        { if ('chr' %in% colnames(.)) { dplyr::rename(., 'seqnames'=chr) } else { . } } %>% 
        as_granges()
    # granges.P1; granges.P2
    # Now find all overlapping pairs of regions between P1, P2, 
    # and return only the overlaps
    join_overlap_intersect(
        granges.P1,
        granges.P2,
        minoverlap=resolution,
        suffix=c('.P1', '.P2')
    ) %>% 
    as_tibble() %>%
    # mutate(width=width-1) %>% 
    dplyr::rename('chr'=seqnames) %>% 
    # https://link.springer.com/article/10.1186/s13059-018-1596-9#Sec9
    # "Assessment of TAD calller performance"
    # # F_ij^2 / (P_i * Q_j) 
    mutate(width=width-1) %>% 
    mutate(moc.inner=((width**2) / (length.P1 * length.P2))) %>% 
    summarize(
        n.Overlaps=n(),
        MoC=sum(moc.inner) # (-1 + sum_ij moc.inner_ij) / sqrt(|P1| + |P2| - 1)
    ) %>%
    add_column(
        n.regions.P1=nrow(regions.df.P1),
        n.regions.P2=nrow(regions.df.P2)
    ) %>% 
    # mutate(MoC=(MoC - 1) / (sqrt(n.regions.P1 * n.regions.P2) - 1))
    mutate(MoC=(MoC) / (sqrt(n.regions.P1 * n.regions.P2)))
}

calculate_all_MoCs <- function(
    region.comparisons.df,
    ...){
    region.comparisons.df %>% 
    # Finally compute all MoCs for all listed pairs of annotations
    mutate(
        MoCs=
            # pmap(
            future_pmap(
                .l=.,
                .f=calculate_MoC,
                .progress=TRUE
            )
    ) %>%
    select(-c(regions.df.P1, regions.df.P2)) %>%
    unnest(MoCs)
}

######################################################################
# Handle pairs of samples and sample groups
######################################################################
list_all_sample_group_comparisons <- function(
    merge_status='individual',
    comparison.groups=NULL,
    delim='.',
    suffixes=c('Numerator', 'Denominator')) {
    comparison.col.names <- paste('Sample.Group', suffixes, sep=delim)
    # get info + filepaths for all contact matrices
    list_all_mcool_files(merge_status=merge_status) %>%
    mutate(Sample.Group.Copy=Sample.Group) %>% 
    # select(Sample.Group, filepath) %>% 
    nest(samples.df=-c(isMerged, Sample.Group.Copy)) %>% 
    dplyr::rename('Sample.Group'=Sample.Group.Copy) %>%
    full_join(
        .,
        .,
        by=join_by(isMerged),
        suffix=paste0(delim, suffixes)
    ) %>% 
    {
        if (!is.null(comparison.groups)) {
            inner_join(
                .,
                comparison.groups %>% set_names(comparison.col.names),
                by=comparison.col.names
            )
        } else {
            .
        }
    }
}

set_foldchange_direction_as_factor <- function(
    results.df,
    sample_group_priority_fnc,
    group1_colname='Sample.Group.Numerator', 
    group2_colname='Sample.Group.Denominator',
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
        group1_colname :=
            case_when(
                .data[[group1_priority]] <  .data[[group2_priority]] ~ .data[[group1_colname]],
                .data[[group1_priority]] >  .data[[group2_priority]] ~ .data[[group2_colname]],
                .data[[group1_priority]] == .data[[group2_priority]] ~ sort(c(.data[[group1_colname]], .data[[group2_colname]]))[[1]],
                TRUE                                                 ~ NA
            ),
        group2_colname :=
            case_when(
                .data[[group1_priority]] <  .data[[group2_priority]] ~ .data[[group2_colname]],
                .data[[group1_priority]] >  .data[[group2_priority]] ~ .data[[group1_colname]],
                .data[[group1_priority]] == .data[[group2_priority]] ~ sort(c(.data[[group1_colname]], .data[[group2_colname]]))[[1]],
                TRUE                                                 ~ NA
            )
    ) %>% 
    # clean up unneded columns since we now have explicity numerator/denominator labels
    select(-c(ends_with('.Priority')))
}

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
    info.format,
    suffixes,
    info.colnames=NULL,
    delim='.',
    ...){
    # Separate IDs of 2 matrices being compared for each results file
    # Extract sample metadata from IDs
    SampleIDs.df %>% 
    # Split SampleID into separate metadata columns specified as input
    parse_metadata_from_names(
        info.format=info.format,
        info.colname=
            ifelse(
                is.null(info.colnames),
                glue('{info.format}{delim}{suffixes[[1]]}'),
                info.colnames[[1]]
            ),
        suffix=suffixes[[1]],
        delim=delim,
        ...
        # include_merged_col=TRUE,
        # keep_id=TRUE,
    ) %>% 
    parse_metadata_from_names(
        info.format=info.format,
        info.colname=
            ifelse(
                is.null(info.colnames[[2]]),
                glue('{info.format}{delim}{suffixes[[2]]}'),
                info.colnames[[2]]
            ),
        suffix=suffixes[[2]],
        delim=delim,
        ...
        # include_merged_col=TRUE,
        # keep_id=TRUE,
    )
}

tidy_pair_metadata <- function(
    sampleID.pairs.df,
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
        delim=delim,
        suffixes=c('P1', 'P2'),
        ...
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
    # sample.group.comparisons=ALL_SAMPLE_GROUP_COMPARISONS; SampleID.fields=NULL; sampleID_col='Sample.Group'; suffixes=c('.Numerator', '.Denominator'); pair_grouping_cols=c('TAD.method', 'TAD.params', 'TAD.metric', 'resolution', 'chr'); delim='.'; suffixes=c('Numerator', 'Denominator')
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
            extract_sample_pair_metadata(
                data.df=.,
                SampleID.cols=paste0(sampleID_col, delim, suffixes),
                delim=delim,
                suffixes=suffixes
            )
        } else {
            .
        }
    }
}

