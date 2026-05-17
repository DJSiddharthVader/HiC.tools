###################################################
# Dependencies
###################################################
library(hictkR)

###################################################
# Genome Binning
###################################################
generate_genome_binning_cmds <- function(
     resolution,
     ...){
    # bin the genome at a specific resolution
    output_dir <- 
        file.path(
            GENOME_BINS_FILES_DIR,
            glue("resolution_{resolution}")
        )
    bins.filepath   <- glue("{output_dir}/genome.bins.tsv")
    # Compose command to generate TAD for this set of inputs + params
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    bins.cmd <- glue("cooltools genome binnify {CHROMSIZES_FILE} {resolution} >| {bins.filepath}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=bins.filepath,
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    bins.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_genome_binning_cmds <- function(
    resolutions,
    cmds.output.filepath=NULL,
    force_redo=FALSE,
    ...){
    pmap(
        .l=list(resolutions),
        .f=generate_genome_binning_cmds
    ) %>%
    bind_rows() %>% 
    save_cmds_to_file(
        cmds.output.filepath=cmds.output.filepath,
        force_redo=force_redo
    )
}

list_all_genome_bin_files <- function(){
    GENOME_BINS_FILES_DIR %>% 
    parse_results_filelist(suffix='-genome.bins.tsv') %>% 
    dplyr::rename('genomic.bins.filepath'=filepath)
}

###################################################
# Genome Phasing Tracks
###################################################
generate_phasing_track_computation_cmds <- function(
    genomic.bins.filepath,
    track.type,
    resolution,
    ...){
    # Create filepaths
    output_dir <- 
        file.path(
            GENOME_TRACK_FILES_DIR,
            glue("track.type_{track.type}"),
            glue("resolution_{resolution}")
        )
    track.filepath   <- glue("{output_dir}/{GENOME_NAME}-genome.track.tsv")
    # Compose command to generate TAD for this set of inputs + params
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    track.cmd <- glue("cooltools genome {track.type} {genomic.bins.filepath} {GENOME_NAME} >| {track.filepath}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=track.filepath,
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    track.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_phasing_track_computation_cmds <- function(
    track.types,
    cmds.output.filepath=NULL,
    force_redo=FALSE,
    ...){
    list_all_genome_bin_files() %>% 
    cross_join(tibble(track.type=track.types)) %>% 
    pmap(
        .l=.,
        .f=generate_phasing_track_computation_cmds,
    ) %>%
    bind_rows() %>% 
    save_cmds_to_file(
        cmds.output.filepath=cmds.output.filepath,
        force_redo=force_redo
    )
}

list_all_phasing_track_files <- function(){
    GENOME_TRACK_FILES_DIR %>% 
    parse_results_filelist(suffix='-genome.track.tsv') %>%
    dplyr::rename('phasing.track.filepath'=filepath)
}

###################################################
# Generate Genome Marginal Coverage 
###################################################
generate_distance_expectation_calculation_cmds <- function(
    threads,
    resolution,
    normalization,
    contact.type,
    ignore.diags,
    mcool.filepath,
    MatrixID,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("contact.type_{contact.type}"),
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}")
        )
    # Create filepaths
    mcool.uri     <- glue("{mcool.filepath}::resolutions/{resolution}")
    output.filepath <- glue("{output_dir}/{MatrixID}-expected.tsv")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    main.cmd  <- glue("cooltools expected-{contact.type} --smooth --aggregate-smoothed --ignore-diags {ignore.diags} -p {threads} -o {output.filepath} {mcool.uri}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=output.filepath,
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    main.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_distance_expectation_calculation_cmds <- function(
    hyper.params.df,
    cmds.output.filepath=NULL,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # list contacts matrices for all samples to generate compartments for
    list_all_mcool_files(merge_status=merge_status) %>%
    dplyr::rename('mcool.filepath'=filepath) %>% 
    # Map any other hyper-params to sets of related files
    cross_join(hyper.params.df) %>% 
    add_column(output_dir=DISTANCE_EXPECTED_CONTACTS_DIR) %>% 
    # build commands from relevant params + input files
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=generate_distance_expectation_calculation{contact.type} _cmds,
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    save_cmds_to_file(
        cmds.output.filepath=cmds.output.filepath,
        force_redo=force_redo
    )
}

list_all_distance_expectation_files <- function(){
    DISTANCE_EXPECTED_CONTACTS_DIR %>% 
    parse_results_filelist(suffix='-expected.tsv') %>%
    dplyr::rename('distance.expectation.filepath'=filepath)
}

###################################################
# Distance-Expectation Calculation
###################################################
generate_marginal_coverage_calculation_cmds <- function(
    threads,
    resolution,
    normalization,
    contact.type,
    ignore.diags,
    mcool.filepath,
    MatrixID,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}")
        )
    # Create filepaths
    mcool.uri       <- glue("{mcool.filepath}::resolutions/{resolution}")
    output.filepath <- glue("{output_dir}/{MatrixID}-coverage.tsv")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd <- glue("mkdir -p {output_dir}")
    main.cmd  <- glue("cooltools coverage {weight_flag} --nproc {threads} --output {output.filepath} {mcool.uri}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=output.filepath,
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    main.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_marginal_coverage_calculation_cmds <- function(
    hyper.params.df,
    cmds.output.filepath=NULL,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # list contacts matrices for all samples to generate compartments for
    list_all_mcool_files(merge_status=merge_status) %>%
    dplyr::rename('mcool.filepath'=filepath) %>% 
    # Map any other hyper-params to sets of related files
    cross_join(hyper.params.df) %>% 
    add_column(output_dir=MARGINAL_COVERAGE_DIR) %>% 
    # build commands from relevant params + input files
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=generate_marginal_coverage_calculation_cmds,
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    save_cmds_to_file(
        cmds.output.filepath=cmds.output.filepath,
        force_redo=force_redo
    )
}

list_all_marginal_coverage_files <- function(){
    # List input files (generated by cooltools coverage)
    MARGINAL_COVERAGE_DIR %>% 
    parse_results_filelist(
        suffix='-coverage.tsv',
        filename.column.name='MatrixID',
        param_delim='_'
    ) %>% 
    dplyr::rename('marginal.coverage.filepath'=filepath)
}

###################################################
# Load Genome Marginal Coverage
###################################################
load_coverage_data <- function(
    filepath,
    ...){
    filepath %>% 
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
    dplyr::rename_with(~ str_remove(.x, '_(raw|weight)')) %>% 
    dplyr::rename(
        'chr'=chrom,
        'coverage.cis'=cov_cis,
        'coverage.total'=cov_tot,
    ) %>% 
    pivot_longer(
        starts_with('coverage'),
        names_to='metric',
        names_prefix='coverage.',
        values_to='coverage'
    )
}

load_all_coverage_data <- function(resolutions=NULL){
    list_all_coverage_data(resolutions=resolutions) %>% 
    # Only need raw cis coverage
    filter(weight == 'raw') %>% 
    filter(count.type == 'cis') %>% 
    # load each coverage file and count % of bins > 0 and > 1000 contacts
    mutate(
        coverage.data=
            pmap(
                .l=.,
                .f=load_coverage_data,
                .progress=TRUE
            )
    ) %>%
    unnest(coverage.data) %>%
    select(-c(filepath))
}

post_process_coverage_summaries <- function(
    results.df,
    sample.metadata.df){
    results.df %>% 
    filter(
        ReadFilter == 'mapq_30',
        count.type == 'cis',
        weight == 'raw',
        metric %in% c(
            'bins.n.covered',  
            'bins.n.detected', 
            'bins.n.nz',       
            'bins.n.total',    
            'bins.pct.covered',
            'bins.pct.nz',     
            'coverage.mean',   
            'coverage.median', 
            # 'coverage.min',    
            # 'coverage.q25',    
            # 'coverage.q75',    
            'coverage.total'
        )
    ) %>% 
    filter(
        !(weight == 'balanced' & metric == 'bins.n.covered'),
        !(weight == 'balanced' & metric == 'bins.pct.covered') 
    ) %>% 
    left_join(
        sample.metadata.df %>% select(SampleID, FlowcellID),
        by=join_by(SampleID)
    ) %>%
    mutate(FlowcellID=ifelse(is.na(FlowcellID), 'Merged', FlowcellID))
}

###################################################
# Marginal Coverage Summary Stats
###################################################
compute_summary_stats <- function(df){
    df %>% 
    summarize(
        coverage.min=min(coverage),
        coverage.q25=quantile(coverage, 0.25, na.rm=TRUE),
        coverage.median=median(coverage, na.rm=TRUE),
        coverage.mean=mean(coverage, na.rm=TRUE),
        coverage.q75=quantile(coverage, 0.75, na.rm=TRUE),
        coverage.max=max(coverage),
        coverage.total=sum(coverage, na.rm=TRUE),
        bins.n.total=n(),
        bins.n.detected=sum(!is.nan(coverage)),
        bins.n.nz=sum(coverage > 0, na.rm=TRUE),
        # 1000 is harcoded based on Rao et al. 2014 paper
        bins.n.covered=sum(coverage > 1000, na.rm=TRUE) 
    ) %>%
    ungroup() %>%
    mutate(
        bins.pct.nz=bins.n.nz / bins.n.total,
        bins.pct.covered=bins.n.covered / bins.n.total
    )
}

compute_coverage_summary <- function(
    filepath,
    ...){
    # Load bin-wise Interaction Frequency data
    # filepath=tmp$filepath[[1]]
    coverage.df <- 
        filepath %>%
        load_coverage_data() # %>% filter(!is.nan(coverage))
    # compute summary stats of IF across the entire genome
    genome.df <- 
        coverage.df %>% 
        group_by(metric) %>% 
        compute_summary_stats() %>% 
        add_column(chr='Genome.Wide')
    # compute summary stats of IF per chr
    chrs.df <- 
        coverage.df %>% 
        group_by(chr, metric) %>% 
        compute_summary_stats()
    # Bind + melt summaries
    bind_rows(
        chrs.df, 
        genome.df
    ) %>% 
    dplyr::rename('count.type'=metric) %>% 
    pivot_longer(
        matches('^(coverage|bins)\\.'),
        names_to='metric',
        values_to='value'
    )
}

compute_all_coverage_summaries <- function(resolutions=NULL){
    # List input files (generated by cooltools coverage)
    list_all_coverage_data(resolutions=resolutions) %>% 
        # {.} -> tmp; tmp
    # load each coverage file and count % of bins > 0 and > 1000 contacts
    mutate(
        coverage.summary=
            pmap(
                .l=.,
                .f=compute_coverage_summary,
                .progress=TRUE
            )
    ) %>%
    unnest(coverage.summary) %>%
    select(-c(filepath))
}

