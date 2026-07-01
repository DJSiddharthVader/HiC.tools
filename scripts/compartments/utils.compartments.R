######################################################################
# Depdendencies
######################################################################

######################################################################
# Generate Cooltools Results
######################################################################
generate_cooltools_calling_cmd <- function(
    threads,
    normalization,
    resolution,
    MatrixID,
    mcool.filepath,
    track.type,
    phasing.track.filepath,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("track.type_{track.type}"),
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}")
        )
    # Create filepaths
    input.filepath     <- glue("{mcool.filepath}::resolutions/{resolution}")
    compartment.prefix <- glue("{output_dir}/{MatrixID}-")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd       <- glue("mkdir -p {output_dir}")
    compartment.cmd <- glue("cooltools eigs-cis --phasing-track {phasing.track.filepath} --n-eigs 3 {weight_flag}  -o {compartment.prefix} {input.filepath}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=glue("{compartment.prefix}.cis.vecs.tsv"),
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    compartment.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_cooltools_calling_cmds <- function(
    hyper.params.df,
    cmds.output.filepath=NULL,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # list all hyper-params and corresponding input files (i.e. phasing track files)
    # to call compartments with
    GENOME_TRACK_FILES_DIR %>%
    parse_results_filelist(suffix='-genome.track.tsv') %>%
    dplyr::rename('phasing.track.filepath'=filepath) %>% 
    inner_join(
        hyper.params.df,
        by=join_by(track.type, resolution)
    ) %>% 
    # list contacts matrices for all samples to generate compartments for
    cross_join(
        list_all_mcool_files(merge_status=merge_status) %>%
        dplyr::rename('mcool.filepath'=filepath),
    ) %>% 
    mutate(output_dir=file.path(COMPARTMENTS_RESULTS_DIR, glue("Comp.method_{Comp.method}"))) %>% 
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(Comp.method, ...) {
                        case_when(
                            Comp.method == 'cooltools' ~ generate_cooltools_calling_cmd(...),
                            .unmatched='error'
                        )
                    },
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    save_cmds_to_file(
        cmds.output.filepath=cmds.output.filepath,
        force_redo=force_redo
    )
}

######################################################################
# Load cooltools results
######################################################################
quantize_compartment_scores <- function(
    scores.df,
    n.compartment.lvls,
    ...){
    # scores.df=tmp2$scores.df[[1]]; n.compartment.lvls=tmp2$n.compartment.lvls[[1]]
    scores.df %>%
    pivot_longer(
        -c(start, end),
        names_to='score.source',
        values_to='score'
    ) %>% 
    mutate(compartment.category=ifelse(score > 0, 'A', 'B')) %>% 
    group_by(score.source) %>% 
    mutate(
        # bin.abs.zscore=(abs.score - mean(abs.score, na.rm=TRUE)) / sd(abs.score, na.rm=TRUE),
        compartment.strength.lvl=
            cut(
                x=abs(score),
                breaks=n.compartment.lvls,
                labels=seq(1, n.compartment.lvls)
            ) %>%
            as.integer()
    ) %>%
    # ungroup() %>%
    mutate(
        compartment.switch=
            case_when(
                is.na(compartment.category) & is.na(lead(compartment.category)) ~ NA,
                compartment.category == lead(compartment.category) ~ 'no.switch',
                .default=glue('{compartment.category}->{lead(compartment.category, n=1L)}'),
            ),
        bin.label=
            case_when(
                is.na(compartment.category) & is.na(compartment.strength.lvl) ~ NA,
                .default=glue("{compartment.category}.{compartment.strength.lvl}")
            ),
        bin.transition=
            case_when(
                is.na(bin.label) & is.na(lead(bin.label)) ~ NA,
                .default=glue('{bin.label}->{lead(bin.label, n=1L)}'),
            ),
        transition.category=
            case_when(
                compartment.strength.lvl <  lag(compartment.strength.lvl)       ~ 'weaker',
                compartment.strength.lvl >  lag(compartment.strength.lvl)       ~ 'stronger',
                compartment.category     != lead(compartment.category)          ~ 'switch',
                compartment.category == lead(compartment.category)             ~ 'no.switch',
                is.na(compartment.category) | is.na(lead(compartment.category)) ~ NA
            )
    )
}

load_cooltools_compartment_results <- function(
    filepath,
    n.compartment.lvls.list,
    ...){
    # filepath=tmp$filepath[[1]]
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
    dplyr::rename('chr'=chrom) %>% 
    # select(chr, start, end, E1, E2, E3) %>% 
    select(chr, start, end, E1) %>% 
    nest(scores.df=-c(chr)) %>% 
    cross_join(tibble(n.compartment.lvls=n.compartment.lvls.list)) %>% 
    mutate(
        compartment.labels=
            pmap(
                .l=.,
                .f=quantize_compartment_scores
            )
    ) %>% 
    select(-c(scores.df)) %>% 
    unnest(compartment.labels)
}

load_all_cooltools_compartment_results <- function(
    resolutions=NULL,
    n.compartment.lvls.list=c(5)){
    # check_cached_results(
    #     results_file=ALL_COMPARTMENTS_RESULTS_FILE,
    #     force_redo=TRUE,
    #     results_fnc=load_all_cooltools_compartment_results,
    #     resolutions=NULL,
    #     n.compartment.lvls.list=c(5)
    # )
    # force_redo=parsed.args$force.redo; resolutions=parsed.args$resolutions; n.compartment.lvls.list=c(5, 10)
    COMPARTMENTS_RESULTS_DIR %>%
    parse_results_filelist(
        filename.column.name='MatrixID',
        suffix='-.cis.vecs.tsv'
    ) %>%  
    convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
    {
        if (!is.null(resolutions)) {
            filter(., resolution %in% resolutions)
        } else {
            .
        }
    } %>% 
    mutate(
        compartments=
            future_pmap(
                 .l=.,
                 .f=load_cooltools_compartment_results,
                 n.compartment.lvls.list=n.compartment.lvls.list,
                 .progress=TRUE
            )
    ) %>%
    unnest(compartments) %>% 
    select(-c(filepath)) 
}

post_process_cooltools_compartment_results <- function(results.df){
    results.df
}

merge_bins_into_compartments <- function(binwise.df){
    binwise.df %>% 
    dplyr::rename('seqnames'=chr) %>% 
    nest(
        compartments=
            -c(
                resolution, 
                normalization, track.type, score.source,
                SampleID, compartment.category
            )
    ) %>% 
    mutate(
        compartments=
            pmap(
                .l=.,
                .f=
                    function(compartments, ...) {
                        compartments %>% 
                        as_granges() %>% 
                        reduce_ranges(total.score=sum(score)) %>% 
                        as_tibble() %>% 
                        dplyr::rename('chr'=seqnames)
                    }
            )
    ) %>% 
    unnest(compartments) %>% 
    select(-c(strand))
}

######################################################################
# Generate Saddle Plot Data
######################################################################
generate_saddle_data_calculation_cmds <- function(
    normalization,
    resolution,
    contact.type,
    n.bins,
    qrange,
    track.col.name,
    expected.col.name,
    MatrixID,
    mcool.filepath,
    track.filepath,
    expected.path,
    output_dir,
    ...){
    output_dir <- 
        file.path(
            output_dir,
            glue("normalization_{normalization}"),
            glue("resolution_{resolution}"),
            glue("contat.type_{contact.type}"),
            glue("n.bins_{n.bins}"),
            glue("expectation.metric_{expected.col.name}"),
            glue("saddle.metric_{track.col.name}")
        )
    # Create filepaths
    mcool.uri     <- glue("{mcool.filepath}::resolutions/{resolution}")
    track.uri     <- glue("{track.filepath}::{track.col.name}")
    expected.uri  <- glue("{expected.path}::{expected.col.name}")
    output.prefix <- glue("{output_dir}/{MatrixID}-")
    # Compose command to generate TAD for this set of inputs + params
    weight_flag <- 
        case_when(
            normalization == 'balanced' ~ '--clr-weight-name weight',
            normalization == 'raw'      ~ '',
            .unmatched="error"
        )
    mkdir.cmd       <- glue("mkdir -p {output_dir}")
    saddle.cmd <- glue("cooltools saddle --qrange {qrange} --strength {weight_flag} -t {contact.type} --n-bins {n.bins} -o {output.prefix} {mcool.uri} {track.uri} {expected.path}")
    # Paste  all commands together in one line to run in bash
    tibble_row(
        output.filepath=glue("{output.prefix}signals.tsv"),
        cmd=
            paste(
                c(
                    mkdir.cmd,
                    saddle.cmd
                ),
                collapse='; '
            )
    )
}

generate_all_saddle_data_calculation_cmds <- function(
    hyper.params.df,
    merge_status='merged',
    force_redo=FALSE,
    ...){
    # merge_status='merged'; force_redo=FALSE;
    # list all binwise eigenvector results files
    compartment.results.files.df <- 
        COMPARTMENTS_RESULTS_DIR %>%
        parse_results_filelist(
            filename.column.name='MatrixID',
            suffix='.cis.vecs.tsv'
        ) %>% 
        convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
        dplyr::rename('track.filepath'=filepath)
    # list all distance expected contact files
    distance.expectation.results.files.df <- 
        DISTANCE_EXPECTED_CONTACTS_DIR %>% 
        parse_results_filelist(
            filename.column.name='MatrixID',
            suffix='-expected.tsv'
        ) %>% 
        convert_MatrixID_to_SampleID_and_SampleGroup() %>% 
        dplyr::rename('expected.path'=filepath)
    # Contact matrix files
    matrix.files.df <- 
        list_all_mcool_files(merge_status=merge_status) %>%
        dplyr::rename('mcool.filepath'=filepath)
    # Map all inputs together by matching param valuies
    compartment.results.files.df %>% 
    inner_join(
        matrix.files.df,
        by=join_by(isMerged, Sample.Group, SampleID)
    ) %>% 
    inner_join(
        distance.expectation.results.files.df,
        by=join_by(isMerged, Sample.Group, SampleID, normalization, resolution)
    ) %>% 
    # list contacts matrices for all samples to generate compartments for
    # Map any other hyper-params to sets of related files
    inner_join(
        hyper.params.df,
        by=
            join_by(
                normalization,
                resolution,
                contact.type
            )
    ) %>% 
    add_column(output_dir=COMPARTMENT_SADDLE_FILES_DIR) %>% 
    # build commands from relevant params + input files
    mutate(
        cmd.data=
            pmap(
                .l=.,
                .f=
                    function(Comp.method, ...) {
                        case_when(
                            Comp.method == 'cooltools' ~ generate_saddle_data_calculation_cmds(...),
                            .unmatched='error'
                        )
                    },
                .progress=TRUE
            )
    ) %>%
    unnest(cmd.data) %>% 
    # Only include cmds generating outputfiles that dont exist
    {
        if (!force_redo) {
            filter(., !file.exists(output.filepath))
        } else {
            .
        }
    }
}

