###################################################
# Dependencis
###################################################
library(stringi)
library(furrr)
# library(idr2d)
# library(plyranges)

###################################################
# Handle JASPAR 2022 CTCF data
###################################################
load_CTCF_sites <- function(force.redo=FALSE){
    # https://dozmorovlab.github.io/CTCF/
    # annotation object
    check_cached_results(
        results_file=CTCF_SITE_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                AnnotationHub() %>% 
                # subset specific set of CTCF sites predicted in humans by JASPAR 2022
                subset(
                    species == 'Homo sapiens' &
                    genome == 'hg38' &
                    dataprovider == 'JASPAR 2022' &
                    preparerclass == "CTCF"
                ) %>%
                {.[['AH104727']]} %>%
                # make tidy tibble
                as_tibble() %>% 
                dplyr::rename(
                    'motif'=name,
                    'length'=width,
                    'chr'=seqnames
                ) %>% 
                select(
                    -c(
                        # db,
                        # db.set,
                        # bioset,
                        # genomic.feature,
                        sequence
                    )
                )
            }
    )
}

split_CTCF_sites <- function(CTCF.df){
    CTCF.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    nest(positions.df=-c(motif)) %>%
    mutate(
        output_dir=
            file.path(
                ANNOTATIONS_BED_DIR,
                glue('annotation_CTCF'),
                glue('annotation.type_{motif}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{motif}_JASPAR.2022.CTCF.sites.bed')
            )
    ) %>%
    pwalk(
        .l=.,
        .f=
            function(positions.df, results_file, ...){
                dir.create(dirname(results_file), showWarnings=FALSE, recursive=TRUE)
                write_tsv(
                    positions.df,
                    results_file,
                    col_names=TRUE
                )
            }
    )
}

count_CTCFs <- function(force.redo=FALSE){
    check_cached_results(
        results_file=CTCF_COUNTS_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                load_CTCF_sites() %>%
                dplyr::rename('annotation.type'=motif) %>% 
                count(
                    chr,
                    annotation.type,
                    name='total.annotations'
                )
            }
    )
}

###################################################
# Handle ENCODE cCRE data
###################################################
load_encode_ccres <- function(force.redo=FALSE){
    check_cached_results(
        results_file=ENCODE_CCRE_ANNOTATIONS_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                # Graphical + explicit definitions of of cCRE types
                # https://screen-v4.wenglab.org/about
                # cCRE.descriptions <- 
                #     tribble(
                #         ~cCRE.type,   ~cCRE.description,
                #         'CA',         'Chromatin Accessible Only',
                #         'CA-CTCF',    'Chromatin Accessible with CTCF',
                #         'CA-H3K4me3', 'Chromatin Accessible with H4K4me3',
                #         'CA-TF',      'Chromatin Accessible with TF',
                #         'PLS',        'Promoter-like',
                #         'TF',         'TF Only',
                #         'dELS',       'Distal enhancer-like'
                #         'pELS',       'Proximal enhancer-like'
                #     )
                ccres.df <- 
                    file.path(ENCODE_CCRE_DIR, 'GRCh38-cCREs.bed') %>%
                    read_tsv(
                        col_names=
                            c(
                                'chr',
                                'start', 'end',
                                # 'cCREID.P1', 'cCREID.P2',
                                'cCREID', 'cCREID.2',
                                'cCRE.type'
                            )
                    ) %>% 
                    select(-c(cCREID.2))
                    # left_join(cCRE.descriptions, by='cCRE.type')
            }
    )
}

split_cCRE_annotations <- function(cCRES.df){
    # Split large file into 1 filer per cCRE type 
    # easier to count per type with bedtools this way (multiple databases)
    cCRES.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    nest(positions.df=-c(cCRE.type)) %>% 
    mutate(
        output_dir=
            file.path(
                ANNOTATIONS_BED_DIR,
                glue('annotation_cCRE'),
                glue('annotation.type_{cCRE.type}')
            ),
        results_file=
            file.path(
                output_dir,
                glue('{cCRE.type}_ENCODE.cCREs.bed')
            )
    ) %>%
    pwalk(
        .l=.,
        .f=
            function(positions.df, results_file, ...){
                dir.create(dirname(results_file), showWarnings=FALSE, recursive=TRUE)
                write_tsv(
                    positions.df,
                    results_file,
                    col_names=TRUE
                )
            }
    )
}

count_cCREs <- function(force.redo=FALSE){
    check_cached_results(
        results_file=ENCODE_CCRE_COUNTS_FILE,
        force_redo=force.redo,
        results_fnc=
            function(){
                load_encode_ccres() %>%
                dplyr::rename('annotation.type'=cCRE.type) %>% 
                count(
                    chr,
                    annotation.type,
                    name='total.annotations'
                )
            }
    )
}

###################################################
# Generate BED Files
###################################################
write_all_bed_files <- function(
    features.df,
    force_redo=FALSE){
    features.df %>% 
    {
        if (!force_redo) {
            filter(., !file.exists(results_file))
        } else{
            .
        }
    } %>% 
    # pmap(
    future_pwalk(
        .l=.,
        .f=
            function(features, results_file, ...){
                dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
                features %>%
                arrange(across(colnames(.)[1:3])) %>% 
                write_tsv(
                    str_replace_all(results_file, ' ' , '.'),
                    col_names=TRUE
                )
            },
        .progress=TRUE
    )
}
# Handle TAD Boundaries
format_TAD_anchors_for_bed_files <- function(TADs.df){
    TADs.df %>% 
    # The end is actually the end of the last bin, transformt so its the start of the last bin inside the TAD
    mutate(end=end - resolution) %>% 
    pivot_TADs_to_boundaries()
}

generate_all_TAD_bed_files <- function(force_redo=FALSE){
    # Load all TAD annotations
    message('Loading TADs...')
    TADs.df <- 
        ALL_TAD_RESULTS_FILE %>%
        read_tsv() %>% 
        mutate(
            param_dir=
                file.path(
                    glue('method_{method}'),
                    glue('TAD.params_{TAD.params}'),
                    glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                    # glue('region_{chr}')
                )
        )
    message('Generating TAD BED files...')
    # Make bed files for the span of each TAD in a nested directory structure for bedtools queries
    TADs.df %>% 
        dplyr::rename('# chr'=chr) %>% 
        nest(
            features=
                c(
                    `# chr`,
                    start, 
                    end,
                    starts_with('TAD.'),
                    -TAD.params
                )
        ) %>% 
        mutate(
            results_file=
                file.path(
                    TAD_BED_FILES_DIR,
                    param_dir,
                    'feature.type_TAD.spans',
                    glue('{Sample.Group}-TADs.bed')
                )
        ) %>% 
        write_all_bed_files(force_redo=force_redo)
    # Make bed files for just the TAD anchors in a nested directory structure for bedtools queries
    TADs.df %>% 
        format_TAD_anchors_for_bed_files() %>% 
        dplyr::rename('# chr'=chr) %>% 
        nest(
            features=
                c(
                    `# chr`,
                    boundary.start, 
                    boundary.end,
                    boundary.side,
                    boundary.score,
                    TAD.bins,
                    TAD.length

                )
        ) %>% 
        mutate(
            results_file=
                file.path(
                    TAD_BED_FILES_DIR,
                    param_dir,
                    'feature.type_TAD.boundaries',
                    glue('{Sample.Group}-TAD.boundaries.bed')
                )
        ) %>% 
        write_all_bed_files(force_redo=force_redo)
    message('Finished generating TAD BED files')
}
# Handle Loops interiors as bed
make_loop_span_bed_files <- function(loops.df) {
    loops.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    mutate(anchor.right=anchor.right + resolution) %>%
    nest(
        features=
            c(
                `# chr`,
                anchor.left,
                anchor.right,
                length,
                count,
                enrichment,
                log10.qval
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'feature.type_loop.spans',
                glue('{SampleID}-loops.bed')
            )
    )
}
# Handle Loops as beddpe
make_loop_bedpe_files <- function(loops.df) {
    loops.df %>% 
    dplyr::rename(
        'anchor.left.start'=anchor.left,
        'anchor.right.start'=anchor.right,
        '# chr.left'=chr
    ) %>% 
    mutate(
        chr.right=`# chr.left`,
        anchor.left.end=anchor.left.start   + resolution,
        anchor.right.end=anchor.right.start + resolution
    ) %>% 
    nest(
        features=
            c(
                `# chr.left`, anchor.left.start,  anchor.left.end,
                   chr.right, anchor.right.start, anchor.right.end,
                length,
                count,
                enrichment,
                log10.qval
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'feature.type_loops',
                glue('{SampleID}-loops.bedpe')
            )
    )
}
# Handle Loop Anchors 
make_loop_anchor_bed_files <- function(loop.valency.df) {
    loop.valency.df %>% 
    mutate(
        param_dir=
            file.path(
                'method_cooltools',
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                # glue('feature.type_{chr}')
            )
    ) %>% 
    # pivot stats so every row represents a unique loop
    pivot_wider(
        names_from=feature,
        names_glue="loops_{feature}_{.value}",
        values_from=c(min, mean, median, max)
    ) %>% 
    # Prepare for writing to bed files
    dplyr::rename(
        'anchor.start'=anchor.position,
        '# chr'=chr
    ) %>% 
    mutate(anchor.end=anchor.start + resolution) %>% 
    # add so bed file can incluide column titles in the header
    nest(
        features=
            c(
                `# chr`,
                anchor.start,
                anchor.end,
                valency,
                starts_with('loops_')
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'feature.type_loop.anchors',
                glue('{SampleID}-loop.anchors.bed')
            )
    )
}
# Handle Loop Nesting
make_loop_nesting_bed_files <- function(loop.nesting.df) {
    loop.nesting.df %>% 
    mutate(
        param_dir=
            file.path(
                'method_cooltools',
                glue('kernel_{kernel}'),
                glue('type_{type}'),
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
            )
    ) %>% 
    # Prepare for writing to bed files
    dplyr::rename('# chr'=chr) %>% 
    mutate(nest.end=nest.end + resolution) %>% 
    # add so bed file can incluide column titles in the header
    nest(
        features=
            c(
                `# chr`,
                nest.start,
                nest.end,
                nesting.lvl,
                starts_with('metric.')
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                LOOP_BED_FILES_DIR,
                param_dir,
                'feature.type_loop.nesting',
                glue('{SampleID}-loop.nesting.bed')
            )
    )
}

generate_all_Loop_bed_files <- function(force_redo=FALSE){
    # Load all loop annotations
    message('Loading loops...')
    loops.df <- 
        ALL_COOLTOOLS_LOOPS_RESULTS_FILE %>%
        read_tsv() %>% 
        post_process_cooltools_dots_results() %>% 
        filter_loop_results() %>% 
        mutate(
            param_dir=
                file.path(
                    'method_cooltools',
                    glue('kernel_{kernel}'),
                    glue('type_{type}'),
                    glue('weight_{weight}'),
                    glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                    # glue('region_{chr}')
                )
        )
    message('Generating all loop BED files...')
    # Make bed files for the span of each loop in a nested directory structure for bedtools queries
    loops.df %>% 
        make_loop_span_bed_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)
    # Make bed files for just the loop anchors in a nested directory structure for bedtools queries
    loops.df %>% 
        make_loop_bedpe_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)
    # Make bed files for just the loop anchors in a nested directory structure for bedtools queries
    # We already generated this data with the loop valency analysis, so we can use that and convert it to structured bed files
    ALL_LOOP_VALENCY_RESULTS_FILE %>%
        read_tsv() %>% 
        post_process_loop_valency_results() %>% 
        make_loop_anchor_bed_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)
    # Also make bed files for the loop nesting regions, where every row represents a
    # contiguous set of bins that are overlapped by the exact same set of loops.
    ALL_LOOP_NESTING_RESULTS_FILE %>%
        read_tsv() %>% 
        make_loop_nesting_bed_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)

    message('Finished generating all loop BED files')
}

generate_all_Compartment_bed_files <- function(force_redo=FALSE){
# Load compartment data
    # compartments.df <- 
}
# Handle DIRs interiors
make_DIR_span_bed_files <- function(DIRs.df){
    DIRs.df %>% 
    dplyr::rename('# chr'=chr) %>% 
    nest(
        features=
            c(
                `# chr`,
                region1,
                region2,
                logFC,
                p.adj.gw,
                log.p.adj.gw
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                MULTIHICCOMPARE_BED_FILES_DIR,
                param_dir,
                'feature.type_DIR.spans',
                glue('{SampleID.Numerator}-{SampleID.Denominator}-DIRs.bed')
            )
    )
}
# Handle DIRs as bedpe
make_DIR_bedpe_files <- function(DIRs.df) {
    DIRs.df %>% 
    dplyr::rename(
        'region1.start'=region1,
        'region2.start'=region2,
        '# chr.left'=chr
    ) %>% 
    mutate(
        chr.right=`# chr.left`,
        region1.end=region1.start + resolution,
        region2.end=region2.start + resolution
    ) %>% 
    nest(
        features=
            c(
                `# chr.left`, region1.start, region1.end,
                   chr.right, region2.start, region2.end,
                logFC,
                p.adj.gw,
                log.p.adj.gw
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                MULTIHICCOMPARE_BED_FILES_DIR,
                param_dir,
                'feature.type_DIRs',
                glue('{SampleID.Numerator}-{SampleID.Denominator}-DIRs.bedpe')
            )
    )
}
# Handle DIR anchors
make_DIR_anchor_bed_files <- function(DIRs.df) {
    DIRs.df %>% 
    dplyr::rename(
        '# chr'=chr,
        'anchor.left'=region1,
        'anchor.right'=region2,
    ) %>% 
    pivot_longer(
        c(anchor.left, anchor.right),
        names_to='anchor.side',
        values_to='anchor.start'
    ) %>%
    mutate(anchor.end=anchor.start + resolution) %>% 
    nest(
        features=
            c(
                `# chr`,
                anchor.start,
                anchor.end,
                anchor.side,
                logFC,
                p.adj.gw,
                log.p.adj.gw
            )
    ) %>% 
    mutate(
        results_file=
            file.path(
                MULTIHICCOMPARE_BED_FILES_DIR,
                param_dir,
                'feature.type_DIR.anchors',
                glue('{SampleID.Numerator}-{SampleID.Denominator}-DIR.anchors.bed')
            )
    )
}

generate_all_DIR_bed_files <- function(force_redo=FALSE){
    # Load all DIRs
    message('Loading DIRs...')
    DIRs.df <- 
        ALL_MULTIHICCOMPARE_RESULTS_FILE %>%
        read_tsv() %>% 
        post_process_multiHiCCompare_results() %>% 
        select(-c(bin.pair.idx, distance.bp)) %>% 
        mutate(
            param_dir=
                file.path(
                    glue('merged_{merged}'),
                    glue('resolution_{scale_numbers(resolution, force_numeric=TRUE)}')
                )
        )
    message('Generating DIR BED files...')
    # Make bed files for the region between DIR anchors in a nested directory structure for bedtools queries. Unlike with TADs or loops, I do not necessairily expect the "inside" of a DIR to show signal unless it overlps with some other features.
    DIRs.df %>% 
        make_DIR_span_bed_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)
    # Make bedpe files, which store both anchors on the same line, to preserve the association
    DIRs.df %>% 
        make_DIR_bedpe_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)
    # Make bed files for just the DIR anchors in a nested directory structure for bedtools queries
    DIRs.df %>% 
        make_DIR_anchor_bed_files() %>% 
        write_all_bed_files(force_redo=FORCE_REDO)
    message('Finished generating DIR BED files')
}

###################################################
# Generate bedtools intersect cmds for enrichment
###################################################
list_bed_files <- function(
    bed_dir,
    suffix='-.*.bed',
    filename.column.name='SampleID',
    ...){
    bed_dir %>% 
    parse_results_filelist(
        suffix=suffix,
        filename.column.name=filename.column.name
    ) %>%
    mutate(
        param_dir=dirname(str_remove(filepath, glue('{bed_dir}/'))),
        file.type=str_extract(filepath, '\\.([^\\.]+)$', group=1),
    ) %>%
    select(
        resolution,
        feature.type,
        file.type,
        SampleID,
        param_dir,
        filepath
    )
}

list_all_annotation_bed_files <- function(){
    ANNOTATIONS_BED_DIR %>%
    parse_results_filelist(suffix='.bed') %>%
    mutate(filepath=str_replace_all(filepath, ' ' , '.')) %>% 
    group_by(annotation) %>% 
    summarize(
        annotation.filepaths=paste(filepath, collapse=' '),
        annotation.types=paste(annotation.type, collapse=' ')
    )
}

make_bedtools_cmds <- function(
    feature.bed.files.df,
    bed_cmds_file,
    force_redo=FALSE){
    # bed_cmds_file=ALL_ENRICHMENT_CMDS_FILE; force_redo=FALSE
    dir.create(dirname(bed_cmds_file), showWarnings=FALSE, recursive=TRUE)
    feature.bed.files.df %>% 
    # {
    #     if (!force_redo) {
    #         filter(., !file.exists(results_file))
    #     } else{
    #         .
    #     }
    # } %>% 
    # generate bedtools intersect command
    # bedtools intersect 
    #     -a features.bed     # all loops called on the same chr at the same resolution
    #     -b annotation.fiels # list of all cCRE database files to intersect against
    #     -wa                 # keep all overalpping region info
    #     -names              # include cCREs names in counting
    rowwise() %>% 
    mutate(
        # remove spaces in filenames
        across(
            c(filepath, results_file),
            ~ str_replace_all(.x, ' ', '.')
        ),
        mkdir.cmd=
            glue('mkdir -p {dirname(results_file)}'),
        features.colnames=
            filepath %>%
            read_tsv(n_max=1, show_col_types=FALSE, progress=FALSE) %>%
            colnames() %>% 
            str_remove('^# +') %>% 
            list(),
        annotation.colnames=
            annotation.filepaths %>%
            str_split(' ') %>% 
            first() %>% 
            read_tsv(n_max=1, show_col_types=FALSE, progress=FALSE) %>%
            colnames() %>% 
            str_remove('^# +') %>% 
            list(),
        bed.colnames=
            list(
                c(
                    features.colnames,
                    'annotation.type',
                    paste0(annotation.colnames, '.annotation'),
                    'overlap.length'
                ) %>%
                paste(collapse=fixed('\\t'))
            ),
        add.column.names.cmd=
            # glue("echo -e '# {bed.colnames}' >| {results_file}"),
            glue("echo -e '{bed.colnames}' >| {results_file}"),
        bedtools.intersect.cmd=
            glue('bedtools intersect -a {filepath} -b {annotation.filepaths} -names {annotation.types} -wo >> {results_file}'),
        bed.cmd=
            glue('{mkdir.cmd}; {add.column.names.cmd}; {bedtools.intersect.cmd}')
    ) %>% 
    select(bed.cmd) %>%
    write_tsv(
        bed_cmds_file,
        col_names=FALSE
    )
}

###################################################
# Load Overlap Results
###################################################
load_and_summarize_CTCF_overlaps <- function(filepath){
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
    mutate(log10_qvalue=-log10(qvalue.annotation)) %>% 
    dplyr::rename('score'=score.annotation) %>% 
    select(-c(overlap.length, qvalue.annotation)) %>% 
    group_by(across(-c('score', 'log10_qvalue', ends_with('.annotation')))) %>%
    summarize(
        n.overlaps=n(),
        across(
            .cols=
                c(
                    score,
                    # pvalue,
                    # qvalue,
                    log10_qvalue
                ),
            .fn=
                list(
                    'min'=min,
                    'q25'=partial(stats::quantile, probs=0.25, na.rm=TRUE),
                    'mean'=mean,
                    'median'=median,
                    'var'=var,
                    'q75'=partial(stats::quantile, probs=0.75, na.rm=TRUE),
                    'max'=max,
                    'total'=sum
                ),
            .names="CTCF.{.col}.{.fn}"
        )
    ) %>%
    ungroup()
}

load_and_summarize_cCRE_overlaps <- function(filepath){
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>% 
    select(-c(overlap.length)) %>% 
    group_by(across(-ends_with('.annotation'))) %>% 
    count(name='n.overlaps')
}

fetch_overlap_file_locations <- function(feature.type){
    # print(feature.type)
    case_when(
        feature.type == 'binwise'                ~ BINWISE_ENRICHMENT_DIR,
        feature.type == 'TAD.spans'              ~ TAD_ENRICHMENTS_DIR,
        feature.type == 'TAD.boundaries'         ~ TAD_ENRICHMENTS_DIR,
        feature.type == 'loop.anchors'           ~ LOOP_ENRICHMENTS_DIR,
        feature.type == 'loop.nesting'           ~ LOOP_ENRICHMENTS_DIR,
        feature.type == 'loop.spans'             ~ LOOP_ENRICHMENTS_DIR,
        feature.type == 'compartment.boundaries' ~ COMPARTMENTS_ENRICHMENTS_DIR,
        feature.type == 'compartments'           ~ COMPARTMENTS_ENRICHMENTS_DIR,
        feature.type == 'DIR.anchors'            ~ MULTIHICCOMPARE_ENRICHMENTS_DIR,
        feature.type == 'DIR.spans'              ~ MULTIHICCOMPARE_ENRICHMENTS_DIR,
        # TRUE                                     ~ stop(glue('Invalid feature.type: {feature.type}'))
    )
}

load_all_overlap_results <- function(
    specific.feature.type,
    specific.annotation=NULL,
    ...){
    # specific.feature.type='TAD.boundaries'; specific.annotation.type='CTCF'; 
    specific.feature.type %>% 
    fetch_overlap_file_locations() %>% 
    parse_results_filelist(
        suffix='.tsv',
        filename.column.name='SampleID',
    ) %>% 
    add_column(feature.type=specific.feature.type) %>% 
    separate_wider_delim(
        SampleID,
        delim='-',
        names=c('SampleID', NA)
    ) %>%
    {
        if (!is.null(specific.annotation)) {
            filter(., annotation == specific.annotation)
        } else {
            .
        }
    } %>% 
        # {.} -> tmp; tmp
        # filepath=tmp$filepath[[1]]
    mutate(
        overlaps=
            # pmap(
            future_pmap(
                .l=.,
                .f=
                    function(filepath, annotation, ...){
                        filepath %>% 
                        {
                            if (annotation == 'CTCF') {
                                load_and_summarize_CTCF_overlaps(filepath=.)
                            } else {
                                load_and_summarize_cCRE_overlaps(filepath=.)
                            }
                        }
                    },
                # ...,
                .progress=TRUE
            )
    ) %>%
    unnest(overlaps) %>% 
    dplyr::rename('chr'=chrom) %>% 
    {
        if (specific.annotation == 'CTCF') {
            dplyr::rename(
                .,
                'motif'=annotation.type,
                'CTCF.overlaps.n'=n.overlaps
            )
        } else if (specific.annotation == 'cCRE') {
            dplyr::rename(
                .,
                'cCRE.type'=annotation.type,
                'cCRE.overlaps.n'=n.overlaps
            )
        } else {
            .
        }
    } %>% 
    # only keep columns with at least 2 unique values
    # i.e. exclude all columsn that are only NA
    # select(where(~ n_distinct(.) > 1)) %>% 
    select(-c(filepath))
}

load_specific_overlap_results <- function(
    force.redo=FALSE,
    feature.type,
    annotation,
    ...){
    # feature.type='TAD.boundaries'; annotation.type='CTCF'; force.redo=FALSE
    check_cached_results(
        results_file=
            file.path(
                COALLATED_ENRICHMENTS_DIR,
                glue('annotation_{annotation}'),
                glue('All-{feature.type}-overlaps.tsv')
            ),
        force_redo=force.redo,
        results_fnc=load_all_overlap_results,
        suffix='-overlaps.tsv',
        specific.feature.type=feature.type,
        specific.annotation=annotation
    )
}

post_process_overlap_results <- function(
    results.df,
    annotations.df){
    results.df %>% 
    {
        if ('SampleID' %in% colnames(.)){
            separate_wider_delim(
                .,
                SampleID,
                delim='.',
                names=c('Edit', 'Celltype', 'Genotype'),
                cols_remove=FALSE
            )
        } else {
            .
        }
    } %>% 
    left_join(
        annotations.df,
        by=join_by(chr, annotation.type == annotation.type)
    ) %>% 
    mutate(
        overlaps.per_Kb=n.overlaps / (length / 1000),
        frac.overlaps.per_Kb=overlaps.per_Kb / total.annotations
    )
}

###################################################
# Post Process Overlap Results
###################################################
pivot_CTCF_enrichment_metrics <- function(
    results.df,
    statlist=NULL){
    results.df %>% 
    # Pivot CTCF summary stats to tidy-format
    {
        if ('n.overlaps' %in% colnames(.)) {
            dplyr::rename(., 'CTCF.overlaps.n'=n.overlaps)
        } else {
            .
        }
    } %>% 
    pivot_longer(
         starts_with('CTCF.'),
         names_prefix='CTCF.',
         names_to='sumstat',
         values_to='value'
    ) %>% 
    separate_wider_delim(
        sumstat,
        delim=fixed('.'),
        names=c('metric', 'stat'),
        cols_remove=FALSE
    ) %>%
    {
        if (!is.null(statlist)) {
            filter(., stat %in% statlist)
        } else {
            .
        }
    }
}

assign_nearest_feature_to_bins <- function(bins.df, features.df){
    bins.df %>% 
    # match bins to closest downstream feature 
    left_join(
        features.df,
        relationship='many-to-many',
        by=join_by(chr, closest(start >= feature.end))
    ) %>% 
    # match bins to closest upstream feature
    left_join(
        features.df,
        relationship='many-to-many',
        suffix=c('.ds', '.us'),
        by=join_by(chr, closest(end <= feature.start))
    ) %>% 
    # check if any bins match feature exactly
    left_join(
        features.df %>%
        rename_with(.cols=-c(chr), .fn=~ str_replace(.x, '$', '.match')) %>% 
        add_column(isExactMatch=TRUE),
        relationship='many-to-many',
        by=join_by(chr, start == feature.start.match, end == feature.end.match)
    ) %>% 
    mutate(isExactMatch=ifelse(!is.na(isExactMatch), isExactMatch, FALSE)) %>% 
    # pick the closer TAD boundary for each bin
    mutate(
        # for all bins within each TAD, get distance to nearest feature (start or end)
        dist.to.feature.ds=ifelse(is.na(feature.start.ds), Inf, abs(feature.start.ds - start)),
        dist.to.feature.us=ifelse(is.na(feature.end.us),   Inf, abs(feature.end.us - end)),
        dist.to.nearest=
            case_when(
                isExactMatch                             ~  0,
                dist.to.feature.ds <= dist.to.feature.us ~ -dist.to.feature.ds,
                dist.to.feature.ds >  dist.to.feature.us ~  dist.to.feature.us
            ),
        nearest.boundary=
            case_when(
                isExactMatch                             ~ 'match',
                dist.to.feature.ds <= dist.to.feature.us ~ 'ds',
                dist.to.feature.ds >  dist.to.feature.us ~ 'us'
            ),
        # Create new column with only the data for the closest feature (upstream or downstream)
        across(
            ends_with('.match'), 
            ~ case_when(
                # if ds bin is closest get *.ds column data
                nearest.boundary == 'match' ~ .x,
                nearest.boundary == 'ds'    ~ get(str_replace(cur_column(), '.match', '.ds')), 
                nearest.boundary == 'us'    ~ get(str_replace(cur_column(), '.match', '.us'))
            ),
            .names="{str_remove(.col, '.match$')}"
        )
    ) %>% 
    # Re-categorize bins "close" to HiC Features as HiC Features for testing
    mutate(
        isHiCFeature=
            case_when(
                abs(dist.to.nearest / resolution) <= HiCFeatureRadius.bins ~ 'Feature',
                abs(dist.to.nearest / resolution) >  HiCFeatureRadius.bins ~ 'Not Feature' 
            ) %>% 
            # force isHiCFeature == TRUE  to be x group in t.test
            # so "greater" is testing if TAD > non-TAD in terms of CTCF sites overlaps
            # i.e. is the mean number of CTCF sites in/near TAD boundaries > than
            # all bins within that TAD but not near either boundary
            factor(levels=c('Feature', 'Not Feature'))
    ) %>% 
    select(
        -c(
            nearest.boundary,
            ends_with(c('.match', '.ds', '.us'))
        )
    )
}

compute_CTCF_features_fisher_tests <- function(
    CTCF.df,
    resolution,
    HiCFeatureRadius.bins,
    n.CTCF.min.thresh,
    ...){
    # Count overlap of bins with enough CTCF sites & bins close enough to HiC Features
    contingency.table <- 
        CTCF.df %>% 
        group_by(isHiCFeature) %>% 
        summarize(
            n.bins.wo.CTCF=sum(n.overlaps <  n.CTCF.min.thresh),
            n.bins.w.CTCF=sum(n.overlaps  >= n.CTCF.min.thresh)
        ) %>%
        select(n.bins.w.CTCF, n.bins.wo.CTCF) %>%
        as.matrix()
    # Calculate enrichment stat for test
    bins.overlap      <- contingency.table[1,1]
    bins.w.feature    <- bins.overlap + contingency.table[1,2]
    bins.w.CTCF       <- bins.overlap + contingency.table[2,1]
    bins.total        <- sum(contingency.table)
    # Calculate enrichment stat for test
    expected          <- bins.w.CTCF * (bins.w.feature / bins.total)
    variance_term_1   <- (bins.total - bins.w.CTCF) / (bins.total - 1)
    variance_term_2   <- (bins.total - bins.w.feature) / (bins.total - 1)
    std_dev           <- sqrt(expected * variance_term_1 * variance_term_2)
    enrichment.zscore <- (bins.overlap - expected) / std_dev
    # Calculate fisher pvalue 
    fisher.test.row <- 
        contingency.table %>% 
        fisher.test(alternative='greater') %>% 
        tidy()
    # tidy data into single row tibble
    list(
        contingency.A=contingency.table[1,1],
        contingency.B=contingency.table[1,2],
        contingency.C=contingency.table[2,1],
        contingency.D=contingency.table[2,2],
        enrichment.zscore=enrichment.zscore,
        metric='overlaps',
        stat='n',
        test='fisher'
    ) %>% 
    as_tibble() %>% 
    bind_cols(fisher.test.row)
}

compute_CTCF_features_t_tests <- function(CTCF.df) {
    CTCF.df %>% 
    # pivot metrics so can do 1 test per metric+stat (i.e. per row)
    pivot_CTCF_enrichment_metrics() %>% 
    group_by(metric, stat) %>% 
    summarize(
        # welch's t-test of the mean CTCF metric is > near  TAD boundaries or not
        # results.t.test=
            t.test(
                value ~ isHiCFeature,
                alternative='greater' # only care if TADs are enriched for CTCFs, not depleted
            ) %>%
            tidy()
    ) %>%
    ungroup() %>% 
    # unnest(results.t.test) %>% 
    add_column(test='t.test')
}

compute_CTCF_features_corr_tests <- function(CTCF.df) {
    CTCF.df %>% 
    pivot_CTCF_enrichment_metrics() %>% 
    group_by(metric, stat) %>% 
    # calcualte test results for each CTCF metric + stat combo
    summarize(
        # are TADs more enriched for CTCF signal closer or farther from boundaries
        results.pearson=
            cor.test(
                x=dist.to.nearest,
                y=value,
                method='pearson',
                exact=FALSE,
                alternative='greater' # want closer to TAD ~ more CTCF signal, so +ve signal only
            ) %>% list(),
        results.spearman=
            cor.test(
                x=dist.to.nearest,
                y=value,
                method='spearman',
                exact=FALSE,
                alternative='greater' # want closer to TAD ~ more CTCF signal, so +ve signal only
        ) %>% list(),
    ) %>%
    ungroup() %>% 
    pivot_longer(
        starts_with('results.'),
        names_prefix='results.',
        names_to='test',
        values_to='test.results'
    ) %>% 
    rowwise() %>% 
    mutate(test.results=tidy(test.results)) %>% 
    unnest(test.results)
}

calculate_CTCF_enrichment_tests <- function(
    bins.df,
    features.df,
    test.type,
    ...){
    # paste0('row.idx=1; ', paste(colnames(tmp), '=tmp$', colnames(tmp), '[[row.idx]]', sep='', collapse='; '))
    # row.idx=1; annotation=tmp$annotation[[row.idx]]; feature.type=tmp$feature.type[[row.idx]]; motif=tmp$motif[[row.idx]]; method=tmp$method[[row.idx]]; Sample.Group=tmp$Sample.Group[[row.idx]]; resolution=tmp$resolution[[row.idx]]; bins.df=tmp$bins.df[[row.idx]]; features.df=tmp$features.df[[row.idx]]; HiCFeatureRadius.bins=tmp$HiCFeatureRadius.bins[[row.idx]]; n.CTCF.min.thresh=tmp$n.CTCF.min.thresh[[row.idx]]
    # Re-classify bins as being "HiC Features" if they are close enoughy i.e. 
    # within HiCFeatureRadius.bins bins of the actual feature
    assign_nearest_feature_to_bins(
        bins.df,
        features.df
    ) %>%
    {
        # calculate fisher's exact test pvalue of whether bins that are at/near TAD boundaries
        # are more likely to have > n.CTCF.min.thresh CTCF sites overlapping them than bins
        # that are at least HiCFeatureRadius.bins bins away from a TAD boundary
        if (test.type == 'fisher') {
            compute_CTCF_features_fisher_tests(
                CTCF.df=.,
                ...
            ) %>%
            select(-c(method, alternative))
        } else if (test.type == 't.test') {
            # Directly test if summary stats over CTCF site qvalues/scores are different closer to Features
            compute_CTCF_features_t_tests(CTCF.df=.) %>% 
            select(-c(estimate1, estimate2, parameter, statistic, method, alternative))
        # Test if distance to feature is correlated with CTCF stats
        } else if (test.type == 'corr.test') {
            compute_CTCF_features_corr_tests(CTCF.df=.)
        } else {
            stop(glue('invalid test type: {test.type}'))
        }
    } %>% 
    dplyr::rename_with(.cols=-c('test'), .fn=~str_replace(., '^', 'test.'))
}

calculate_all_feature_CTCF_enrichments <- function(
    overlaps.df,
    p.corr.group.cols=c(),
    ...){
    # overlaps.df=TAD.CTCF.overlaps.df;
    # overlaps.df=TAD.CTCF.overlaps.df %>% cross_join(expand_grid(HiCFeatureRadius.bins=c(0, 1, 2, 3), n.CTCF.min.thresh=c(1, 5, 10, 20, 30, 40) ));
    overlaps.df %>% 
            # {.} -> tmp; tmp
    # calculate enrichment test results across all conditions + params + hyper-params
    mutate(
        test.results=
            pmap(
            # future_pmap(
                .l=.,
                .f=calculate_CTCF_enrichment_tests,
                ...,
                .progress=TRUE
            )
    )  %>% 
    unnest(test.results) %>% 
    select(-c(bins.df, features.df)) %>% 
    # correct tests across resolutions + TAD calling methods + tests
    group_by(
        across(
            all_of(
                intersect(
                    colnames(.), 
                    c(
                        'test',
                        'test.metric',
                        'test.stat',
                        p.corr.group.cols
                    )
                )
            )
        )
    ) %>% 
    mutate(p.adj=p.adjust(test.p.value, method='BH')) %>% 
    ungroup() %>%
    # calculate log pvalues for plotting
    mutate(
        log.p.value=-log10(test.p.value),
        log.p.adj=-log10(p.adj)
    )
}

post_process_CTCF_test_results <- function(sig.expr=NULL){
    test.df %>% 
    {
        if (!is.null(sig.expr)) {
            mutate(., is.Enriched=ifelse(!! rlang::parse_expr(sig.expr), sig.expr, 'N.S.'))
        } else {
            .
        }
    } %>% 
    {
        if ('motif' %in% colnames(.)) {
            mutate(., motif=str_remove(motif, '_CORE_vertebrates'))
        } else {
            .
        }
    } %>% 
    {
        if ('HiCFeatureRadius.bins' %in% colnames(.)) {
            mutate(
                .,
                HiCFeatureRadius.bins.fct=
                    glue('<= {HiCFeatureRadius.bins} away') %>% 
                    fct_reorder(HiCFeatureRadius.bins)
            )
        } else {
            .
        }
    } %>% 
    {
        if ('n.CTCF.min.thresh' %in% colnames(.)) {
            mutate(
                .,
                n.CTCF.min.thresh.fct=
                    glue('>= {n.CTCF.min.thresh} CTCFs') %>% 
                    fct_reorder(n.CTCF.min.thresh)
            )
        } else {
            .
        }
    } %>% 
    rename(
        'metric'=test.metric,
        'stat'=test.stat
    )
}

###################################################
# Misc
###################################################
nest_data_for_plotting <- function(
    results.df,
    plot.file.suffix,
    nest_cols=NULL){
    # results.df=TAD.boundary.CTCFs.df; nest_cols=c('method', 'TAD.params', 'resolution'); plot.file.suffix='-overlap.enrichment.pdf'
    nest_cols <- 
        c(
            'feature.type',
            'annotation',
            nest_cols
        )
    # print(nest_cols)
    dir_glue_str <- 
        paste0(
            nest_cols,
            '_',
            '{', nest_cols, '}',
            collapse='/'
        )
    # print(dir_glue_str)
    results.df %>% 
    nest(plot.df=-all_of(nest_cols)) %>% 
    mutate(
        plot_dir=
            file.path(
                PLOT_DIR,
                glue(dir_glue_str)
            ),
        results_file=
            file.path(
                plot_dir,
                glue('{annotation}-{feature.type}-{plot.file.suffix}')
            )
    ) %>%
    select(plot.df, results_file)
}

