###################################################
# Depdendencies
###################################################

###################################################
# Load cooltools results
###################################################
quantize_track <- function(
    scores,
    nbins=4){
    cut(
        x=scores,
        breaks=
            quantile(
                E1,
                probs=
                    c(
                        seq(0.0, 0.5, length.out=nbins+1),
                        seq(0.5, 1.0, length.out=nbins+1),
                    ) %>%
                    unique()
                na.rm=TRUE
            ),
        labels=
            c(
              paste
            )
    )
}

load_cooltools_compartment_results <- function(
    filepath,
    ...){
    # filepath=tmp$filepath[[1]]
    filepath %>%
    read_tsv(
        show_col_types=FALSE,
        progress=FALSE
    ) %>%
        {.} -> e1.df; e1.df
    e1.df %>% 
    # classify bins by first eigenvector
    mutate(
        compartment.n2=
            case_when(
                E1 >  0   ~ 'A',
                E1 <  0   ~ 'B',
                E1 == 0   ~ '0',
                is.na(E1) ~ NA,
                TRUE      ~ as.character(E1)
            )
    ) %>%
    group_by(chrom, compartment.n2) %>% 
    mutate(
        across(
            .fn=,
            .names=''
        )
    )
        compartment.n6=
            cut(
                x=E1,
                breaks=
                    quantile(
                        E1,

                        probs=c(0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0),
                        na.rm=TRUE
                    ),
                labels=
                    c(
                        'Strong A',
                        'Med A',
                        'Weak A',
                        # 'Near 0',
                        'Weak B',
                        'Med B',
                        'Strong B'
                    )
            )
        compartment.n20=
            cut(
                x=E1,
                breaks=
                    quantile(
                        E1,
                        probs=seq(0, 1, 0.05),
                        na.rm=TRUE
                    ),
                labels=
                    c(
                    )
            )
    ) %>% 
    arrange(chrom, start, end) %>% 
    group_by(chrom) %>% 
    mutate(
        across(
            .cols=starts_with('compartment.'),
            .fn=\(x, ...) glue('{x}->{lead(x, n=1L)}'), # need the ... for dplyr reasons idk
            .names='switch.{str_remove(.col, "^compartment.")}'
        )
    ) %>% 
    ungroup() %>% 
        count(compartment.n6, switch.n6) %>% print(n=Inf)
    select(
        chrom, start, end,
        E1,
        starts_with('compartment.', 'switch.')
    ) %>% 
    filter(!is.na(E1)) %>% 
    dplyr::rename('chr'=chrom)
}

load_all_cooltools_compartment_results <- function(resolutions=NULL){
    COMPARTMENTS_RESULTS_DIR %>%
    parse_results_filelist(
        filename.column.name='SampleID',
        suffix='-.cis.vecs.tsv'
    ) %>%  
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
                 .progress=TRUE
            )
    ) %>%
    unnest(compartments) %>% 
    select(-c(filepath)) 
}

post_process_cooltools_compartment_results <- function(results.df){
    results.df %>%
    mutate(
        compartment.binary=
            factor(
                compartment.binary,
                levels=
                    c(
                        'A',
                        'B'
                    )
            ),
        compartment.quantile=
            factor(
                compartment.quantile,
                levels=
                    c(
                        'Strong A',
                        '  Weak A',
                        '  Near 0',
                        '  Weak B',
                        'Strong B'
                    )
            ),
        compartment.switch=
            factor(
                compartment.switch,
                levels=
                    c(
                        'A->A'
                        'B->B'
                        'A->B'
                        'B->A'
                    )
            )
    )
}

