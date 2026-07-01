######################################################################
# Dependencies
######################################################################
library(ggplot2)
library(gghic)

######################################################################
# Load track data
######################################################################
load_TADs_for_gghic <- function(){
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>%
    filter(method == 'hiTAD') %>% 
    post_process_hiTAD_TAD_results() %>% 
    # add_column(weight='balanced') %>% 
    dplyr::select(-c(length)) %>% 
    mutate(chr2=chr) %>% 
    nest(TADs=c(chr, start, end)) %>% 
    mutate(SampleID=str_remove(SampleID, '.Merged.Merged')) %>% 
    dplyr::rename(
        'chr'=chr2,
        'TAD.method'=method
    )
}

load_loops_for_gghic <- function(){
    COOLTOOLS_LOOPS_RESULTS_FILE %>% 
    read_tsv(show_col_types=FALSE) %>%
    post_process_cooltools_dots_results() %>% 
    filter_loop_results() %>% 
    dplyr::rename(
        'start.P1'=anchor.left,
        'start.P2'=anchor.right
    ) %>% 
    mutate(
        end.P1=start.P1 + resolution,
        end.P2=start.P2 + resolution,
        chr.P1=chr,
        chr.P2=chr,
    ) %>% 
    dplyr::select(
        weight, resolution, SampleID, chr,
        chr.P1, start.P1, end.P1, 
        chr.P2, start.P2, end.P2,
    ) %>% 
    nest(
        loops=
            c(
                chr.P1, start.P1, end.P1, 
                chr.P2, start.P2, end.P2,
            )
    )
}

load_tracks_for_gghic <- function(){
    return(NULL)
}

######################################################################
# Set up plot obects + params
######################################################################
make_gghic_plot_obj <- function(
    cooler_path, resolution, 
    focus,
    TADs, loops, tracks=NULL,
    ...){
    # paste(colnames(tmp), '=tmp$', colnames(tmp), '[[1]]', sep='', collapse=';')
    # cooler_path=tmp$cooler_path[[1]];Edit=tmp$Edit[[1]];Celltype=tmp$Celltype[[1]];Genotype=tmp$Genotype[[1]];SampleID=tmp$SampleID[[1]];resolution=tmp$resolution[[1]];focus=tmp$focus[[1]];TAD.method=tmp$TAD.method[[1]];TADs=tmp$TADs[[1]];weight=tmp$weight[[1]];loops=tmp$loops[[1]]
    # make object with contacts
    cc <- 
        ChromatinContacts(
            cooler_path=cooler_path,
            focus=focus,
            resolution=as.integer(resolution)
        ) %>%
        import()
    # add TAD annotations
    if (!is.null(TADs)){
        features(cc, 'TADs') <- GRanges(TADs)
    }
    # add loop annotations
    if (!is.null(loops)){
        features(cc, 'loops') <- 
            GInteractions(  
                anchor1=
                    loops %>% 
                    dplyr::select(ends_with('.P1')) %>% 
                    rename_with(~ str_remove(.x, '.P1')) %>% 
                    GRanges(),
                anchor2=
                    loops %>% 
                    dplyr::select(ends_with('.P2')) %>% 
                    rename_with(~ str_remove(.x, '.P2')) %>% 
                    GRanges()
            ) %>% 
            pairs() %>% 
            makeGInteractionsFromGRangesPairs()
    }
    # add 1D track data
    if (!is.null(tracks)){
        features(cc, 'tracks') <- GRanges(tracks)
    } 
    # return annotated contacts plot object
    cc
}

enumerate_all_param_combinations <- function(
    samples.df,
    regions.df,
    hyper.params.df,
    tads.df=NULL,
    loops.df=NULL,
    tracks.df=NULL,
    ...){
    # Specify regions to plot + all hyper params
    regions.df %>% 
    cross_join(hyper.params.df) %>% 
    # Specify which samples data to plot
    cross_join(samples.df)  %>% 
    # Join TADs
    {
        if (!is.null(tads.df)) {
            left_join(
                .,
                tads.df,
                by=
                    join_by(
                        weight,
                        resolution,
                        SampleID,
                        chr
                    )
            )
        } else {
            add_column(., tads.df=NULL)
        }
    } %>% 
    # Join Loops
    {
        if (!is.null(loops.df)) {
            left_join(
                .,
                loops.df,
                by=
                    join_by(
                        weight,
                        resolution,
                        SampleID,
                        chr,
                    )
            )
        } else {
            add_column(., loops.df=NULL)
        }
    } %>% 
    # Join tracks
    {
        if (!is.null(tracks.df)) {
            left_join(
                .,
                tracks.df,
                by=
                    join_by(
                        resolution,
                        SampleID,
                        chr
                    )
            )
        } else {
            add_column(., tracks=list(NULL))
        }
    } 
}

set_up_gghic_plot_param_sets <- function(...){
    # Set up plot params
    params.df <- 
        enumerate_all_param_combinations(...) %>% 
        rename('resolution.int'=resolution) %>% 
        mutate(
            resolution=scale_numbers(resolution.int),
            panel.title=glue('{SampleID} @ {resolution}'),
            # region.title=glue('{SampleID} @ {resolution} in {region.title} ({focus})'),
            plot.title=glue('{region.title} ({focus})')
        )
        # arrange(weight, resolution, focus) %>% 
    # Group specific plots together by  condition
    # Plot all genotypes for each condition together in a single plot
    group.by.celltype.params.df <- 
        params.df %>% 
        add_column(group='Celltype') %>% 
        nest(
            plot.df=
                c(
                    # Edit, Celltype, Genotype, SampleID,
                    Genotype, 
                    SampleID,
                    cooler_path,
                    resolution.int,
                    TADs, loops, tracks,
                    panel.title
                )
        )
    # Plot all genotypes for each condition together in a single plot
    group.by.genotype.params.df <- 
        params.df %>% 
        add_column(group='Genotype') %>% 
        nest(
            plot.df=
                c(
                    Celltype, 
                    SampleID,
                    cooler_path,
                    resolution.int,
                    TADs, loops, tracks,
                    panel.title
                )
        )
    # Set output file path t obe descriptive
    bind_rows(
        group.by.celltype.params.df,
        group.by.genotype.params.df
    ) %>% 
    mutate(
        results_file=
            file.path(
                GGHIC_DIR,
                glue('weight_{weight}'),
                glue('resolution_{scale_numbers(resolution)}'),
                glue('grouping_{group}'),
                glue('{region.title}-gghic.plot.pdf')
            )
    )
}

######################################################################
# Make plots
######################################################################
plot_gghic <- function(
    cooler_path,
    resolution.int,
    focus,
    weight,
    TADs,
    loops,
    tracks,
    panel.title,
    gtf_path=NULL,
    ...){
    # paste(colnames(plots.df), '=plots.df$', colnames(plots.df), '[[1]]', sep='', collapse=';')
    # cooler_path=plots.df$cooler_path[[1]];Edit=plots.df$Edit[[1]];Celltype=plots.df$Celltype[[1]];Genotype=plots.df$Genotype[[1]];SampleID=plots.df$SampleID[[1]];chr=plots.df$chr[[1]];resolution.int=plots.df$resolution.int[[1]];weight=plots.df$weight[[1]];TAD.method=plots.df$TAD.method[[1]];TADs=plots.df$TADs[[1]];loops=plots.df$loops[[1]];tracks=plots.df$tracks[[1]];region.title=plots.df$region.title[[1]];focus=plots.df$focus[[1]];resolution=plots.df$resolution[[1]]
    # make plot obj with contacts + annotation data
    make_gghic_plot_obj(
        cooler_path=cooler_path,
        resolution=resolution.int, 
        focus=focus,
        TADs=TADs, 
        loops=loops,
        tracks=tracks
    ) %>% 
    # Make gghic plot with specified annotations
    {
        gghic(
            .,
            scale_column=weight,
            genome="hg38", 
            gtf_path=gtf_path,
            annotation=!is.null(gtf_path),
            ideogram=TRUE, 
            loop=TRUE,
            tad=TRUE,
            tad_is_0based=TRUE,
            track=TRUE,
            ...
        ) +
        ggtitle(panel.title) +
        theme_hic()
    } #%>% 
    # Handle faceting + scaling + theme options
    # post_process_plot(...)
}

plot_all_regions_gghic <- function(
    plot.df,
    focus,
    weight,
    plot.title,
    results_file=NA,
    width=8,
    height=14,
    ...){
    plot.df %>%
    mutate(
        plot_obj=
            pmap(
                .l=.,
                .f=plot_gghic,
                focus=focus,
                weight=weight,
                ...,
                .progress=TRUE
            )
    ) %>%
    pull(plot_obj) %>% 
    # pmap(
    #     .f=post_process_plot,
    #     ...,
    #     .progress=FALSE
    # ) %>% 
    {
        wrap_plots(
            .,
            ncol=1
        ) +
        plot_annotation(title=plot.title)
    } %>%
    {
        if (!is.na(results_file) & !file.exists(results_file)){
            ggsave(
                results_file, 
                plot=.,
                width=width,
                height=height,
                units='in'
            )
        } else {
            .
        }
    }
}

