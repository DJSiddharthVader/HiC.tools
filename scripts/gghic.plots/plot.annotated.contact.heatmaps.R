######################################################################
# Dependencies
######################################################################
library(here)
# here::i_am('scritps/gghic.plots/plot.annotated.contact.heatmaps.R')
here::i_am('scripts/gghic.plots/plot.annotated.contact.heatmaps.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(InteractionSet)
    library(gghic)
    library(patchwork)
    source(file.path(BASE_DIR, 'scripts', 'locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(BASE_DIR, 'scripts', 'constants.R'))
    source(file.path(SCRIPT_DIR, 'utils.plot.R'))
    source(file.path(SCRIPT_DIR, 'gghic.plots', 'utils.gghic.R'))
    source(file.path(SCRIPT_DIR, 'TADs',  'utils.TADs.R'))
    source(file.path(SCRIPT_DIR, 'loops', 'utils.loops.R'))
    library(magrittr)
    library(tidyverse)
})

######################################################################
# Load Annotation data
######################################################################
tads.df   <- load_TADs_for_gghic()
loops.df  <- load_loops_for_gghic()
tracks.df <- NULL
# tracks.df <- 
#     check_cached_results(
#         results_file=GGHIC_TRACK_DATA,
#         # force_redo=TRUE,
#         results_fnc=load_tracks_for_gghic,
#         resolutions=resolutions
#     )

######################################################################
# Set up plotting params
######################################################################
# List of regions to plot annotated heatmaps for
regions.df <- 
    GENOMIC_REGIONS %>%
    filter(region.dist < 6e7, region.chr == 'chr16') %>% 
    head(4) %>% 
    select(region, region.chr, region.UCSC) %>%
    dplyr::rename(
        'chr'=region.chr,
        'focus'=region.UCSC,
        'region.title'=region
    )
# Heatmap params
hyper.params.df <- 
    expand_grid(
        resolution=c(50, 25, 10) * 1e3,
        weight=c('balanced', 'raw')
    )
# List all samples to plot
samples.df <- 
    list_mcool_files(only_use_included_samples=FALSE) %>%
    filter(isMerged) %>% 
    dplyr::select(-c(CloneID, TechRepID,  ReadFilter, isMerged)) %>% 
    dplyr::rename(cooler_path=filepath) %>% 
    mutate(SampleID=str_remove(SampleID, '.Merged.Merged'))

######################################################################
# Make all plots
######################################################################
# Match annotation data <-> SampleID <-> regions <-> plot params 
set_up_gghic_plot_param_sets(
        samples.df=samples.df,
        regions.df=regions.df,
        hyper.params.df=hyper.params.df,
        tads.df=tads.df,
        loops.df=loops.df,
        tracks.df=tracks.df
    ) %>% 
    # gghic specific params
    add_column(
        length_ratio=1,
        # ideogram_width_ratio=1/30,
        # ideogram_fontsize=10,
        # ideogram_colour="red",
        # ideogram_fill="#FFE3E680",
        gtf_path=GENOME_GTF_FILE,
        annotation_style="arrow",
        maxgap=-1,
        # annotation_width_ratio=1/50,
        # annotation_spacing_ratio=1/3,
        annotation_fontsize=10,
        annotation_colour="#48CFCB",
        annotation_fill="#48CFCB",
        include_ncrna=FALSE,
        # track_width_ratio=1/20,
        # track_spacing_ratio=1/2,
        # track_fill="black",
        # track_fontsize=5,
        tad_colour="#00ff83",
        loop_style="arc",
        stroke=1,
        # loop_colour="black",
        # loop_fill=NA,
        # expand_xaxis=FALSE,
        # expand_left=NULL,
        # expand_right=NULL,
    ) %>% 
    # generated pdfs for each plot
    pmap(
        .l=., 
        .f=plot_all_regions_gghic,
        .progress=TRUE
    )

