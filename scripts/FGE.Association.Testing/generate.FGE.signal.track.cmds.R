################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/FGE.Association.Testing/generate.FGE.signal.track.cmds.R')
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'FGE.Association.Testing/utils.enrichment.R'))
    library(tidyverse)
    library(magrittr)
})

################################################################################
# CLI Args
################################################################################
parsed.args <- 
    handle_CLI_args(
        args=c('resolutions', 'threads', 'force'),
        has.positional=FALSE
    )
parsed.args$resolutions <- c(100, 50, 25) * 1e3
parsed.args$threads <- 1
if (parsed.args$threads > 1){
    # plan(multisession, workers=N_WORKERS_FOR_PARLLELIZATION)
    options(future.globals.maxSize=1.23 * 1024**3)
    plan(multisession, workers=parsed.args$threads)
} else {
    plan(sequential)
}

################################################################################
# Load raw FGE data
################################################################################
# List of genomic bins at various resolutions
genomic.bins.df <- 
    GENOME_BINS_FILES_DIR %>% 
    parse_results_filelist(suffix='.tsv') %>%
    select(-c(MatrixID)) %>% 
    filter(resolution %in% parsed.args$resolutions) %>% 
    mutate(
        bins.df=
            pmap(
                .l=list(filepath),
                .f=~ read_tsv(.x, show_col_types=FALSE) %>% dplyr::rename('chr'=chrom)
            )
        )
    # dplyr::rename('bins.filepath'=filepath)
# All FGE data for various signals
all.raw.FGE.data.df <- 
    # CTCF sites
    load_CTCF_sites(force.redo=FALSE) %>% 
    nest(sites.df=-c(motif)) %>%
    dplyr::rename('FGE.subtype'=motif) %>% 
    add_column(FGE.type='CTCF') %>% 
    # Other TF binding sites
    # bind_rows(
    #     load_TF_binding_sites() %>% 
    #     nest(sites.df=-c(TF.geneID)) %>% 
    #     dplyr::rename('FGE.subtype'=TF.geneID) %>% 
    #     dplyr::rename('FGE.type'='TFBS')
    # ) %>% 
    # ENCODE cCREs
    bind_rows(
        load_ENCODE_cCREs() %>% 
        nest(sites.df=-c(cCRE.type)) %>% 
        dplyr::rename('FGE.subtype'=cCRE.type) %>% 
        add_column('FGE.type'='cCRE')
    )

################################################################################
# MAP all FGEs to genomic bins they are inside, calculate summary stats for each bin 
################################################################################
# Now compute summary stats per bin using the overlap data just generated
compute_all_binwise_signals(
    all.raw.FGE.data.df=all.raw.FGE.data.df,
    genomic.bins.df=genomic.bins.df,
    silence=TRUE,
    force.redo=parsed.args$force.redo
)

