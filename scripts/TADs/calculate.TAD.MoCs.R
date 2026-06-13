################################################################################
# Dependencies
################################################################################
library(here)
here::i_am('scripts/TADs/calculate.TAD.MoCs.R')
BASE_DIR <- here()
source(file.path(BASE_DIR,   'scripts/constants.R'))
source(file.path(BASE_DIR,   'scripts/locations.R'))
source(file.path(SCRIPT_DIR, 'utils.data.R'))
source(file.path(SCRIPT_DIR, 'TADs/utils.Comparing.TADs.R'))
library(tidyverse)
library(magrittr)

################################################################################
# Load TADs
################################################################################
nested.TADs.df <- 
    ALL_TAD_RESULTS_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    # select(
    #     resolution,
    #     TAD.method, TAD.params, TAD.metric,
    #     Sample.Group,
    #     chr, start, end,
    #     TAD.length, TAD.start.score, TAD.end.score, starts_with('TAD.inner.')
    # ) %>% 
    convert_SampleID_to_SampleGroup() %>% 
    select(-c(SampleID)) %>% 
    nest(
        TADs=
            -c(
                Sample.Group,
                # SampleID,
                isMerged,
                resolution,
                TAD.method, TAD.params, TAD.metric,
                chr
            )
    )

################################################################################
# Calculate MoC for all pairs of TADs
################################################################################
check_cached_results(
    results_file=ALL_MOC_FILE,
    force_redo=TRUE,
    results_fnc=calculate_all_MoCs,
    nested.TADs.df=nested.TADs.df,
    sample.group.comparisons=
        ALL_SAMPLE_GROUP_COMPARISONS %>% 
        dplyr::rename(
            'Sample.Group.P1'=Sample.Group.Numerator,
            'Sample.Group.P2'=Sample.Group.Denominator
        ),
    sampleID_col='Sample.Group',
    suffixes=
        c(
            'Numerator',
            'Denominator'
        ),
    pair_grouping_cols=
        c(
            'isMerged',
            'resolution',
            'TAD.method', 'TAD.params', 'TAD.metric',
            'chr'
        )
)

