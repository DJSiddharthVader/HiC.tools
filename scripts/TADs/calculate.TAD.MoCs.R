################################################################################
# Dependencies
################################################################################
library(here)
BASE_DIR <- here()
suppressPackageStartupMessages({
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.loading.results.R'))
    # source(file.path(SCRIPT_DIR, 'TADs/utils.Comparing.TADs.R'))
    library(tidyverse)
    library(magrittr)
})

################################################################################
# Load TADs
################################################################################
TAD.sets.df <- 
    load_per_condition_TADs() %>% 
    convert_SampleID_to_SampleGroup() %>% 
    select(
        resolution,
        TAD.method, TAD.params, TAD.metric,
        Sample.Group,
        FeatureID, chr, start, end
        # TAD.length, TAD.start.score, TAD.end.score
    ) %>% 
    mutate(seqnames=chr) %>% 
    nest(
        regions.df=
            c(
              FeatureID, seqnames, start, end
              # TAD.length, TAD.start.score, TAD.end.score
            )
    )

################################################################################
# Calculate MoC for all pairs of TADs
################################################################################
# set up TAD Comparisons
ALL_SAMPLE_GROUP_COMPARISONS %>% 
    left_join(
        TAD.sets.df,
        relationship='many-to-many',
        by=join_by(Sample.Group.Numerator == Sample.Group)
    ) %>% 
    left_join(
        TAD.sets.df,
        suffix=c('.P1', '.P2'),
        relationship='many-to-many',
        by=
            join_by(
                resolution,
                TAD.method, TAD.params, TAD.metric,
                chr,
                Sample.Group.Denominator == Sample.Group
            )
    ) %>% 
    dplyr::rename(
        'Sample.Group.P1'=Sample.Group.Numerator,
        'Sample.Group.P2'=Sample.Group.Denominator
    ) %>%
    check_cached_results(
        results_file=ALL_MOC_FILE,
        force_redo=TRUE,
        results_fnc=calculate_all_MoCs,
        region.comparisons.df=.
    )

