################################################################################
# Dependencies
################################################################################
library(here)
BASE_DIR <- here()
suppressPackageStartupMessages({
    # library(hictkR)
    source(file.path(BASE_DIR,   'scripts/constants.R'))
    source(file.path(BASE_DIR,   'scripts/locations.R'))
    source(file.path(SCRIPT_DIR, 'utils.data.R'))
    source(file.path(SCRIPT_DIR, 'utils.loading.results.R'))
    source(file.path(SCRIPT_DIR, 'compartments/utils.compartments.R'))
    library(tidyverse)
    library(magrittr)
    library(furrr)
})
parsed.args <- 
    handle_CLI_args(
        args=c('threads', 'force', 'resolutions'),
        has.positional=FALSE
    )
plan(multisession, workers=length(availableWorkers()))

################################################################################
# Compute A/B Compartment MoCs between conditions
################################################################################
# Load Compartment regions
compartment.sets.df <- 
    load_per_condition_compartment_regions() %>% 
    select(
        resolution,
        compartment.params,
        Sample.Group,
        compartment,
        FeatureID, chr, start, end
    ) %>% 
    mutate(seqnames=chr) %>% 
    nest(
        regions.df=
            c(
              FeatureID, seqnames, start, end
            )
    )
# set up Comparisons for pre-specified pairs of conditions
ALL_SAMPLE_GROUP_COMPARISONS %>% 
    left_join(
        compartment.sets.df,
        relationship='many-to-many',
        by=join_by(Sample.Group.Numerator == Sample.Group)
    ) %>% 
    left_join(
        compartment.sets.df,
        suffix=c('.P1', '.P2'),
        relationship='many-to-many',
        by=
            join_by(
                resolution,
                compartment.params,
                compartment, # A or B, dont want to comapre A vs B overlap, just A vs A and B vs B
                chr,
                Sample.Group.Denominator == Sample.Group
            )
    ) %>% 
    dplyr::rename(
        'Sample.Group.P1'=Sample.Group.Numerator,
        'Sample.Group.P2'=Sample.Group.Denominator
    ) %>%
    # compute MoC for all pairs of sets of compartments 
    # for each chr, for each comparison, for A/B compartments separately
    check_cached_results(
        results_file=ALL_COMPARTMENT_MOCS_FILE,
        # force_redo=parsed.args$force.redo,
        force_redo=TRUE,
        results_fnc=calculate_all_MoCs,
        region.comparisons.df=.
    )

################################################################################
# Test Differece in PC1 Scores for WT Segments
################################################################################
# nest all bins into rows
nested.bins.df <- 
    ALL_COMPARTMENT_BINWISE_FILE %>%
    read_tsv(show_col_types=FALSE) %>% 
    select(
        resolution,
        compartment.params,
        Sample.Group,
        compartment,
        chr, start, end,
        score
    ) %>% 
    mutate(seqnames=chr) %>% 
    nest(
        segments=
            c(
                compartment,
                seqnames, start, end,
                score
            )
    ) %>%
    mutate(segments=pmap(.l=list(segments), .f=as_granges))
# set up Comparisons for pre-specified pairs of conditions
ALL_SAMPLE_GROUP_COMPARISONS %>% 
    left_join(
        nested.bins.df,
        relationship='many-to-many',
        by=join_by(Sample.Group.Numerator == Sample.Group)
    ) %>% 
    left_join(
        nested.bins.df,
        suffix=c('.Numerator', '.Denominator'),
        relationship='many-to-many',
        by=
            join_by(
                resolution,
                compartment.params,
                chr,
                Sample.Group.Denominator == Sample.Group
            )
    ) %>% 
    # compute differences in binwise PC1 scores between conditions across segments
    # for each chr, for each comparison, for A/B compartments separately
    check_cached_results(
        results_file=ALL_COMPARTMENT_SEGMENT_TEST_RESULTS_FILE,
        # force_redo=parsed.args$force.redo,
        force_redo=TRUE,
        results_fnc=calculate_all_binwise_difference_stats,
        segments.df=.
    )

################################################################################
# Compare distance and PC1 difference of compartment swtiches across Conditions
################################################################################
# load all compartment switches
nested.switches.df <- 
    load_per_condition_compartment_switches() %>% 
    select(
        resolution,
        compartment.params,
        Sample.Group,
        # FeatureID, chr, start, end,
        chr, start, end,
        compartment.change,
        score, 
        score.change, strength.lvl.change, 
    ) %>% 
    # mutate(seqnames=chr) %>% 
    dplyr::rename('seqnames'=chr) %>% 
    nest(
        switches=
            c(
                # FeatureID, seqnames, start, end,
                seqnames, start, end,
                # compartment.change,
                score, 
                score.change, strength.lvl.change, 
            )
    ) %>%
    mutate(switches=pmap(.l=list(switches), .f=as_granges))
# get all nearest Numerator switches to each Denominator switch
ALL_SAMPLE_GROUP_COMPARISONS %>% 
    left_join(
        nested.switches.df,
        relationship='many-to-many',
        by=join_by(Sample.Group.Numerator == Sample.Group)
    ) %>% 
    left_join(
        nested.switches.df,
        suffix=c('.Numerator', '.Denominator'),
        relationship='many-to-many',
        by=
            join_by(
                resolution,
                compartment.params,
                # chr,
                Sample.Group.Denominator == Sample.Group
            )
    ) %>% 
    # for each Denominator (WT) swtich, map the nearest Numeartor (DEL) switch
    # then for all switch pairs, compare distance and FC, 
    # FC for a switch is the difference in PC1 score between the switch bin and the 
    # next adjacent bin (which is annotated as a different compartment)
    check_cached_results(
        results_file=ALL_COMPARTMENT_SWITCH_DIFFERENCE_RESULTS_FILE,
        # force_redo=parsed.args$force.redo,
        force_redo=TRUE,
        results_fnc=match_all_nearest_switches_between_conditions,
        switches.df=.
    )

