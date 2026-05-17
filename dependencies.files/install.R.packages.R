install.packages(
    c(
        'ggplot2',
        'ggpubr',
        'ggh4x',
        'ggridges',
        'GGally',
        'scales',
        'ggpointdensity',
        'viridis',
        # 'ggstatsplot',
        'here',
        'stringi',
        'glue',
        'optparse',
        'tictoc',
        'devtools',
        'tidyverse',
        'dplyr',
        'magrittr',
        'cluster',
        'furrr',
        'future',
        'PRIMME',
        'Matrix',
        # 'HiCcompare',
        # 'stringr',
        # 'RMTstat',
        # 'strawr'
        'BiocManager'
    )
)
# install.packages('diffdomain_0.1.0.tar.gz', repos=NULL,type='source')
BiocManager::install(
    c(
        'TADCompare',
        'multiHiCcompare'
    )
)
pak::pak("jasonwong-lab/gghic")
devtools::install_github("paulsengroup/hictkR")
    # pak::pak("paulsengroup/hictkR")
