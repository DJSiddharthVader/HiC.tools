###################################################
# Dependencies
###################################################
library(ggplot2)
library(ggpubr)
library(ggh4x)
library(ggridges)
library(GGally)
library(scales)
library(ggpointdensity)
library(viridis)
library(furrr)

###################################################
# Transform data for plotting
###################################################
calc_pct <- function(
    count.df,
    cols_exclude=c(),
    col_pct=NULL){
    # count.df=n.TADCompare.df; cols_exclude=c('chr', 'DifferenceType', 'Enriched.Condition'); col_pct=NULL
    # calculate relative frequency from count data
    count.df %>% 
    {
        if ('n' %in% colnames(.)){
            group_by(., across(-c(cols_exclude, 'n'))) %>% 
            summarize(n=sum(n))
        } else {
            group_by(., across(-c(cols_exclude))) %>% 
            count()
        }
    } %>% 
    {
        if (!is.null(col_pct)) {
            ungroup(.) %>%
            group_by(across(-c(col_pct, 'n')))
        } else {
            .
        }
    } %>% 
    mutate(total=sum(n)) %>% 
    ungroup() %>% 
    mutate(pct=n / total) %>% 
    mutate(n.label=glue('n = {n}')) %>% 
    mutate(pct.label=glue('{round(100 * n / total, digits=1)}%')) %>% 
    mutate(n.and.pct.label=glue('{pct.label}\n({n.label})'))
}

copy_data_along_inclusive_intervals <- function(
    plot.df,
    input_colname,
    output_colname,
    thresholds,
    comparison_op='<',
    decreasing=TRUE){
    # plot.df=loops.df; input_colname='log10.qval'; thresholds=c(1, 10, 50, 100, 200); output_colname='sig.band'; comparison_op='>'; decreasing=TRUE
    # turn comparison operator into a function i.e. '<' becomes function(x, y) {x < y}
    comparison.fn <- match.fun(comparison_op)
    # put all thresholds as a tibble column
    thresholds %>% 
    tibble(threshold=.) %>%
    # for every threshold get pretty name and make ordered facet 
    mutate(
        {{output_colname}} :=
            paste(
                input_colname,
                comparison_op,
                threshold
            ) %>%
            fct_reorder(thresholds, .desc=decreasing)
    ) %>% 
    # Now for every threshold, join all rows from the input dataset
    cross_join(plot.df) %>%
    # Only keep rows where the specified column meets the threshold
    filter(comparison.fn(!!sym(input_colname), threshold)) %>%
    select(-c(threshold))
}

###################################################
# Handling/formatting plots
###################################################
make_ggtheme <- function(...){
    theme(
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        axis.line=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks.length=unit(3, "pt"),
        axis.title=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=10
            ),
        axis.text=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=8
            ),
        ...
    )
}

build_axis_fnc <- function(
    figure,
    scale.mode,
    scale.axis,
    limits){
    # fetch axis data to check if data is discerete or not
    scale.data <- rlang::eval_tidy(rlang::quo_squash(figure@mapping[[scale.axis]]), figure@data)
    # Determine axis type
    axis.type <- 
        case_when(
            scale.mode == 'pct'      ~ 'continuous',
            scale.mode == 'm'        ~ 'continuous',
            scale.mode == 'mb'       ~ 'continuous',
            scale.mode == 'log10'    ~ 'log10',
            is.character(scale.data) ~ 'discrete',
            is.factor(scale.data)    ~ 'discrete',
            # scale.mode == 'binned'   ~ 'binned',
            scale.mode == ''         ~ 'continuous'
        )
    # build axis scaling functiona name for ggplot2 and eval to return function object
    glue('scale_{scale.axis}_{axis.type}') %>%
    as.symbol() %>%
    eval()

}

scale_axis <- function(
    figure,
    scale.axis='x',
    scale.mode='',
    log.base=10,
    axis.label.accuracy=0.1,
    n.breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    # figure out which axis function to call
    axis_fnc <- 
        build_axis_fnc(
            figure=figure,
            scale.axis=scale.axis,
            scale.mode=scale.mode,
            limits=limits
        )
    # Set axis scaling/labeling based on scale.mode argument
    # Scaling for percentages
    if (scale.mode == 'pct') {
        figure +
        axis_fnc(
            expand=expand,
            n.breaks=n.breaks,
            limits=limits,
            labels=label_percent(),
            ...
        )
    # Human scaling for large numbers e.g. 1K, 1M, 1B
    } else if (scale.mode == 'm') {
        figure +
        axis_fnc(
            expand=expand,
            n.breaks=n.breaks,
            limits=limits,
            labels=
                label_number(
                    scale_cut=cut_short_scale(),
                    accuracy=axis.label.accuracy
                ),
            ...
        )
    # Human scaling in bytes e.g. 1Kb, 1Mb, 1Gb. 
    # Useful for genomic coordinates + sizes 
    } else if (scale.mode == 'mb') {
        figure +
        axis_fnc(
            expand=expand,
            n.breaks=n.breaks,
            limits=limits,
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis.label.accuracy
                ),
            ...
        )
    # log scale + include minor ticks, uses base 10 by default
    } else if (scale.mode == 'log') {
        figure +
        axis_fnc(
            expand=expand,
            limits=limits,
            guide='axis_logticks',
            labels=
                label_log(
                    base=log.base,
                    digits=max(1, -log10(axis.label.accuracy)),
                    signed=FALSE
                ),
            ...
        )
    # Discrete set of labeld i.e. DEL,WT
    } else if (scale.mode == 'discrete') {
        figure +
        axis_fnc(expand=expand)
    # No scaling, just set limits and label rounding
    } else if (scale.mode == '') {
        if (is.null(limits)) {
            figure + 
            axis_fnc(
                expand=expand,
                labels=
                    function(x) {
                        format(x, digits=max(1, -log10(axis.label.accuracy)))
                    }
            )
        } else {
            figure + 
            axis_fnc(
                limits=limits,
                expand=expand,
                ...
            )
        }
    }
}

add_faceting <- function(
    figure,
    # space='fixed',
    # solo_line=TRUE,
    # independent='',
    # axes=FALSE,
    # trim_blank=TRUE,
    facet.group=NULL,
    facet.col=NULL,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.ncol=NULL,
    ...){
    # Facet as specified
    facet.formula <- 
    if (!is.null(facet.col) & !is.null(facet.row)) {
        figure <- 
            figure +
            facet_nested(
                paste(
                    paste(facet.row, collapse=' + '),
                    paste(facet.col, collapse=' + '),
                    sep=' ~ '
                ) %>%
                formula(),
                # space=space,
                # solo_line=solo_line,
                # independent=independent,
                # axes=axes,
                # trim_blank=trim_blank,
                ...
            )
    } else if (!is.null(facet.row)) {
        figure <- 
            figure +
            facet_nested(
                formula(glue('{paste(facet.row, collapse=" + ")} ~ .')),
                # space=space,
                # solo_line=solo_line,
                # independent=independent,
                # axes=axes,
                # trim_blank=trim_blank,
                ...
            )
    } else if (!is.null(facet.col)) {
        figure <- 
            figure +
            facet_nested(
                formula(glue('~ {paste(facet.col, collapse=" + ")}')),
                # space=space,
                # solo_line=solo_line,
                # independent=independent,
                # trim_blank=trim_blank,
                # axes=axes,
                ...
            )
    } else if (!is.null(facet.group)) {
        figure <- 
            figure +
            facet_wrap2(
                vars(!!sym(facet.group)),
                nrow=facet.nrow,
                ncol=facet.ncol
                # ...
            )
    }
    figure
}

post_process_plot <- function(
    figure,
    plot.elements=NULL,
    theme.obj=NULL,
    scales='fixed',
    space='fixed',
    independent=FALSE,
    drop=TRUE,
    axes='margins',
    margins=FALSE,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.col=NULL,
    facet.ncol=NULL,
    facet.group=NULL,
    x.scale.mode='',
    x.log.base=10,
    x.axis.label.accuracy=0.1,
    x.n.breaks=NULL,
    x.limits=NULL,
    x.expand=c(0.00, 0.00, 0.00, 0.00),
    y.scale.mode='',
    y.log.base=10,
    y.axis.label.accuracy=0.1,
    y.n.breaks=NULL,
    y.limits=NULL,
    y.expand=c(0.00, 0.00, 0.00, 0.00),
    color.scale.mode='',
    color.log.base=10,
    color.axis.label.accuracy=0.1,
    color.n.breaks=NULL,
    color.limits=NULL,
    color.expand=c(0.00, 0.00, 0.00, 0.00),
    fill.scale.mode='',
    fill.log.base=10,
    fill.axis.label.accuracy=0.1,
    fill.n.breaks=NULL,
    fill.limits=NULL,
    fill.expand=c(0.00, 0.00, 0.00, 0.00),
    # linetype.scale.mode='',
    # linetype.log.base=10,
    # linetype.axis.label.accuracy=0.1,
    # linetype.n.breaks=NULL,
    # linetype.limits=NULL,
    # linetype.expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    figure %>% 
    add_faceting(
        scales=scales,
        space=space,
        independent=independent,
        axes=axes,
        drop=drop,
        margins=margins,
        facet.row=facet.row,
        facet.col=facet.col,
        facet.group=facet.group,
        facet.nrow=facet.nrow,
        facet.ncol=facet.ncol
    ) %>% 
    # Set x-axis scaling (log, Mb, percent etc.)
    scale_axis(
        scale.axis='x',
        scale.mode=x.scale.mode,
        log.base=x.log.base,
        axis.label.accuracy=x.axis.label.accuracy,
        n.breaks=x.n.breaks,
        limits=x.limits,
        expand=x.expand
    ) %>%
    # Set y-axis scaling (log, Mb, percent etc.)
    scale_axis(
        scale.axis='y',
        scale.mode=y.scale.mode,
        log.base=y.log.base,
        axis.label.accuracy=y.axis.label.accuracy,
        n.breaks=y.n.breaks,
        limits=y.limits,
        expand=y.expand
    ) %>% 
    # Set color scaling (log, Mb, percent etc.)
    {
        if ('color' %in% names(figure@mapping)) {
            scale_axis(
                figure=.,
                scale.axis='color',
                scale.mode=color.scale.mode,
                log.base=color.log.base,
                axis.label.accuracy=color.axis.label.accuracy,
                n.breaks=color.n.breaks,
                limits=color.limits,
                expand=color.expand
            )
        } else {
            .
        }
    } %>% 
    # Set fill scaling (log, Mb, percent etc.)
    {
        if ('fill' %in% names(figure@mapping)) {
            scale_axis(
                figure=.,
                scale.axis='fill',
                scale.mode=fill.scale.mode,
                log.base=fill.log.base,
                axis.label.accuracy=fill.axis.label.accuracy,
                n.breaks=fill.n.breaks,
                limits=fill.limits,
                expand=fill.expand
            )
        } else {
            .
        }
    } %>% 
    # Set linetype scaling (log, Mb, percent etc.)
    # {
    #     if ('linetype' %in% names(figure@mapping)) {
    #         scale_axis(
    #             figure=.,
    #             scale.axis='linetype',
    #             scale.mode=linetype.scale.mode,
    #             log.base=linetype.log.base,
    #             axis.label.accuracy=linetype.axis.label.accuracy,
    #             n.breaks=linetype.n.breaks,
    #             limits=linetype.limits,
    #             expand=linetype.expand
    #         )
    #     } else {
    #         .
    #     }
    # } %>% 
    # Add theme elements, as either an object or individual args
    {
        if (!is.null(theme.obj)) {
            . + theme.obj
        } else {
            . + make_ggtheme() + theme(...)
        } 
    } %>% 
    # Add extra elements
    {
        if (!is.null(plot.elements)) {
            purrr::reduce(
                .x=plot.elements,
                .f=function(x, y) { x + y },
                .init=.
            )
        } else {
            .
        }
    }
}

###################################################
# Make tabs per plot in Rmd
###################################################
plot_figure_tabs <- function(
    plot.df,
    group.col,
    plot.fnc,
    header.lvl,
    nl.delim,
    figure.output.mode='rmd',
    grob.nrow=1,
    grob.ncol=NULL,
    ...){
    # List all the individual groups, generate 1 plot per group
    group.values <- 
        plot.df[[group.col]] %>% 
        as.factor() %>% 
        droplevels() %>% 
        levels()
    # Generate figures on mutually exclusive subsets of the data
    figures <- 
        # make plot with options for each group of the data
        future_pmap(
            .l=list(group.value=group.values),
            .f=
                function(group.value, group.col, plot.df, plot.fnc){
                    plot.df %>%
                    filter(!!sym(group.col) == group.value) %>%
                    plot.fnc(...)
                },
            group.col=group.col,
            plot.df=plot.df,
            plot.fnc=plot.fnc,
            .progress=FALSE
        )
    # print(group.values)
    # print(length(figures))
    # make a named list of figures
    names(figures) <- group.values
    # Print/combine/return plots as specified
    {
        # Print each figure under a md heading for Rmd notebooks
        if (figure.output.mode == 'rmd') {
            figures %>%
            names() %>% 
            sapply(
                FUN=
                    function(group.value, figures, header.lvl, nl.delim){
                        cat(
                            strrep('#', header.lvl), 
                            group.value,
                            nl.delim
                        )
                        print(figures[[group.value]])
                        cat(nl.delim)
                    },
                figures=figures,
                header.lvl=header.lvl,
                nl.delim=nl.delim
            )
        # merge all the plots into a single figure with labeled panels
        } else if (figure.output.mode == 'merged') {
            cat(
                strrep('#', header.lvl),
                nl.delim
            )
            cowplot::plot.grid(
                plotlist=figures,
                nrow=grob.nrow,
                ncol=grob.ncol,
                labels=group.values,
                axis='tb',
                align='hv'
            ) %>%
            print()
            cat(nl.delim)
        # just return all the figures
        } else if (figure.output.mode == 'return') {
            return(figures)
        } else {
            stop(glue('Invalid arg for figure.output.mode: {figure.output.mode}'))
        }
    }
}

make_tabs_recursive <- function(
    plot.df, 
    group.cols,
    current.header.lvl,
    plot.fnc,
    tabset.format,
    nl.delim,
    figure.output.mode,
    ...){
    if (length(group.cols) > 1) {
        group.col <- group.cols[1]
        group.values <- 
            plot.df[[group.col]] %>% 
            as.factor() %>% 
            droplevels() %>% 
            levels()
        for (group.value in group.values) {
            cat(
                strrep('#', current.header.lvl), group.value, tabset.format,
                nl.delim
            )
            make_tabs_recursive(
                plot.df=plot.df %>% filter(get({{group.col}}) == group.value),
                group.cols=group.cols[2:length(group.cols)],
                current.header.lvl=current.header.lvl + 1,
                plot.fnc=plot.fnc,
                tabset.format=tabset.format,
                nl.delim=nl.delim,
                figure.output.mode=figure.output.mode,
                ...
            )
        }
    } else if (length(group.cols) == 1) {
        plot_figure_tabs(
            plot.df=plot.df, 
            group.col=group.cols[1],
            header.lvl=current.header.lvl,
            plot.fnc=plot.fnc,
            nl.delim=nl.delim,
            figure.output.mode=figure.output.mode,
            ...
        )
    } else if (length(group.cols) == 0) {
        figure <- 
            plot.df %>%
            plot.fnc(...)
        if (figure.output.mode == 'rmd') {
            cat(
                strrep('#', header.lvl),
                nl.delim
            )
            print(figure)
            cat(nl.delim)
        } else if (figure.output.mode == 'return') {
            return(figure)
        } else if (figure.output.mode == 'merged') {
            return(figure)
        } else {
            stop(glue('Invalid arg for figure.output.mode: {figure.output.mode}'))
        }
    }
}

make_nested_plot_tabs <- function(
    plot.df,
    group.cols,
    plot.fnc,
    max.header.lvl=2,
    add.top.layer=FALSE,
    tabset.format="{.tabset}",
    nl.delim="\n\n\n",
    figure.output.mode='rmd',
    ...){
    cat(nl.delim)
    if (add.top.layer) {
        cat(strrep('#', max.header.lvl), tabset.format, nl.delim)
        max.header.lvl <- max.header.lvl + 1
    }
    plot.df %>% 
    make_tabs_recursive(
        group.cols=group.cols,
        current.header.lvl=max.header.lvl,
        plot.fnc=plot.fnc,
        tabset.format=tabset.format,
        nl.delim=nl.delim,
        figure.output.mode=figure.output.mode,
        ...
    )
    cat(nl.delim)
}

###################################################
# Basic Plots
###################################################
plot_barplot <- function(
    plot.df,
    x.var=NULL,
    y.var=NULL,
    fill.var=NULL, 
    label.var=NULL,
    label.color='black',
    label.size=3,
    position='dodge',
    legend.cols=NA,
    ...){
    # Set fill group if specified
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    # make it a boxplot 
    { 
        if (!is.na(legend.cols)) {
            . + 
            geom_col(position=position) +
            guides(fill=guide_legend(ncol=legend.cols))
        } else {
            . + 
            geom_col(position=position)
        }
    } %>% 
    # add text lables to bars (usually pcts)
    {
        if (!(is.null(label.var))) {
            . +
            geom_text(
                aes(label=.data[[label.var]]), 
                position=position_stack(vjust=0.5),
                # position=
                #     ifelse(
                #         position == 'dodge',
                #         position_dodge(),
                #         position_stack(vjust=0.5)
                #     ),
                color=label.color,
                size=label.size,
            )
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_boxplot <- function(
    plot.df,
    x.var=NULL,
    y.var=NULL,
    fill.var=NULL,
    outliers=FALSE,
    outlier.size=1,
    ...){
    # Set fill group if specified
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    { . + geom_boxplot(outliers=TRUE, outlier.size=outlier.size) } %>% 
    # { . + geom_boxplot(outliers=FALSE) } %>% 
    # {
    #     if (outliers) {
    #         {.} + 
    #         geom_jitter(
    #             data=
    #                 .@data %>% 
    #                 filter(
    #                     !!sym(fill.var) < quantile(!!sym(fill.var), 0.25) - 1.5 * IQR(!!sym(fill.var)) | 
    #                     !!sym(fill.var) > quantile(!!sym(fill.var), 0.75) + 1.5 * IQR(!!sym(fill.var))
    #                 ) %>% 
    #                 ungroup(), 
    #             show.legend=FALSE,
    #             size=outlier.size
    #         )
    #     } else {
    #         .
    #     }
    # } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_density <- function(
    plot.df,
    x.var='',
    color.var=NULL, 
    alpha=0.2,
    binwidth=500,
    ...){
    # Set fill group if specified
    {
        # if (!is.null(fill.var)) {
        if (!is.null(color.var)) {
            ggplot(
                plot.df,
                aes(
                    fill=.data[[color.var]],
                    # color=.data[[color.var]],
                    x=.data[[x.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(x=.data[[x.var]])
            )
        }
    } %>% 
    # make it a boxplot 
    { . + geom_density(alpha=alpha) } %>% 
    # { . + geom_density_ridges(alpha=alpha) } %>% 
    # { . + geom_freqpoly(binwidth=binwidth, alpha=alpha) } %>% 
    # { . + geom_histogram(binwidth=binwidth, alpha=alpha) } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_histogram <- function(
    plot.df,
    x.var='',
    fill.var=NULL, 
    alpha=0.7,
    binwidth=0.05,
    ...){
    { 
        ggplot(
            plot.df,
            aes(x=.data[[x.var]])
        ) 
    } %>%
    {
        if (!is.null(fill.var)) {
            . + aes(fill=.data[[fill.var]])
        } else {
            .
        }
    } %>% 
    {
        . +
        geom_histogram(
            position="identity",
            alpha=alpha,
            binwidth=binwidth
        ) 
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_violin <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var=NULL, 
    plot_pts=FALSE,
    jitter.size=0.1,
    quantile.color=NULL,
    draw.quantiles=0L,
    quantile.linewidth=NULL,
    position='dodge',
    adjust=0.5,
    ...){
    # make it a violin plot
    { 
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]],
                fill=.data[[fill.var]],
            )
        )
    } %>%
    {
        .+
        geom_violin(
            quantile.color=quantile.color,
            quantile.linetype=draw.quantiles,
            quantile.linewidth=quantile.linewidth,
            position=position,
            adjust=adjust
        )
    } %>% 
    # plot individual points if specified
    {
        if (plot_pts){
            . +
            geom_jitter(size=jitter.size)
            # geom_jitter(aes(size=jitter.size))
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_heatmap <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var='', 
    label.var=NULL,
    label.size=2,
    label.color='white',
    ...){
    {
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]],
                fill=.data[[fill.var]]
            )
        ) +
        geom_tile()
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        x.scale.mode='discrete',
        y.scale.mode='discrete',
        ...
    ) %>% 
    {
        if (!(is.null(label.var))) {
            . +
            geom_text(
                aes(label=.data[[label.var]]), 
                color=label.color,
                size=label.size,
            )
        } else {
            .
        }
    }
}

plot_jitter <- function(
    plot.df,
    x.var='',
    y.var='',
    color.var=NA, 
    shape.var=NA,
    regression_fnc=NULL,
    add_regression_SE=TRUE,
    alpha=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # Set fill group if specified
    # make it a scatter plot
    { 
        if (!is.na(color.var) & !is.na(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    color=.data[[color.var]],
                    shape=.data[[shape.var]]
                )
            )
        } else if (!is.na(color.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    color=.data[[color.var]]
                )
            )
        } else if (!is.na(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    shape=.data[[shape.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        }
    } %>% 
    {
        . +
        geom_jitter(
            alpha=alpha,
            size=size
        )
    } %>% 
    # add regression line(s) is specified
    {
        if (!is.null(regression_fnc)){
            . +
            geom_smooth(
                method=regression_fnc,
                se=add_regression_SE
            )
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

plot_lineplot <- function(
    plot.df,
    x.var='',
    y.var='',
    group.var='',
    color.var=NULL, 
    linetype.var=NULL,
    alpha=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # Set fill group if specified
    {
        ggplot(
            plot.df,
            aes(
                group=.data[[group.var]],
                x=.data[[x.var]],
                y=.data[[y.var]]
            )
        )
    } %>% 
    {
        if (!is.null(color.var)) {
            . + aes(color=.data[[color.var]])
            # . + aes(color=.data[[color.var]], fill=.data[[color.var]])
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(linetype.var)) {
            . + 
            aes(
                linetype=.data[[linetype.var]],
                shape=.data[[linetype.var]]
            )
        } else {
            .
        }
    } %>% 
    { . + geom_line() } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

plot_pointdensity <- function(
    plot.df,
    x.var='',
    y.var='',
    adjust=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # make it a scatter plot
    { 
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]]
            )
        ) +
        geom_pointdensity(
            adjust=adjust,
            size=size
        ) +
        scale_color_viridis()
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        color.scale.mode='',
        fill.scale.mode='',
        ...
    )
}

plot_contours <- function(
    plot.df,
    x.var='',
    y.var='',
    scales='fixed',
    ...){
    # Set fill group if specified
    {
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]]
            )
        )
    } %>% 
    # make it a scatter plot
    { 
        . + geom_contour_filled(bins=bins)
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

