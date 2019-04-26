assert_that(exists("genome"))


### Umbrella theme for the whole plot
theme_bw_no_border <- function (base_size = 12, base_family = "") 
{
    theme(
        # axis.text = element_text(size = rel(0.8)),
        # axis.ticks = element_line(colour = "black"), 
        # legend.key = element_rect(colour = "grey80"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_line(colour = "grey98", size = 0.5),
        strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
}

### Colors
scale_bg_colors = scale_fill_manual(values = c("1" = "#E9E8CA", 
                                               "2" = "#E3E7D5", 
                                               "3" = "#DCE5DF", 
                                               "4" = "#D6E4EA", 
                                               "5" = "#CEE3F4", 
                                              "insignificant" = "dodgerblue4", 
                                              "significant" = "darkorange"))


#' Given GRanges region, subset global data and combine all the plots into one.
#' Dataset can be "embryo", "adult", or a vector with both in it.
#'
#' A lot of global data needs to be in the environment!
#' 
main_plot <- function(wh, dataset = "embryo", genome = "dm6") 
{
    assert_that(is.character(dataset) && length(dataset) == 1)
    assert_that(dataset == "embryo" || dataset == "adult")
    
    # Data
    gene_data <- subset_gene_data(wh, dataset)
    
    # Vertical lines
    breakpoints <- extract_breakpoints(genome, as.character(seqnames(wh)[1]), start(wh), end(wh), bin_size = 5e3, summarize = "mean_if_range_lt", range_thr = 1e3)
    vlines.vrg <- start(subsetByOverlaps(tad.lines.vrg, wh, minoverlap = 0L))
    vlines.vrg.relaxed <- start(subsetByOverlaps(tad.lines.vrg.relaxed, wh, minoverlap = 0L))
    vlines.bal <- start(subsetByOverlaps(tad.lines.bal, wh, minoverlap = 0L))
    vlines.bal.relaxed <- start(subsetByOverlaps(tad.lines.bal.relaxed, wh, minoverlap = 0L))

    # Graphical parameters
    breakpoint.gpar <- gpar(lty = "dashed", lwd = 2, lineend = "butt", col = "purple")
    vlines.gpar <- gpar(lty = "dotted", lwd = 2, lineend = "butt", col = "dodgerblue4")
    vlines.gpar.relaxed <- gpar(lty = "dotted", lwd = 2, lineend = "butt", col = "#6f6f6f", alpha = 0.6)

    hs = c(3, 4.5, 5.5, 5.5, 5.5, 2.5, 2.5, 2.5, 2, 2, 2.5, 2)
    plts <- list(
                "Region" = ggplot(wh) + 
                        annotation_custom(grob = ggplotGrob(plot_overview(wh, genome))),
                "TADs" = plt_TADcalls(wh),
                "Hi-C wild\u00adtype" = triangle_ggplot("HiC_DB_6-8h_combined_VRG", genome, 5e3, which = wh, use_mult = F, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)),
                "Hi-C balancer" = triangle_ggplot("HiC_DB_6-8h_combined_BAL", genome, 5e3, which = wh, use_mult = F, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)),
                "Hi-C bal/wt" = triangle_ggplot("HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint", genome, 5e3, which = wh,
                    what = "log2FoldChange", colscale = "bluered", # , filter = "significant" _dist_100kb
                    mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)),
                "IS\nwt" = plot_IS_track(wh, "VRG") +
                        ggplot_arbitrary_vertical_lines(wh, vlines.vrg.relaxed, gp = vlines.gpar.relaxed) +
                        ggplot_arbitrary_vertical_lines(wh, vlines.vrg, gp = vlines.gpar),
                "IS\nbal" = plot_IS_track(wh, "BAL") +
                        ggplot_arbitrary_vertical_lines(wh, vlines.bal.relaxed, gp = vlines.gpar.relaxed) +
                        ggplot_arbitrary_vertical_lines(wh, vlines.bal, gp = vlines.gpar),
                # "PC1\nwt" = plot_PC1_track(wh, "VRG"),
                # "PC1\nbal" = plot_PC1_track(wh, "BAL"),
                "SV/Enh." = plt_SV_enhancers(wh),
                "Genes" = autoplot(txdb, exon.rect.h = 0.125, label = F, which = wh),
                plot_gene2column_mapping(gene_data) + scale_bg_colors +
                        ggplot_arbitrary_vertical_lines(wh, breakpoints, df = gene_data, ymax = 6.25, gp = breakpoint.gpar) +
                        ggplot_arbitrary_vertical_lines(wh, vlines.vrg, df = gene_data, ymax = 3.5, gp = vlines.gpar),
                "DE" = plot_lfc(gene_data, c(-3, 3)) + scale_bg_colors + 
                        guides(fill = 'none'),
                "Expr." = plot_gene_expression(gene_data, ylims = c(0.01, 300)) +
                        scale_bg_colors)

    message("[main] creating plot")
    title <- paste0(as.character(seqnames(wh))[1], ":", round(start(wh)*1.0/1e6, 1), 
                    "-", round(end(wh)*1.0/1e6, 1), "Mb on ", genome, " in ", dataset, "s")
    t <- tracks(plts, xlim = wh, heights = hs, theme = theme_bw_no_border(), title = title, label.text.cex = 0.8)

    return(t + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10), expand = c(0, 0)))
}


