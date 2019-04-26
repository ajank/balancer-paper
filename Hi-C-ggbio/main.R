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
main_plot <- function(wh, dataset = "embryo", genome = "dm6", add_lines = numeric()) 
{
    assert_that(is.character(dataset) && length(dataset) == 1)
    assert_that(dataset == "embryo" || dataset == "adult")
    
    # Data
    gene_data <- subset_gene_data(wh, dataset)
    
    # Vertical lines
    breakpoints <- extract_breakpoints(genome, as.character(seqnames(wh)[1]), start(wh), end(wh), bin_size = 5e3, summarize = "mean_if_range_lt", range_thr = 1e3)
    vlines <- add_lines[start(wh) <= add_lines & add_lines <= end(wh)]

    # Graphical parameters
    breakpoint.gpar <- gpar(lty = "dashed", lwd = 2, lineend = "butt", col = "purple")
    vlines.gpar <- gpar(lty = "dotted", lwd = 2, lineend = "butt", col = "purple")

    hs = c(3, 6, 5, 6, 6, 6, 2.5, 1.5, 1.5, 1.5, 1.5, 2, 2.5, 2, 2)
    plts <- list(
                "Region" = ggplot(wh) + 
                        annotation_custom(grob = ggplotGrob(plot_overview(wh, genome))),
                "Hi-C all reads" = triangle_ggplot("HiC_DB_6-8h_combined_All", genome, 5e3, which = wh, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, df = NULL, gp = vlines.gpar),
                "TADs" = plt_TADcalls(wh) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = breakpoints, gp = breakpoint.gpar) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, ymax = 4.5, gp = vlines.gpar),
                "Hi-C wild\u00adtype" = triangle_ggplot("HiC_DB_6-8h_combined_VRG", genome, 5e3, which = wh, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)) + 
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "Hi-C balancer" = triangle_ggplot("HiC_DB_6-8h_combined_BAL", genome, 5e3, which = wh, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)) + 
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "Hi-C balancer/wild\u00adtype" = triangle_ggplot("HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint_dist_100kb", genome, 5e3, which = wh,
                    what = "log2FoldChange", colscale = "bluered", # , filter = "significant"
                    mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)) + 
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "SVs/Enhancers" = plt_SV_enhancers(wh) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = breakpoints, gp = breakpoint.gpar) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "DHS" = plt_DHS_track(wh) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = breakpoints, gp = breakpoint.gpar) +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "RNA\nall" = plot_RNA_track(wh, dataset, "all") +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = breakpoints, gp = breakpoint.gpar) + 
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "RNA\nwt" = plot_RNA_track(wh, dataset, "vrg") +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = breakpoints, gp = breakpoint.gpar) + 
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "RNA\nbal" = plot_RNA_track(wh, dataset, "bal") +
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = breakpoints, gp = breakpoint.gpar) + 
                        ggplot_arbitrary_vertical_lines(which = wh, xintercepts = vlines, gp = vlines.gpar),
                "Genes" = autoplot(txdb, exon.rect.h = 0.125, label = F, which = wh),
                plot_gene2column_mapping(gene_data) + scale_bg_colors +
                        ggplot_arbitrary_vertical_lines(wh, breakpoints, df = gene_data, ymax = 1.8, gp = breakpoint.gpar) +
                        ggplot_arbitrary_vertical_lines(wh, vlines, df = gene_data, ymax = 1.8, gp = vlines.gpar),
                "DE" = plot_lfc(gene_data, c(-3, 3)) + scale_bg_colors + 
                        guides(fill = 'none'),
                "Expr." = plot_gene_expression(gene_data, ylims = c(0.01, 300)) +
                        scale_bg_colors )

    message("[main] creating plot")
    title <- paste0(as.character(seqnames(wh))[1], ":", round(start(wh)*1.0/1e6, 1), 
                    "-", round(end(wh)*1.0/1e6, 1), "Mb on ", genome, " in ", dataset, "s")
    t <- tracks(plts, xlim = wh, heights = hs, theme = theme_bw_no_border(), title = title, label.text.cex = 0.8)

    return(t + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10), expand = c(0, 0)))
}


