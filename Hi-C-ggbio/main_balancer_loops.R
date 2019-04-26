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


### Annotation tracks

DHS_dt <- fread("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/annotations/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed", header = F, sep = "\t", col.names = c("chrom", "start", "end"))
DHS_dt <- DHS_dt[grepl("^chr[23]", chrom), ]
DHS_dt[, start := start + 1L]
DHS_gr <- GRanges(DHS_dt)
DHS_gr$dhs_id <- seq_along(DHS_gr)

# combine insulator binding into clusters
insulator.dir <- "/g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData"
insulator.files <- c(
  "CTCF" = "ChipSeq.CTCF.embryo_stage2.Boyle_2014.bed",
  "su(Hw)" = "ChipChip.SuHw.embryo_0-12h.White_K_27.bed",
  "Zw5" = "ChipChip.ZW5.embryo_2-4h.Karpen_G_5265.bed",
  "BEAF-32" = "ChipChip.BEAF-32.embryo_0-12h.White_K_21.bed",
  "Cp190" = "ChipChip.CP190.embryo_0-12h.White_K_22.bed"
)

insulator_dt <- NULL
for (n in names(insulator.files))
{
  dt <- fread(paste0(insulator.dir, "/", insulator.files[n]), header = F, sep = "\t")
  dt <- dt[, 1:3]
  names(dt) <- c("chrom", "start", "end")
  dt[, start := start + 1L]
  dt <- dt[grepl("^chr[23]", chrom), ]
  dt$source <- n
  insulator_dt <- rbind(insulator_dt, dt)
}

insulator_gr <- GRanges(insulator_dt)

unmetadata <- function(gr)
{
  elementMetadata(gr) <- NULL
  return(gr)
}

CAD4_dt <- fread("/g/furlong/project/CAD/CAD4_corrected_names/CAD4_plus_vienna_minus_inactive_corrected_names_dm6.bed",
header = F, sep = "\t", col.names = c("chrom", "start", "end", "name", "score", "strand"))
CAD4_dt <- CAD4_dt[grepl("^chr[23]", chrom), ]
CAD4_dt[, start := start + 1L]
CAD4_gr <- GRanges(CAD4_dt)

# 8008 CRMs, Zinzen et al. 2009
mesoCRM_gr <- rtracklayer::import("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/annotations/MesoCRM_dm6_Nature_Zinzen2009.gff")

enhancer_gr <- c(unmetadata(CAD4_gr), unmetadata(mesoCRM_gr))
enhancer_dt <- as.data.table(enhancer_gr)
setnames(enhancer_dt, "seqnames", "chrom")


SV.dt <- rbind(
    data.table(DHS_dt, name = "DNase I hypersensitive sites"),
    # data.table(enhancer_dt, name = "Enhancers"),
    # data.table(insulator_dt, name = "Insulators"),
    fill = T
)
SV.dt[, start := start - 1L]
SV.dt[, name := factor(name)]


plt_SV_enhancers <- function(wh) 
{
    # Manual color assignment:
    bed.scale_cols <- c("DNase I hypersensitive sites" = "#6f6f6f",
                        "Enhancers" = "#984EA3",
                        "Insulators" = "#984EA3")
    SV.dt.wh <- SV.dt[chrom == as.character(seqnames(wh))[1] & start >= start(wh)[1] & end <= end(wh)[1], ]
    # SV.ovlps.wh <- SV.ovlps[chrom == as.character(seqnames(wh))[1] & start >= start(wh)[1] & end <= end(wh)[1], ]
    pSV <- ggplot(SV.dt.wh) + 
        aes(x = start, xend = end, y = name, yend = name, col = name) +
        geom_segment(size=2) +
        guides(color = "none") +
        ylab(NULL) + 
        scale_y_discrete(position = "right") +
        scale_color_manual(values = bed.scale_cols) + 
        # geom_point(data = SV.ovlps.wh,
        #            aes( x = (end+start)/2, y = name, col = del),
        #            inherit.aes = F, size=2, col = "darkorange", alpha = 0.8) +
        coord_cartesian(expand = 0)
    return(pSV)
}


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

    hs = c(0.75, 1.75, 1.75, 3.5, 2)
    plts <- list(
                if (genome == "dm6") plt_SV_enhancers(wh) else NULL,
                "Genes" = autoplot(txdb, exon.rect.h = 0.125, label = F, which = wh),
                plot_gene2column_mapping(gene_data) + scale_bg_colors, # +
                        # ggplot_arbitrary_vertical_lines(wh, breakpoints, df = gene_data, ymax = 6.25, gp = breakpoint.gpar) +
                        # ggplot_arbitrary_vertical_lines(wh, vlines.vrg, df = gene_data, ymax = 3.5, gp = vlines.gpar),
                "DE" = plot_lfc(gene_data, c(-3, 3)) + scale_bg_colors + 
                        guides(fill = 'none'),
                "Expr." = plot_gene_expression(gene_data, ylims = c(0.01, 300)) +
                        scale_bg_colors)

    message("[main] creating plot")
    title <- paste0(as.character(seqnames(wh))[1], ":", round(start(wh)*1.0/1e6, 1), 
                    "-", round(end(wh)*1.0/1e6, 1), "Mb on ", genome, " in ", dataset, "s")
    t <- tracks(plts, xlim = wh, heights = hs, theme = theme_bw_no_border(), title = NULL, label.text.cex = 0.8)

    return(t + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 3), expand = c(0, 0)))
}


