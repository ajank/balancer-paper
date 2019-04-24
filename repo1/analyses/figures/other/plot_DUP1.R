# Hack for Sascha to have the "good" packages  
# Thanks to Alek for this script!
.libPaths(c(
    "/g/furlong/jankowsk/R-lib/3.4.0-foss-2016b",
    "/g/easybuild/x86_64/CentOS/7/nehalem/software/R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0",
    "/g/easybuild/x86_64/CentOS/7/nehalem/software/R/3.4.0-foss-2016b/lib64/R/library"
))
# END

if (!exists("ggbio_loaded") || ggbio_loaded == FALSE) {
  options(warn = 1)
  library(data.table, quietly = T)
  library(ggbio, quietly = T)
  library(scales, quietly = T) # label=comma
  library(dplyr, quietly = T)
  library(tidyr, quietly = T)
  library(assertthat, quietly = T)

  genome <- "dm6"
  Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"
  Alek.dir <- "/g/furlong/project/33_Hi-C"
  source(paste0(Alek.dir, "/src/R/plot_map_binned_functions.R"))
  ggbio.dir <- "/g/furlong/project/33_Hi-C/Hi-C-ggbio"

  source(paste0(ggbio.dir, "/files.", genome, ".R"))
  source(paste0(ggbio.dir, "/genes.R"))
  source(paste0(ggbio.dir, "/tads.R"))
  ggbio_loaded <- T
}




################################################################################
# Adding bi-allelic frequency here !
if (!exists("bi_freq_loaded") || bi_freq_loaded == FALSE) {
  library(rtracklayer)
  library(AnnotationDbi)
  library(GenomicRanges)
  source("../../common/R/biallelicFreq_v02.R") 
  
  file_snv       <- TabixFile("../../../data/variants/SNVs2/wgs.freebayes-k.filter.norm.vcf.gz")
  file_mapp_new  <- TabixFile("../../../data/dm6_annotation/mappability/dm6_split_90.exact.bedGraph.gz")
  file_cov_vrg   <- TabixFile("../../cov/wgs.VRG.100.GC_norm.bed.gz")
  file_cov_cross <- TabixFile("../../cov/wgs.CROSS.100.GC_norm.bed.gz")
  
  dd = loadDb("../../common/data/dmelanogaster_gene_ensembl.txdb")
  
  bi_freq_loaded <- TRUE
}
################################################################################



### Umbrella theme for the whole plot
theme_bw_no_border <- function (base_size = 12, base_family = "") {
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


#Should I add      vertical_dashed_lines     ??


main_plot <- function(wh) {
  # Data
  gene_data <- subset_gene_data(wh, "embryo")

  #hs = c(6, 6, 6, 2,2,3,2.5, 1.5)
  hs = c(6,3)
  plts <- list(
    #"Hi-C wild-type" = triangle_ggplot("HiC_DB_6-8h_combined_VRG", genome, 5e3, which = wh, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0), zmax = 0.01),
    #"Hi-C balancer" = triangle_ggplot("HiC_DB_6-8h_combined_BAL", genome, 5e3, which = wh, mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0), zmax = 0.01),
    "Hi-C balancer/wild-type" = triangle_ggplot("HiC_DB_6-8h_combined_BAL_vs_VRG_naive", genome, 5e3, which = wh,
        what = "log2FoldChange", colscale = "bluered", #, filter = "significant"
        mark_breakpoints = T, lwd.breakpoints = 2, col.breakpoints = "purple", expand = c(0, 0)),
    #"Mappability"   = svview_plot_mappability(wh, file_mapp_new, normalize=180),
    "Bi-allelic"    = svview_biallelicfreq(wh, file_snv, mapp_bedGraph=file_mapp_new)  # ,
    #"Read coverage" = svview_plot_coverage(wh, cross=file_cov_cross, vrg=file_cov_vrg),
    #"TADs" = plt_TADcalls(wh),
    #"Genes" = autoplot(txdb, exon.rect.h = 0.125, label = F, which = wh)
  )

  message("[main] creating plot")
  title <- paste0(as.character(seqnames(wh))[1], ":", round(start(wh)*1.0/1e6, 1), 
    "-", round(end(wh)*1.0/1e6, 1), "Mb on ", genome)
  t <- tracks(plts, xlim = wh, heights = hs, 
              theme = theme_bw_no_border(), title = title, label.text.cex = 0.8)

  return(t + scale_x_continuous(labels = comma, breaks = scales::pretty_breaks(n = 10), expand = c(0, 0)))
}




DUP <- data.table(chrom = "chr2L", start = 8957e3, end = 9215e3)
SPAN = 350e3

wh <- GRanges("chr2L", IRanges(DUP$start - SPAN, DUP$end + SPAN))
p <- main_plot(wh)

cairo_pdf("plot_DUP1.pdf", width = 10, height = 6)
print(p)
dev.off()
