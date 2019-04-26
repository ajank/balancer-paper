source("src/R/plot_map_binned_functions.R")

pdf.dim <- 10
pdf("analysis/balancer/plot_balancer_genome.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)
square_plot("HiC_DB_6-8h_VRG", "HiC_DB_6-8h_BAL", "dm6", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
square_plot("HiC_DB_6-8h_VRG", "HiC_DB_6-8h_BAL", "dm6bal3", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
dev.off()


source("src/R/plot_map_binned_functions.R")
pdf.dim <- 10
for (bin_size in as.integer(c(200e3, 100e3)))
{
  chrom_sep_bins <- 1e6 / bin_size
  breakpoint_summarize <- "mean"
  for (dataset in c("HiC_DB_6-8h_BAL", "HiC_DB_6-8h_VRGdownBAL"))
  {
    ensure_loaded(dataset, "dm6", bin_size)
    c <- chrom_maps[[paste0(dataset, "_dm6")]][[as.character(bin_size)]]
    v1 <- c$map_start[match("chr2L", c$chrom)]:c$map_end[match("chr2R", c$chrom)]
    v2 <- c$map_start[match("chr3L", c$chrom)]:c$map_end[match("chr3R", c$chrom)]
    maps[[paste0(dataset, "_dm6")]][[as.character(bin_size)]]$genome$observed[v1, v2] <- NA
  }
  pdf(paste0("analysis/balancer/plot_balancer_genome_VRGdownBAL_", bin_size / 1e3, "kb.pdf"), width = pdf.dim, height = pdf.dim - 1, onefile = T)
  square_plot("HiC_DB_6-8h_BAL", "HiC_DB_6-8h_VRGdownBAL", "dm6", "genome", NA, NA, bin_size, mark_breakpoints = T, background = "gray",
    zmax = ifelse(bin_size == 200e3, 0.01, 0.005))
  # square_plot("HiC_DB_6-8h_BAL", "HiC_DB_6-8h_VRGdownBAL", "dm6", "chr2L", 1, 3e6, bin_size, mark_breakpoints = T, background = "gray", zmax = 0.01)
  dev.off()
  system(paste0("ls -lh analysis/balancer/plot_balancer_genome_VRGdownBAL_", bin_size / 1e3, "kb.pdf"))
}
source("src/R/plot_map_binned_functions.R")


pdf.dim <- 10
pdf("analysis/balancer/plot_balancer_genome_VRGdownBAL.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)
square_plot("HiC_DB_6-8h_VRGdownBAL", "HiC_DB_6-8h_BAL", "dm6", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
square_plot("HiC_DB_6-8h_VRGdownBAL", "HiC_DB_6-8h_BAL", "dm6bal3", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
dev.off()


pdf.dim <- 14
cairo_pdf("analysis/balancer/plot_balancer_genome_triangle.pdf", width = pdf.dim, height = pdf.dim / 2 - 0.1, onefile = T)
triangle_plot("HiC_DB_6-8h_VRG", "dm6", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
# triangle_plot("HiC_DB_6-8h_VRGdownBAL", "dm6", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
triangle_plot("HiC_DB_6-8h_BAL", "dm6", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)

triangle_plot("HiC_DB_6-8h_VRG", "dm6bal3", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
# triangle_plot("HiC_DB_6-8h_VRGdownBAL", "dm6bal3", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
triangle_plot("HiC_DB_6-8h_BAL", "dm6bal3", "genome", NA, NA, 50e3, mark_breakpoints = F, background = "gray", zmax = 0.003)
dev.off()


pdf.dim <- 14
cairo_pdf("analysis/balancer/plot_balancer_genome_DESeq2.pdf", width = pdf.dim, height = pdf.dim / 2 - 0.1, onefile = T)
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG_naive", "dm6", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG", "dm6", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG", "dm6", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  filter = "significant", mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint", "dm6", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint", "dm6", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  filter = "significant", mark_breakpoints = F, background = "gray")


triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG_naive", "dm6bal3", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG", "dm6bal3", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG", "dm6bal3", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  filter = "significant", mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint", "dm6bal3", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  mark_breakpoints = F, background = "gray")
triangle_plot("HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint", "dm6bal3", "genome", NA, NA, 50e3, what = "log2FoldChange", colscale = "bluered",
  filter = "significant", mark_breakpoints = F, background = "gray")
dev.off()
