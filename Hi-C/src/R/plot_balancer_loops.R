options(warn = 1)

require(data.table)
source("src/R/plot_map_binned_functions.R")
source("src/R/functions_balancer_genes.R")

dt <- fread("analysis/balancer/HiC_DB_6-8h_VRG_dm6_annotations.tsv")
dt[, xm := as.integer((x1 + x2) / 2)]
dt[, ym := as.integer((y1 + y2) / 2)]

dt1_gr <- with(dt, GRanges(chr1, IRanges(xm, xm)))
dt[, int1 := findOverlaps(dt1_gr, bp_int_gr, select = "arbitrary")]
stopifnot(!is.na(dt$int1))

dt2_gr <- with(dt, GRanges(chr2, IRanges(ym, ym)))
dt[, int2 := findOverlaps(dt2_gr, bp_int_gr, select = "arbitrary")]
stopifnot(!is.na(dt$int2))


# annotate genes with names of nearby genes

d1 <- 10e3
dt[, loop_id := seq_len(nrow(dt))]
g_all_gr <- GRanges(genes)

for (i in seq_len(nrow(dt)))
{
  int_gr <- GRanges(dt$chr1[i], IRanges(dt$xm[i] - d1 + 1, dt$xm[i] + d1))
  g_gr <- subsetByOverlaps(GRanges(genes[signf != "n"]), int_gr)
  dt$x.within10kb.expr.genes[i] <- paste(g_gr$gene_symbol, collapse = ", ")
  if (length(g_gr) > 0)
    dt$x.anno.gene[i] <- g_gr[nearest(int_gr, g_gr)]$gene_symbol
  else
    dt$x.anno.gene[i] <- g_all_gr[nearest(int_gr, g_all_gr)]$gene_symbol

  int_gr <- GRanges(dt$chr2[i], IRanges(dt$ym[i] - d1 + 1, dt$ym[i] + d1))
  g_gr <- subsetByOverlaps(GRanges(genes[signf != "n"]), int_gr)
  dt$y.within10kb.expr.genes[i] <- paste(g_gr$gene_symbol, collapse = ", ")
  if (length(g_gr) > 0)
    dt$y.anno.gene[i] <- g_gr[nearest(int_gr, g_gr)]$gene_symbol
  else
    dt$y.anno.gene[i] <- g_all_gr[nearest(int_gr, g_all_gr)]$gene_symbol
}

dt_anno <- with(dt[int1 != int2, ], data.table(loop_id, chr1, midpoint1 = xm, chr2, midpoint2 = ym, anno.gene1 = x.anno.gene, anno.gene2 = y.anno.gene,
  within10kb.expr.genes1 = x.within10kb.expr.genes, within10kb.expr.genes2 = y.within10kb.expr.genes))
print(dt_anno)


liftOver_to_dm6bal3 <- function(d)
{
  ftmp1 <- tempfile()
  ftmp2 <- tempfile()
  ftmp3 <- tempfile()
  write.table(d, file = ftmp1, sep = "\t", quote = F, row.names = F, col.names = F)
  system(paste("Rscript src/R/balancer_split_and_liftOver_bed.R", ftmp1, "dm6", "dm6bal3", ftmp2, ftmp3))
  cat("\nShould be only deletions:\n")
  system(paste("cat", ftmp3))
  cat("\n")

  dl <- fread(ftmp2, sep = "\t", header = F)
  unlink(c(ftmp1, ftmp2, ftmp3))
  dl[, V4 := as.numeric(stringr::str_split_fixed(V4, "_", 2)[, 1])]
  if (length(names(d)) == 3)
    dl[, V4 := NULL]
  names(dl) <- names(d)[seq_along(names(dl))]
  return(dl)
}

dt_bal <- cbind(
  liftOver_to_dm6bal3(dt[, c("chr1", "x1", "x2")]),
  liftOver_to_dm6bal3(dt[, c("chr2", "y1", "y2")])
)
dt_bal[, xm := as.integer((x1 + x2) / 2)]
dt_bal[, ym := as.integer((y1 + y2) / 2)]


pdf.dim <- 6
cairo_pdf("analysis/balancer/plot_loops_VRG_vs_BAL.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)

for (i in seq_len(nrow(dt)))
{
  # # plot only the interesting cases of chromatin loops across breakpoints
  # if (dt$int1[i] == dt$int2[i]) next

  try(square_plot("HiC_DB_6-8h_VRG", NA, "dm6", mark_breakpoints = T, zmax = 0.0015,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 250e3, end1 = dt$xm[i] + 250e3, bin_size1 = 5e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 250e3, end2 = dt$ym[i] + 250e3, bin_size2 = 5e3
  ))
  panel.abline(h = current.panel.limits()$ylim[2] - dt$ym[i] + current.panel.limits()$ylim[1])
  panel.abline(v = dt$xm[i])

  try(square_plot("HiC_DB_6-8h_BAL", NA, "dm6bal3", mark_breakpoints = T, zmax = 0.0015,
    chrom1 = dt_bal$chr1[i], start1 = dt_bal$xm[i] - 250e3, end1 = dt_bal$xm[i] + 250e3, bin_size1 = 5e3,
    chrom2 = dt_bal$chr2[i], start2 = dt_bal$ym[i] - 250e3, end2 = dt_bal$ym[i] + 250e3, bin_size2 = 5e3
  ))
  panel.abline(h = current.panel.limits()$ylim[2] - dt_bal$ym[i] + current.panel.limits()$ylim[1])
  panel.abline(v = dt_bal$xm[i])
}

dev.off()


pdf.dim <- 4
cairo_pdf("analysis/balancer/plot_loops_VRG_vs_BAL_across_breakpoint_5kb.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)

for (i in seq_len(nrow(dt)))
{
  # plot only the interesting cases of chromatin loops across breakpoints
  if (dt$int1[i] == dt$int2[i]) next

  try(square_plot("HiC_DB_6-8h_All", NA, "dm6", mark_breakpoints = T, zmax = 0.0015,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 250e3, end1 = dt$xm[i] + 250e3, bin_size1 = 5e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 250e3, end2 = dt$ym[i] + 250e3, bin_size2 = 5e3
  ))

  try(square_plot("HiC_DB_6-8h_VRG", NA, "dm6", mark_breakpoints = T, zmax = 0.0015,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 250e3, end1 = dt$xm[i] + 250e3, bin_size1 = 5e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 250e3, end2 = dt$ym[i] + 250e3, bin_size2 = 5e3
  ))

  try(square_plot("HiC_DB_6-8h_BAL", NA, "dm6bal3", mark_breakpoints = T, zmax = 0.0015,
    chrom1 = dt_bal$chr1[i], start1 = dt_bal$xm[i] - 250e3, end1 = dt_bal$xm[i] + 250e3, bin_size1 = 5e3,
    chrom2 = dt_bal$chr2[i], start2 = dt_bal$ym[i] - 250e3, end2 = dt_bal$ym[i] + 250e3, bin_size2 = 5e3
  ))
}

dev.off()


pdf.dim <- 6
cairo_pdf("analysis/balancer/plot_loops_VRG_vs_BAL_across_breakpoint_20kb_with_diagonal.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)

for (i in seq_len(nrow(dt)))
{
  # plot only the interesting cases of chromatin loops across breakpoints
  if (dt$int1[i] == dt$int2[i]) next
  stopifnot(dt$chr1[i] == dt$chr2[i])

  s <- max((dt$ym[i] - dt$xm[i]) * 0.25, 0.5e6)
  square_plot("HiC_DB_6-8h_VRG", "HiC_DB_6-8h_BAL", "dm6", mark_breakpoints = T, zmax = 0.005,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - s, end1 = dt$ym[i] + s, 20e3
  )
}

dev.off()


pdf.dim <- 4
cairo_pdf("analysis/balancer/plot_loops_VRG_vs_BAL_across_breakpoint_2kb.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)

for (i in seq_len(nrow(dt)))
{
  # plot only the interesting cases of chromatin loops across breakpoints
  if (dt$int1[i] == dt$int2[i]) next

  try(square_plot("HiC_DB_6-8h_All", NA, "dm6", mark_breakpoints = T, zmax = 0.001,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 100e3, end1 = dt$xm[i] + 100e3, bin_size1 = 2e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 100e3, end2 = dt$ym[i] + 100e3, bin_size2 = 2e3
  ))

  try(square_plot("HiC_DB_6-8h_VRG", NA, "dm6", mark_breakpoints = T, zmax = 0.001,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 100e3, end1 = dt$xm[i] + 100e3, bin_size1 = 2e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 100e3, end2 = dt$ym[i] + 100e3, bin_size2 = 2e3
  ))

  try(square_plot("HiC_DB_6-8h_BAL", NA, "dm6bal3", mark_breakpoints = T, zmax = 0.001,
    chrom1 = dt_bal$chr1[i], start1 = dt_bal$xm[i] - 100e3, end1 = dt_bal$xm[i] + 100e3, bin_size1 = 2e3,
    chrom2 = dt_bal$chr2[i], start2 = dt_bal$ym[i] - 100e3, end2 = dt_bal$ym[i] + 100e3, bin_size2 = 2e3
  ))
}

dev.off()

pdf.dim <- 4
cairo_pdf("analysis/balancer/plot_loops_VRG_vs_BAL_across_breakpoint_1kb.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)

for (i in seq_len(nrow(dt)))
{
  # plot only the interesting cases of chromatin loops across breakpoints
  if (dt$int1[i] == dt$int2[i]) next

  try(square_plot("HiC_DB_6-8h_All", NA, "dm6", mark_breakpoints = T, zmax = 0.0005,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 50e3, end1 = dt$xm[i] + 50e3, bin_size1 = 1e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 50e3, end2 = dt$ym[i] + 50e3, bin_size2 = 1e3
  ))

  try(square_plot("HiC_DB_6-8h_VRG", NA, "dm6", mark_breakpoints = T, zmax = 0.0005,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 50e3, end1 = dt$xm[i] + 50e3, bin_size1 = 1e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 50e3, end2 = dt$ym[i] + 50e3, bin_size2 = 1e3
  ))

  try(square_plot("HiC_DB_6-8h_BAL", NA, "dm6bal3", mark_breakpoints = T, zmax = 0.0005,
    chrom1 = dt_bal$chr1[i], start1 = dt_bal$xm[i] - 50e3, end1 = dt_bal$xm[i] + 50e3, bin_size1 = 1e3,
    chrom2 = dt_bal$chr2[i], start2 = dt_bal$ym[i] - 50e3, end2 = dt_bal$ym[i] + 50e3, bin_size2 = 1e3
  ))
}

dev.off()

pdf.dim <- 4
cairo_pdf("analysis/balancer/plot_loops_VRG_vs_BAL_across_breakpoint_20kb.pdf", width = pdf.dim, height = pdf.dim - 1, onefile = T)

for (i in seq_len(nrow(dt)))
{
  # plot only the interesting cases of chromatin loops across breakpoints
  if (dt$int1[i] == dt$int2[i]) next

  try(square_plot("HiC_DB_6-8h_VRG", NA, "dm6", mark_breakpoints = T, zmax = 0.001,
    chrom1 = dt$chr1[i], start1 = dt$xm[i] - 0.5e6, end1 = dt$xm[i] + 0.5e6, bin_size1 = 20e3,
    chrom2 = dt$chr2[i], start2 = dt$ym[i] - 0.5e6, end2 = dt$ym[i] + 0.5e6, bin_size2 = 20e3
  ))

  try(square_plot("HiC_DB_6-8h_BAL", NA, "dm6bal3", mark_breakpoints = T, zmax = 0.001,
    chrom1 = dt_bal$chr1[i], start1 = dt_bal$xm[i] - 0.5e6, end1 = dt_bal$xm[i] + 0.5e6, bin_size1 = 20e3,
    chrom2 = dt_bal$chr2[i], start2 = dt_bal$ym[i] - 0.5e6, end2 = dt_bal$ym[i] + 0.5e6, bin_size2 = 20e3
  ))
}

dev.off()


annotate_dist <- function(dt)
{
  dt[, dist := NA]

  sel <- dt$chr1 == dt$chr2
  dt$dist[sel] <- (dt$ym - dt$xm)[sel]

  sel <- dt$chr1 == "chr2L" & dt$chr2 == "chr2R"
  dt$dist[sel] <- (with(dm6bal3_chrom_sizes, chrom_size[match("chr2L", chrom)]) - dt$xm + dt$ym)[sel]

  sel <- dt$chr1 == "chr2R" & dt$chr2 == "chr2L"
  dt$dist[sel] <- (with(dm6bal3_chrom_sizes, chrom_size[match("chr2L", chrom)]) - dt$ym + dt$xm)[sel]

  sel <- dt$chr1 == "chr3L" & dt$chr2 == "chr3R"
  dt$dist[sel] <- (with(dm6bal3_chrom_sizes, chrom_size[match("chr3L", chrom)]) - dt$xm + dt$ym)[sel]

  sel <- dt$chr1 == "chr3R" & dt$chr2 == "chr3L"
  dt$dist[sel] <- (with(dm6bal3_chrom_sizes, chrom_size[match("chr3L", chrom)]) - dt$ym + dt$xm)[sel]

  return(dt)
}

print(annotate_dist(dt)[int1 != int2, ])
print(annotate_dist(dt_bal)[dt$int1 != dt$int2, ])
