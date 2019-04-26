options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 1)
  stop("Usage:  Rscript extract_DESeq2_tracks_balancer_cap2.R outputDir")

outputDir <- options[1]

require(Chicago)
require(GenomicRanges)
require(rtracklayer)
options(warn = 1)
eff_genome_size <- 130e6

#
#  read CHiCAGO data
#

load(paste0(outputDir, "/all_interactions.Rdata")) # pm_all
load(paste0(outputDir, "/DESeq2_interactions.Rdata")) # baits, pm, pmd

suppressWarnings(dir.create(paste0(outputDir, "/tracks")))
suppressWarnings(dir.create(paste0(outputDir, "/igv_sessions")))
suppressWarnings(dir.create(paste0(outputDir, "/igv_sessions_differential")))
suppressWarnings(dir.create(paste0(outputDir, "/igv_sessions_differential_CAD4_or_MesoCRM")))
suppressWarnings(dir.create(paste0(outputDir, "/igv_sessions_differential_DHS_or_CAD4_or_MesoCRM")))

#
#  save bedGraph tracks
#

process <- function(d, suffix)
{
  for (i in seq_len(nrow(baits)))
  {
    name <- baits$name[i]
    ID <- baits$baitID[i]
    track <- d[baitID == ID, c("otherEndChr", "otherEndStart", "otherEndEnd", "N", "N.1", "N.2"), with = F]
    track[, otherEndStart := otherEndStart - 1L]
    track[, N := N * eff_genome_size / sum(N.1 + N.2) / 1e3]
    track[, N.1 := NULL]
    track[, N.2 := NULL]
    write.table(track, file = gzfile(paste0(outputDir, "/tracks/", name, "_", suffix, ".RPGC.bedGraph.gz")),
      row.names = F, col.names = F, quote = F, sep = "\t")
  }
}

pm_all[, N := N.1.VRG + N.2.VRG]
process(pm_all, "VRG")
pm_all[, N := N.1.BAL + N.2.BAL]
process(pm_all, "BAL")
pm_all[, N := N.1 + N.2]
process(pm_all, "All")

#
#  save GFF files with called interactions, and IGV tracks
#

export_gff <- function(dt, filename)
{
  setnames(dt, "otherEndChr", "chrom")
  setnames(dt, "otherEndStart", "start")
  setnames(dt, "otherEndEnd", "end")
  removed_cols <- intersect(names(dt), c("baitID", "otherEndID", "s_j", "s_i", "s_j.VRG", "s_i.VRG", "s_j.BAL", "s_i.BAL", "baitChr", "baitStart", "baitEnd", "baitName", "distance", "value", "distance_other", "value_other"))
  dt[, (removed_cols) := NULL]
  export.gff3(if (nrow(dt) > 0) GRanges(dt) else GRanges(), filename)
}

export_bedGraph <- function(dt, filename, score.column = NULL)
{
  dte <- data.table(chrom = dt$otherEndChr, start = dt$otherEndStart - 1L, end = dt$otherEndEnd)
  dte$score <- dt[, score.column, with = F]
  write.table(dte, file = filename, row.names = F, col.names = F, quote = F, sep = "\t")
}

igv <- readLines("analysis/balancer_cap2/igv_session_template.xml")
igv <- gsub('="/Users/ajank/balancer_cap2/igv_sessions/', '="../../', igv)
# igv <- gsub('="../../contacts_all/', '="../', igv)
igv <- gsub('="../../igv_session.xml"', '="igv_session.xml"', igv)

#
#  save IGV sessions for all the viewpoints (also the ones that are not in "baits")
#

load("analysis/balancer_cap2/viewpoints.Rdata") # lib

pm$color <- "204,204,204" # interaction not differential
pm$color[pm$int.padj < 0.05 & pm$int.log2FoldChange < log2(1.5)] <- "77,175,74" # differential interaction stronger in wild-type
pm$color[pm$int.padj < 0.05 & pm$int.log2FoldChange > log2(1.5)] <- "55,126,184" # differential interaction stronger in balancer

for (i in seq_along(lib))
{
  name <- names(lib)[i]
  ID <- lib$frag_id[i]
  # export_gff(pm_all[baitID == ID & score > 5], paste0(outputDir, "/tracks/", name, "_All_interactions.gff"))
  # export_bedGraph(pm_all[baitID == ID & score > 0],
  #   gzfile(paste0(outputDir, "/tracks/", name, "_All_interactions_score.bedGraph.gz")), "score")
  # export_bedGraph(pm_all[baitID == ID],
  #   gzfile(paste0(outputDir, "/tracks/", name, "_All_interactions_log.p.bedGraph.gz")), "log.p")

  export_gff(pm[baitID == ID, ], paste0(outputDir, "/tracks/", name, "_differential.gff"))
  # export_gff(pmd[baitID == ID & int.log2FoldChange > 0, ], paste0(outputDir, "/tracks/", name, "_differential_BAL.gff"))
  # export_gff(pmd[baitID == ID & int.log2FoldChange < 0, ], paste0(outputDir, "/tracks/", name, "_differential_VRG.gff"))

  chrom <- as.character(seqnames(lib)[i])
  mid <- as.integer((start(lib)[i] + end(lib)[i]) / 2)
  this_igv <- gsub(' locus="[^"]*"', paste0(' locus="', chrom, ':', mid - 50000L, '-', mid + 50000L, '"'),
    gsub("Glut1", name, igv))
  writeLines(this_igv, paste0(outputDir, "/igv_sessions/", name, "_igv_session.xml"))
  if (any(pmd$baitID == ID))
    writeLines(this_igv, paste0(outputDir, "/igv_sessions_differential/", name, "_igv_session.xml"))
  if (any(pmd$baitID == ID & (pmd$oeCAD4 | pmd$oeMesoCRM)))
    writeLines(this_igv, paste0(outputDir, "/igv_sessions_differential_CAD4_or_MesoCRM/", name, "_igv_session.xml"))
  if (any(pmd$baitID == ID & (pmd$oeDHS | pmd$oeCAD4 | pmd$oeMesoCRM)))
    writeLines(this_igv, paste0(outputDir, "/igv_sessions_differential_DHS_or_CAD4_or_MesoCRM/", name, "_igv_session.xml"))
}
