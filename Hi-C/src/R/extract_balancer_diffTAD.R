options(warn = 1)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_IS.R")
source("src/R/logFC_vs_dist.R")

# read TAD boundaries
# b.VRG <- read_boundaries("HiC_DB_6-8h_combined_VRG", only23 = T)
b.VRG.relaxed <- read_boundaries("HiC_DB_6-8h_combined_VRG", domain_source = "filtered_5000_thresholdComparisons_0.1", only23 = T)
b.VRG.strict <- read_boundaries("HiC_DB_6-8h_combined_VRG", domain_source = "filtered_5000", only23 = T)
# b.VRG.strict <- b.VRG.strict[b.VRG.strict$delta > 0.1]
# b.BAL <- read_boundaries("HiC_DB_6-8h_combined_BAL", genome = "dm6bal3", only23 = T)
b.BAL.relaxed <- read_boundaries("HiC_DB_6-8h_combined_BAL", genome = "dm6bal3", domain_source = "filtered_5000_thresholdComparisons_0.1", only23 = T)
b.BAL.strict <- read_boundaries("HiC_DB_6-8h_combined_BAL", genome = "dm6bal3", domain_source = "filtered_5000", only23 = T)
# b.BAL.strict <- b.BAL.strict[b.BAL.strict$delta > 0.1]


# match TAD boundaries between datasets
annotate_closest <- function(gr1, gr2)
{
  gr1$midpoint <- start(gr1)
  d <- distanceToNearest(gr1, gr2)
  stopifnot(!anyDuplicated(queryHits(d)))
  gr1$match <- NA_integer_
  gr1$match[queryHits(d)] <- subjectHits(d)
  gr1$distance <- NA_integer_
  gr1$distance[queryHits(d)] <- elementMetadata(d)$distance
  gr1$match_midpoint <- gr2$midpoint[gr1$match]
  return(gr1)
}

# split TAD boundaries into categories
classify_domains <- function(gr1, gr2.strict, gr2.relaxed, match_distance = 25e3, shift_distance = 25e3)
{
  gr1.strict <- annotate_closest(gr1, gr2.strict)
  gr1.relaxed <- annotate_closest(gr1, gr2.relaxed)

  gr1$class_detailed <- factor(
    ifelse(!is.na(gr1.strict$distance) & gr1.strict$distance <= match_distance, "matched (to strict)",
      ifelse(!is.na(gr1.relaxed$distance) & gr1.relaxed$distance <= match_distance, "matched (to relaxed)",
        ifelse(!is.na(gr1.strict$distance) & gr1.strict$distance <= shift_distance, "shifted (to strict)",
          ifelse(!is.na(gr1.relaxed$distance) & gr1.relaxed$distance <= shift_distance, "shifted (to relaxed)",
            "deleted"
          )
        )
      )
    ),
  c("matched (to strict)", "matched (to relaxed)", "shifted (to strict)", "shifted (to relaxed)", "deleted"))

  gr1$class <- ifelse(gr1$class_detailed %in% c("matched (to strict)", "matched (to relaxed)"), "matched",
    ifelse(gr1$class_detailed %in% c("shifted (to strict)", "shifted (to relaxed)"), "shifted", "deleted"))
  return(as.data.table(gr1))
}

b.VRG_to_BAL <- classify_domains(GRanges(b.VRG.strict), GRanges(b.BAL.strict), GRanges(b.BAL.relaxed))
b.VRG_to_BAL$class[b.VRG_to_BAL$class == "deleted"] <- "wild-type-specific"
# print(summary(b.VRG_to_BAL))

b.BAL_to_VRG <- classify_domains(GRanges(b.BAL.strict), GRanges(b.VRG.strict), GRanges(b.VRG.relaxed))
b.BAL_to_VRG$class[b.BAL_to_VRG$class == "deleted"] <- "balancer-specific"
# print(summary(b.BAL_to_VRG))

b.strict <- rbind(data.table(b.VRG_to_BAL, allele = "w"), data.table(b.BAL_to_VRG, allele = "b"))
b.strict$class <- factor(b.strict$class, c("matched", "shifted", "wild-type-specific", "balancer-specific"))
names(b.strict)[names(b.strict) == "seqnames"] <- "chrom"

# take the matched boundaries in wild type as scaffolding
b.strict <- b.strict[!(allele == "b" & class == "matched")]
ext <- b.strict[allele == "w" & class == "matched"]
ext[, allele := "b"]
# ext[, delta := NA]
# ext[, pvalue := NA]
ext[, tad_sep := NA]
ext[, tad_sep_scaled := NA]
ext[, value_scaled := NA]
ext[, prominence := NA]
b.strict <- rbind(b.strict, ext)
setkey(b.strict, chrom, start, allele)

# recalculate prominence, take allele-specific only with pronounced prominence
bb1 <- b.strict[allele == "b"]
# convert to 0-based coordinates
bb1[, start := start - 1L]

bb1[, start := start - 1L]
bb1[, end := end + 1L]
bb1[, width_left := NULL]
bb1[, width_right := NULL]
bb1[, width_smaller := NULL]
bb1[, width_bigger := NULL]
bb2 <- liftOver_to_dm6bal3(bb1)
setkey(bb2, chrom, start)
bb3 <- annotate_boundaries(bb2, "HiC_DB_6-8h_combined_BAL", genome = "dm6bal3", domain_source = "filtered_5000")
b.strict <- rbind(b.strict[allele == "w"], bb3)

setkey(b.strict, chrom, start, allele)
write.table(b.strict, file = "analysis/balancer/differential_TADs.tab", sep = "\t", quote = F, row.names = F, col.names = T)
print(with(b.strict, table(class, allele)))


genes_to_VRG <- logFC_vs_dist(genes, GRanges(b.VRG_to_BAL), bp_gr, mode = "closest", max_dist = Inf)
genes_to_VRG$class <- b.VRG_to_BAL$class[genes_to_VRG$center_id]
# with(genes_to_VRG[signf == "s"], wilcox.test(abs(log2FoldChange[class == "matched"]), abs(log2FoldChange[class == "wild-type-specific"])))
# fisher.test(with(genes_to_VRG[signf != "n"], table(as.character(signf), sign(log2FoldChange))))$p.value
# fisher.test(with(genes_to_VRG[signf != "n" & class %in% c("matched", "wild-type-specific")], table(as.character(class), sign(log2FoldChange))))$p.value

genes_to_BAL <- logFC_vs_dist(genes, GRanges(b.BAL_to_VRG), bp_gr, mode = "closest", max_dist = Inf)
genes_to_BAL$class <- b.BAL_to_VRG$class[genes_to_BAL$center_id]
# with(genes_to_BAL[signf == "s"], wilcox.test(abs(log2FoldChange[class == "matched"]), abs(log2FoldChange[class == "balancer-specific"])))

# genes_to_strict <- rbind(genes_to_VRG, genes_to_BAL[class == "balancer-specific", ])
genes_to_strict <- logFC_vs_dist(genes, GRanges(b.strict), bp_gr, mode = "closest", max_dist = Inf)
genes_to_strict$class <- b.strict$class[genes_to_strict$center_id]
genes_to_strict[, class_simplified := factor(ifelse(grepl("specific", class), "allele\u00adspecific\nboundaries", "matched\nboundaries"),
  c("matched\nboundaries", "allele\u00adspecific\nboundaries"))]
print(fisher.test(with(genes_to_strict[signf != "n"], table(signf, class_simplified))))
