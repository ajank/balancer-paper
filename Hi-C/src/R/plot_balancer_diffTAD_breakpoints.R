options(warn = 1)

require(cowplot)
require(ggplot2)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))
theme_update(
  strip.text.x = element_text(margin = margin(t = 11 / 2, b = 11 / 2), size = 11),
  strip.text.y = element_text(margin = margin(l = 11 / 2, r = 11 / 2), size = 11)
)


levels_id_trimmed <- c("chr2_1", "chr3_1", "chr3_1a", "chr3_1b", "chr3_6", "chr3_Rend",
  "chr2_2", "chr2_3", "chr2_4", "chr2_5", "chr2_6", "chr3_2", "chr3_3", "chr3_4", "chr3_5", "chr3_7", "chr3_8", "chr3_9")


suffix <- "_thresh12.5kb"

# add pseudo-boundaries at ends of chromosomes
chrom_sel <- unique(diffTADs_dt$chrom)
chrom_len <- chrom_map$length[match(chrom_sel, chrom_map$chrom)]
d1 <- data.table(chrom = chrom_sel, start = 1L, end = 0L, midpoint = 1L, class = "matched", delta = 1)
d2 <- data.table(chrom = chrom_sel, start = chrom_len + 1L, end = chrom_len, midpoint = chrom_len + 1L, class = "matched", delta = 1)
db <- rbind(
  diffTADs_dt,
  data.table(d1, allele = "w"), data.table(d1, allele = "b"),
  data.table(d2, allele = "w"), data.table(d2, allele = "b"),
  fill = T
)
setkey(db, allele, chrom, end)

# make sure that there are no boundaries in deleted/duplicated intervals...
bp_int_deldup <- bp_int_dt[end - start <= max(abs(bp_gr$duplication_size))]
bp_int_deldup_gr <- GRanges(bp_int_deldup$chrom, IRanges(bp_int_deldup$start + 1L, bp_int_deldup$end))
stopifnot(length(findOverlaps(bp_int_deldup_gr, GRanges(diffTADs_dt))) == 0)
# ...and that each balancer block has a boundary
bp_int_block <- bp_int_dt[end - start > max(abs(bp_gr$duplication_size))]
bp_int_block_gr <- GRanges(bp_int_block$chrom, IRanges(bp_int_block$start + 1L, bp_int_block$end))
stopifnot(length(unique(queryHits(findOverlaps(bp_int_block_gr, GRanges(diffTADs_dt))))) == length(bp_int_block_gr))

# identify TADs disrupted by breakpoints
extract_disrupted_TADs <- function(bp_gr, db, delta_thr = 0, prominence_thr = 0, allele = NA)
{
  this_allele <- allele
  db_subset <- db[class == "matched" | (delta > delta_thr & prominence > prominence_thr), ]
  if (!is.na(this_allele))
    db_subset <- db_subset[allele == this_allele, ]

  sel <- with(db_subset, which(end == 0 | delta != 1)) # db_subset$delta == 1 for chromosome start/end positions
  n_perturbed_boundaries <- length(sel)

  stopifnot(db_subset$chrom[sel] == db_subset$chrom[sel + 1])
  TAD_gr <- GRanges(db_subset$chrom[sel], IRanges(db_subset$end[sel] + 1L, db_subset$end[sel + 1]))
  TAD_gr$class_start <- db_subset$class[sel]
  TAD_gr$class_end <- db_subset$class[sel + 1]
  stopifnot(length(TAD_gr) == length(reduce(TAD_gr, min.gapwidth = 0)))

  ov <- findOverlaps(bp_gr, TAD_gr + 1L)
  # stopifnot(queryHits(ov) == seq_along(bp_gr))
  bp_to_TAD <- data.table(bp_id = bp_gr[queryHits(ov)]$id, chrom = as.character(seqnames(bp_gr[queryHits(ov)])),
    bp_start = start(bp_gr[queryHits(ov)]), bp_end = end(bp_gr[queryHits(ov)]), bp_duplication_size = bp_gr[queryHits(ov)]$duplication_size,
    TAD_start = start(TAD_gr[subjectHits(ov)]), TAD_start_class = TAD_gr[subjectHits(ov)]$class_start,
    TAD_end = end(TAD_gr[subjectHits(ov)]), TAD_end_class = TAD_gr[subjectHits(ov)]$class_end
  )
  bp_to_TAD[, distance_start := bp_start - TAD_start]
  bp_to_TAD[, distance_end := TAD_end - bp_end]
}

# add a row for chr3_Rend needed after liftOver to the balancer assembly
end_gr <- GRanges("chr3R", IRanges(dm6_chrom_sizes[chrom == "chr3R"]$chrom_size + 1L, dm6_chrom_sizes[chrom == "chr3R"]$chrom_size))
end_gr$id <- "chr3_Rend"
end_gr$duplication_size <- 0

# add an inversion (INV2)
INV_dt <- fread(evalq(FN.bed.sv.inv, env.dm6), col.names = c("chrom", "start", "end", "name"))
INV_dt <- INV_dt[grepl("^chr[23]", chrom), ]
INV_dt[, start := start + 1L]
INV_gr <- GRanges(INV_dt)

ov <- findOverlaps(GRanges(INV_dt), GRanges(INV_dt))
ov_dt <- data.table(query = queryHits(ov), query_width = width(GRanges(INV_dt))[queryHits(ov)],
  subject = subjectHits(ov), subject_width = width(GRanges(INV_dt))[subjectHits(ov)],
  overlap_width = width(pintersect(GRanges(INV_dt)[queryHits(ov)], GRanges(INV_dt)[subjectHits(ov)])))
ov_dt <- ov_dt[query != subject & grepl("5to5,", INV_dt$name[query]) & grepl("3to3,", INV_dt$name[subject]) &
  abs(overlap_width - query_width) <= 100 & abs(overlap_width - subject_width) <= 100 & overlap_width > 10e3]
INV_selected_dt <- rbind(
  data.table(INV_dt[ov_dt$query], id = seq_along(nrow(ov_dt))),
  data.table(INV_dt[ov_dt$subject], id = seq_along(nrow(ov_dt)))
)

stopifnot(nrow(INV_selected_dt) == 2)
inv2_gr <- GRanges(INV_selected_dt$chrom, IRanges(c(min(INV_selected_dt$start), min(INV_selected_dt$end)),
  c(max(INV_selected_dt$start), max(INV_selected_dt$end))))
inv2_gr$id <- c("chr3_1a", "chr3_1b")
inv2_gr$duplication_size <- c(INV_selected_dt$start[1] - INV_selected_dt$start[2] + 1L, INV_selected_dt$end[1] - INV_selected_dt$end[2] + 1L)

bp_ext_gr <- c(bp_gr, inv2_gr, end_gr)

# TAD halves in the wild type
bp_to_TAD_w <- extract_disrupted_TADs(bp_ext_gr, db, delta_thr = 0, prominence_thr = 0, allele = "w")
halves_w <- with(bp_to_TAD_w, rbind(
  data.table(id = paste0(bp_id, "_left"), other_id = paste0(bp_id, "_right"), bp_id, chrom, start = TAD_start, end = bp_start - 1L,
    class = TAD_start_class, bp_distance = distance_start, bp_duplication_size),
  data.table(id = paste0(bp_id, "_right"), other_id = paste0(bp_id, "_left"), bp_id, chrom, start = bp_end + 1L, end = TAD_end,
    class = TAD_end_class, bp_distance = distance_end, bp_duplication_size)
))
setkey(halves_w, bp_id, id)

# TAD halves in the balancer, but in the reference assembly
bp_to_TAD_b <- extract_disrupted_TADs(bp_ext_gr, db, delta_thr = 0, prominence_thr = 0, allele = "b")
halves_b_ref <- with(bp_to_TAD_b, rbind(
  data.table(id = paste0(bp_id, "_left"), bp_id, chrom, start = TAD_start, end = bp_start - 1L,
    class = TAD_start_class, bp_distance = distance_start, bp_duplication_size),
  data.table(id = paste0(bp_id, "_right"), bp_id, chrom, start = bp_end + 1L, end = TAD_end,
    class = TAD_end_class, bp_distance = distance_end, bp_duplication_size)
))
setkey(halves_b_ref, bp_id, id)

# liftOver to dm6bal3
# remove chr3_Rend_right
dt <- halves_b_ref[!(id == "chr3_Rend_right" & bp_distance == 0), ]
dt[, numeric_id := seq_len(nrow(dt))]
co <- c("chrom", "start", "end", "numeric_id")
setcolorder(dt, c(co, setdiff(names(dt), co)))
dt[, start := start - 1L] # convert to 0-based coordinates
halves_b <- liftOver_to_dm6bal3(dt)
# add a row for chr3_Rend_right in the balancer assembly
halves_b <- rbind(halves_b, data.table(chrom = "chr3R", start = dm6bal3_chrom_sizes[chrom == "chr3R"]$chrom_size,
  end = dm6bal3_chrom_sizes[chrom == "chr3R"]$chrom_size, numeric_id = NA, id = "chr3_Rend_right", bp_id = "chr3_Rend", class = "matched",
  bp_distance = 0, bp_duplication_size = 0
))
halves_b[, start := start + 1L] # convert the coordinates back
stopifnot(dt$numeric_id == seq_len(nrow(dt)))
halves_b[, numeric_id := NULL]

# identify the other_id for TAD halves in the balancer
mt <- NULL
for (i in seq_len(nrow(halves_b)))
  for (j in seq_len(nrow(halves_b)))
  {
    d <- halves_b$end[i] - halves_b$start[j] + 1L
    if (-1000 < d & d <= 0 & i != j &
      !(halves_b$id[i] == "chr3_7_left" & halves_b$id[j] == "chr3_6_right") &
      !grepl("^chr3_1[ab]", halves_b$id[i]) & !grepl("chr3_1[ab]", halves_b$id[j])
    )
      mt <- rbind(mt, data.table(i, j, d))
  }
# manually add matching for INV2
mt <- rbind(mt, data.table(
  i = c(grep("chr3_1a_left", halves_b$id), grep("chr3_1a_right", halves_b$id)),
  j = c(grep("chr3_1b_left", halves_b$id), grep("chr3_1b_right", halves_b$id)),
  d = NA
))
stopifnot(2 * nrow(mt) == nrow(halves_b))
halves_b[, other_id := NA]
for (k in seq_len(nrow(mt)))
{
  halves_b$other_id[mt$i[k]] <- halves_b$id[mt$j[k]]
  halves_b$other_id[mt$j[k]] <- halves_b$id[mt$i[k]]
}

# aggregate together all the information on TAD halves
halves <- with(halves_w, data.table(id, w_chrom = chrom, w_start = start, w_end = end, w_class = class, w_bp_distance = bp_distance, wo_id = other_id))
halves[, wo_start := halves_w$start[match(wo_id, halves_w$id)]]
halves[, wo_end := halves_w$end[match(wo_id, halves_w$id)]]
halves[, wo_class := halves_w$class[match(wo_id, halves_w$id)]]
halves[, wo_bp_distance := halves_w$bp_distance[match(wo_id, halves_w$id)]]
halves[, w_size := pmax(w_end, wo_end) - pmin(w_start, wo_start) + 1L]

halves[, b_chrom := halves_b$chrom[match(id, halves_b$id)]]
halves[, b_start := halves_b$start[match(id, halves_b$id)]]
halves[, b_end := halves_b$end[match(id, halves_b$id)]]
halves[, b_class := halves_b$class[match(id, halves_b$id)]]
halves[, b_bp_distance := halves_b$bp_distance[match(id, halves_b$id)]]
halves[, bo_id := halves_b$other_id[match(id, halves_b$id)]]

halves[, bo_start := halves_b$start[match(bo_id, halves_b$id)]]
halves[, bo_end := halves_b$end[match(bo_id, halves_b$id)]]
halves[, bo_class := halves_b$class[match(bo_id, halves_b$id)]]
halves[, bo_bp_distance := halves_b$bp_distance[match(bo_id, halves_b$id)]]
halves[, b_size := pmax(b_end, bo_end) - pmin(b_start, bo_start) + 1L]

# manually fix b_size for INV2
halves[, bo_bp_distance := ifelse(grepl("^chr3_1[ab]", id), w_bp_distance[match(bo_id, id)], bo_bp_distance)]
stopifnot(halves$b_size[grepl("^chr3_1[ab]_left", halves$id)] == 90869L)
stopifnot(halves$b_size[grepl("^chr3_1[ab]_right", halves$id)] == 146797L)
halves[, b_size := ifelse(grepl("^chr3_1[ab]_left", id), 69078L, b_size)]
halves[, b_size := ifelse(grepl("^chr3_1[ab]_right", id), 130958L, b_size)]

# add reference balancer coordinates for balancer TAD halves
stopifnot(halves$id == halves_b_ref$id)
halves[, b_ref_start := halves_b_ref$start]
halves[, b_ref_end := halves_b_ref$end]

# take only the TADs where breakpoint is >25 kb or >12.5 kb from a TAD boundary
message("Distances to the closest TAD boundary:")
print(sort(with(halves, pmin(ifelse(w_bp_distance == 0, Inf, w_bp_distance), ifelse(wo_bp_distance == 0, Inf, wo_bp_distance)))))
if (suffix == "_thresh25kb")
  halves <- halves[(w_bp_distance == 0 | w_bp_distance > 25e3) & (wo_bp_distance == 0 | wo_bp_distance > 25e3), ]
else if (suffix == "_thresh12.5kb")
  halves <- halves[(w_bp_distance == 0 | w_bp_distance > 12.5e3) & (wo_bp_distance == 0 | wo_bp_distance > 12.5e3), ]
write.table(halves, file = paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_halves.tab"), quote = F, row.names = F)


halves_gr <- with(halves, GRanges(w_chrom, IRanges(w_start, w_end)))

# selection of columns for printing
hs <- halves[, c("id", "w_chrom", "w_start", "w_end", "w_bp_distance", "wo_id", "wo_start", "wo_end", "wo_bp_distance",
  "b_chrom", "b_start", "b_end", "b_bp_distance", "bo_id", "bo_start", "bo_end", "bo_bp_distance")]
# message("Considered TAD halves:")
# print(hs)
message("\nThe breakpoints split the ", nrow(halves) / 2, " disrupted TADs into ", nrow(halves), " halves, which are rearranged differently in the balancer.")
message("In ", sum(with(halves, b_size - w_size > 50e3)), " cases, the halves arrive in a TAD at least 50 kb bigger, and in ", sum(with(halves, b_size - w_size < -50e3)), " cases — at least 50 kb smaller in the balancer.")

# aggregate with DE genes, subset == "non-breakpoint"
genes_subset <- GRanges(genes[!grepl("breakpoint", as.character(trivial_class)), ])
start(genes_subset) <- genes_subset$tss
end(genes_subset) <- genes_subset$tss

# quick statistics: are DE genes enriched in disrupted TADs?
genes_disrupted <- subsetByOverlaps(genes_subset, halves_gr)
# print(table(genes_disrupted$signf))
genes_not_disrupted <- genes_subset[!genes_subset$gene_id %in% genes_disrupted$gene_id]
stopifnot(nrow(genes_disrupted) + nrow(genes_not_disrupted) == nrow(genes_subset))

t <- rbind(table(genes_disrupted[genes_disrupted$signf != "n"]$signf == "s"),
  table(genes_not_disrupted[genes_not_disrupted$signf != "n"]$signf == "s"))
message("Genes in disrupted TADs are preferentially DE (",
  sum(genes_disrupted$signf == "s"), "/", sum(genes_disrupted$signf != "n"), " = ",
  sum(genes_disrupted$signf == "s") / sum(genes_disrupted$signf != "n"),
  ") comparing to genes in non-disrupted TADs (",
  sum(genes_not_disrupted$signf == "s"), "/", sum(genes_not_disrupted$signf != "n"), " = ",
  sum(genes_not_disrupted$signf == "s") / sum(genes_not_disrupted$signf != "n"),
  "), p-value: ", fisher.test(t)$p.value
)

t <- rbind(table(genes_disrupted[genes_disrupted$signf != "n"]$signf == "s"
    & abs(genes_disrupted[genes_disrupted$signf != "n"]$log2FoldChange) > log2(1.5)),
  table(genes_not_disrupted[genes_not_disrupted$signf != "n"]$signf == "s"
    & abs(genes_not_disrupted[genes_not_disrupted$signf != "n"]$log2FoldChange) > log2(1.5)))
message("Genes in disrupted TADs are preferentially strongly DE (",
  sum(genes_disrupted$signf == "s" & abs(genes_disrupted$log2FoldChange) > log2(1.5)), "/", sum(genes_disrupted$signf != "n"), " = ", 
  sum(genes_disrupted$signf == "s" & abs(genes_disrupted$log2FoldChange) > log2(1.5)) / sum(genes_disrupted$signf != "n"),
  ") comparing to genes in non-disrupted TADs (",
  sum(genes_not_disrupted$signf == "s" & abs(genes_not_disrupted$log2FoldChange) > log2(1.5)), "/", sum(genes_not_disrupted$signf != "n"), " = ", 
  sum(genes_not_disrupted$signf == "s" & abs(genes_not_disrupted$log2FoldChange) > log2(1.5)) / sum(genes_not_disrupted$signf != "n"),
  "), p-value: ", fisher.test(t)$p.value
)

# prepare data.table with all genes from the disrupted TADs
ov <- findOverlaps(genes_subset, halves_gr) # here, one gene overlaps two TADs!
halves_to_genes <- as.data.table(genes_subset[queryHits(ov)])
halves_to_genes[, id := halves$id[subjectHits(ov)]]
stopifnot(halves_to_genes[, list(v = anyDuplicated(gene_id)), by = "id"]$v == 0)
write.table(halves_to_genes, file = paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_genes.tab"), quote = F, row.names = F)
write.table(halves_to_genes[signf == "s"]$gene_id,
  file = paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_DE_genes.txt"), quote = F, row.names = F, col.names = F)
write.table(halves_to_genes[signf != "n"]$gene_id,
  file = paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_testable_genes.txt"), quote = F, row.names = F, col.names = F)

halves_to_genes_sum <- halves_to_genes[,
  list(genes = .N, genes_DE = sum(signf == "s"), genes_DE1.5 = sum(signf == "s" & abs(log2FoldChange) > log2(1.5))), by = "id"]


v <- sort(with(halves, b_bp_distance - w_bp_distance))
message("Out of ", nrow(halves) / 2, " disrupted TADs, ", sum(v != 0), " have an allele-specific boundary on one side: ", sum(v < 0), " balancer-specific, and ", sum(v > 0), " wild-type-specific.")
message("The boundary changes on one side by ", min(v) / 1e3, " kb to ", max(v) / 1e3, " kb:")
print(v[v != 0])

# # are the genes in the wild-type-specific part of the TAD affecting it too much?
# print(halves[w_start < b_ref_start | b_ref_end < w_end])
# print(sort(with(halves[w_start < b_ref_start | b_ref_end < w_end], b_bp_distance - w_bp_distance)))
# supress "Each of the 2 combined objects has sequence levels not in the other"
halves_diff_b_gr <- suppressWarnings(c(
  with(halves[w_start < b_ref_start], GRanges(w_chrom, IRanges(w_start, b_ref_start - 1L))),
  with(halves[b_ref_end < w_end], GRanges(w_chrom, IRanges(b_ref_end + 1L, w_end)))
))
genes_disrupted_diff_b <- subsetByOverlaps(genes_subset, halves_diff_b_gr)
cat("Genes in disrupted TADs affected by balancer-specific boundaries:")
print(table(genes_disrupted_diff_b$signf))

# # are the genes in the balancer-specific part of the TAD affecting it too much?
# print(halves[b_ref_start < w_start | w_end < b_ref_end])
# print(sort(with(halves[b_ref_start < w_start | w_end < b_ref_end], b_bp_distance - w_bp_distance)))
# supress "Each of the 2 combined objects has sequence levels not in the other"
halves_diff_w_gr <- suppressWarnings(c(
  with(halves[b_ref_start < w_start], GRanges(w_chrom, IRanges(b_ref_start, w_start - 1L))),
  with(halves[w_end < b_ref_end], GRanges(w_chrom, IRanges(w_end + 1L, b_ref_end)))
))
genes_disrupted_diff_w <- subsetByOverlaps(genes_subset, halves_diff_w_gr)
cat("Genes in disrupted TADs affected by wild-type-specific boundaries:")
print(table(genes_disrupted_diff_w$signf))


hgs <- merge(halves, halves_to_genes_sum, all.x = T, by = "id")
hgs$genes[is.na(hgs$genes)] <- 0L
hgs$genes_DE[is.na(hgs$genes_DE)] <- 0L
hgs$genes_DE1.5[is.na(hgs$genes_DE1.5)] <- 0L

lim <- c(0, with(hgs, max(w_size, b_size))) / 1e3
message("TAD size in wild type and balancer is not correlated: r = ",
  with(hgs, cor(w_size, b_size)), ", n = ", nrow(hgs))

hg <- merge(halves_to_genes, halves, all.x = T, by = "id")
hg[, bp_pos := NA]
hg$bp_pos[grepl("_left$", hg$id)] <- start(bp_ext_gr)[match(sub("_left$", "", hg$id[grepl("_left$", hg$id)]), bp_ext_gr$id)]
hg$bp_pos[grepl("_right$", hg$id)] <- end(bp_ext_gr)[match(sub("_right$", "", hg$id[grepl("_right$", hg$id)]), bp_ext_gr$id)]
hg[, bp_dist := abs(bp_pos - pos)]
message("Distributions of DE and non-DE genes around breakpoints differ, p-value: ",
  with(hg, ks.test(bp_dist[signf == "s"], bp_dist[signf == "i"]))$p.value)

# what are the distances between breakpoint and DE gene within disrupted TADs?
d <- median(hg[signf == "s", ]$bp_dist)
dt <- unique(hg[bp_dist <= d][, c("gene_id", "signf", "signfs", "log2FoldChange")])
cat("Median distance between breakpoint and DE gene within disrupted TADs: ", d, ", DE genes: ",
  sum(dt$signf == "s"), "/", sum(dt$signf != "n"), " (", sum(dt$signf == "s") / sum(dt$signf != "n"), ")")
if ("signfs" %in% names(genes))
  print(table(dt$signfs))

dt <- unique(hg[bp_dist <= 25e3][, c("gene_id", "signf", "signfs", "log2FoldChange")])
cat("As above, withing 25 kb: ", d, ", DE genes: ",
  sum(dt$signf == "s"), "/", sum(dt$signf != "n"), " (", sum(dt$signf == "s") / sum(dt$signf != "n"), ")")
if ("signfs" %in% names(genes))
  print(table(dt$signfs))

message("FPKM distributions for DE and non-DE genes in disrupted TADs do not differ, p-value: ",
  with(hg[signf != "n", ], wilcox.test(FPKM[signf == "s"], FPKM[signf == "i"]))$p.value)


# what about the genes at similar distances, but outside the disrupted TAD?
# consider extended TADs and then subtract
rays_gr <- with(halves, GRanges(w_chrom, IRanges(ifelse(grepl("left", id), 0L, w_start), ifelse(grepl("left", id), w_end, 99e6))))
ov <- findOverlaps(genes_subset, rays_gr)
rays_to_genes <- as.data.table(genes_subset[queryHits(ov)])
rays_to_genes[, id := halves$id[subjectHits(ov)]]
stopifnot(rays_to_genes[, list(v = anyDuplicated(gene_id)), by = "id"]$v == 0)

rg <- merge(rays_to_genes, halves, all.x = T, by = "id")
rg[, bp_pos := NA]
rg$bp_pos[grepl("_left$", rg$id)] <- start(bp_ext_gr)[match(sub("_left$", "", rg$id[grepl("_left$", rg$id)]), bp_ext_gr$id)]
rg$bp_pos[grepl("_right$", rg$id)] <- end(bp_ext_gr)[match(sub("_right$", "", rg$id[grepl("_right$", rg$id)]), bp_ext_gr$id)]
rg[, bp_dist := abs(bp_pos - pos)]

# here only the genes outside the disrupted TAD
rhg <- rg[!paste(id, gene_id) %in% paste(hg$id, hg$gene_id), ]
stopifnot(nrow(rg) == nrow(hg) + nrow(rhg))
message("Distributions of DE and non-DE genes around breakpoints outside disrupted TADs do not differ, p-value: ",
  with(rhg[bp_dist <= max(hg$bp_dist)], ks.test(bp_dist[signf == "s"], bp_dist[signf == "i"]))$p.value)

dt <- unique(rhg[bp_dist <= max(hg$bp_dist), c("gene_id", "signf", "log2FoldChange")])
message("Genes outside disrupted TADs are not preferentially DE (",
  sum(dt$signf == "s"), "/", sum(dt$signf != "n"), " = ",
  sum(dt$signf == "s") / sum(dt$signf != "n"), ")"
)
message("Genes outside disrupted TADs are not preferentially strongly DE (",
  sum(dt$signf == "s" & abs(dt$log2FoldChange) > log2(1.5)), "/", sum(dt$signf != "n"), " = ", 
  sum(dt$signf == "s" & abs(dt$log2FoldChange) > log2(1.5)) / sum(dt$signf != "n"), ")"
)


pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints.pdf"), width = 2.5, height = 2.5)
p <- ggplot(hgs, aes(w_size / 1e3, b_size / 1e3)) +
  geom_abline(slope = 1, lty = 2) +
  geom_point(alpha = 0.6) +
  coord_fixed(xlim = lim, ylim = lim) +
  xlab("TAD size in wild type (kb)") +
  ylab("TAD size in balancer (kb)") +
  background_grid(major = "xy", minor = "xy") +
  NULL
print(p)
dev.off()

pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_genes.pdf"), width = 4, height = 3)

p <- ggplot(hg[signf != "n", ], aes((b_size - w_size) / 1e3, log2FoldChange, color = signf)) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(alpha = 0.6) +
  xlab("Change in TAD size (kb)") +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  scale_y_continuous(limits = c(-1, 1) * max(abs(hg$log2FoldChange))) +
  scale_color_manual(NULL, values = gene.colors, labels = gene.labels) +
  theme(legend.justification = c(0.5, 1), legend.position = "top", legend.direction = "horizontal") +
  background_grid(major = "xy", minor = "xy") +
  NULL
print(p)

p <- ggplot(hg[signf != "n", ], aes(w_bp_distance / 1e3, log2FoldChange, color = signf)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(alpha = 0.6) +
  xlab("Size of a TAD half with the gene (kb)") +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  scale_y_continuous(limits = c(-1, 1) * max(abs(hg$log2FoldChange))) +
  scale_color_manual(NULL, values = gene.colors, labels = gene.labels) +
  theme(legend.justification = c(0.5, 1), legend.position = "top", legend.direction = "horizontal") +
  background_grid(major = "xy", minor = "xy") +
  NULL
print(p)

p <- ggplot(hg[signf != "n", ], aes(wo_bp_distance / 1e3, log2FoldChange, color = signf)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(alpha = 0.6) +
  xlab("Size of the other TAD half in the wild type (kb)") +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  scale_y_continuous(limits = c(-1, 1) * max(abs(hg$log2FoldChange))) +
  scale_color_manual(NULL, values = gene.colors, labels = gene.labels) +
  theme(legend.justification = c(0.5, 1), legend.position = "top", legend.direction = "horizontal") +
  background_grid(major = "xy", minor = "xy") +
  NULL
print(p)

p <- ggplot(hg[signf != "n", ], aes(bo_bp_distance / 1e3, log2FoldChange, color = signf)) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(alpha = 0.6) +
  xlab("Size of the other TAD half in the balancer (kb)") +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  scale_y_continuous(limits = c(-1, 1) * max(abs(hg$log2FoldChange))) +
  scale_color_manual(NULL, values = gene.colors, labels = gene.labels) +
  theme(legend.justification = c(0.5, 1), legend.position = "top", legend.direction = "horizontal") +
  background_grid(major = "xy", minor = "xy") +
  NULL
print(p)

dev.off()


pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_distance.pdf"), width = 3.5, height = 2.75)
p <- ggplot(hg[signf != "n", ], aes(bp_dist / 1e3, fill = signf)) +
  facet_grid(signf ~ ., labeller = as_labeller(gene.labels), scales = "free") +
  geom_histogram(binwidth = 5, boundary = 0) +
  scale_fill_manual(guide = F, values = gene.colors) +
  xlab("Distance from the breakpoint (kb)") +
  ylab("Gene count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  background_grid(major = "xy", minor = "x") +
  NULL
print(p)

p <- ggplot(hg[signf != "n", ], aes(x = 0, y = bp_dist / 1e3, fill = signf)) +
  facet_grid(signf ~ ., labeller = as_labeller(gene.labels)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(guide = F, values = gene.colors) +
  xlab("Gene density") +
  ylab("Distance from the breakpoint (kb)") +
  # scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_flip() +
  NULL
print(p)
dev.off()


pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_FPKM.pdf"), width = 5, height = 2)
p <- ggplot(hg[signf != "n", ], aes(log10(FPKM), fill = signf)) +
  facet_grid(signf ~ ., labeller = as_labeller(c(s = paste0("DE\ngenes\n(", sum(hg$signf == "s"), ")"), i = paste0("non\u00adDE\ngenes\n(", sum(hg$signf == "i"), ")"))), scales = "free") +
  geom_histogram(binwidth = 0.2, boundary = 0, alpha = 0.6) +
  scale_fill_manual(guide = F, values = ggbio.gene.colors) +
  xlab(expression(paste(log[10] * " FPKM"))) +
  ylab("Count") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  background_grid(major = "xy", minor = "x") +
  NULL
print(p)
dev.off()


pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_distance_outside_disruptedTAD.pdf"), width = 3.5, height = 2.75)
p <- ggplot(rhg[bp_dist <= max(hg$bp_dist) * 1.1], aes(bp_dist / 1e3, fill = signf)) +
  facet_grid(signf ~ ., labeller = as_labeller(gene.labels), scales = "free") +
  geom_histogram(binwidth = 5, boundary = 0) +
  scale_fill_manual(guide = F, values = gene.colors) +
  xlab("Distance from the breakpoint (kb)") +
  ylab("Gene count") +
  coord_cartesian(xlim = c(0, max(hg$bp_dist) / 1e3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  background_grid(major = "xy", minor = "x") +
  NULL
print(p)

p <- ggplot(rhg[bp_dist <= max(hg$bp_dist) * 1.1], aes(x = 0, y = bp_dist / 1e3, fill = signf)) +
  facet_grid(signf ~ ., labeller = as_labeller(gene.labels)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(guide = F, values = gene.colors) +
  xlab("Gene density") +
  ylab("Distance from the breakpoint (kb)") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_flip(ylim = c(0, max(hg$bp_dist) / 1e3)) +
  NULL
print(p)
dev.off()


pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_misc.pdf"), width = 8, height = 6)

p <- ggplot(hgs, aes(w_size / 1e3, b_size / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  coord_fixed(xlim = lim, ylim = lim) +
  NULL
print(p)

# p <- ggplot(hgs, aes(w_size / 1e3, b_size / 1e3, color = paste(w_bp_distance <= 25e3 | wo_bp_distance <= 25e3, b_bp_distance <= 25e3 | bo_bp_distance <= 25e3))) +
#   geom_point() +
#   coord_fixed(xlim = lim, ylim = lim) +
#   scale_color_discrete("close to bp in wt, close to wt in bal") +
#   NULL
# print(p)

# p <- ggplot(hgs, aes(w_bp_distance / 1e3, (b_size - w_size) / 1e3, color = paste(w_bp_distance <= 25e3 | wo_bp_distance <= 25e3, b_bp_distance <= 25e3 | bo_bp_distance <= 25e3))) +
#   geom_point() +
#   scale_color_discrete("close to bp in wt, close to wt in bal") +
#   NULL
# print(p)

# p <- ggplot(hgs, aes(b_bp_distance / 1e3, (b_size - w_size) / 1e3, color = paste(w_bp_distance <= 25e3 | wo_bp_distance <= 25e3, b_bp_distance <= 25e3 | bo_bp_distance <= 25e3))) +
#   geom_point() +
#   scale_color_discrete("close to bp in wt, close to wt in bal") +
#   NULL
# print(p)

p <- ggplot(hgs, aes(w_size / 1e3, b_size / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  coord_fixed(xlim = lim, ylim = lim) +
  NULL
print(p)

p <- ggplot(hgs, aes(w_bp_distance / 1e3, w_size / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  NULL
print(p)

p <- ggplot(hgs, aes(w_bp_distance / 1e3, (b_size - w_size) / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  NULL
print(p)

p <- ggplot(hgs, aes(pmin(w_bp_distance, wo_bp_distance) / 1e3, (b_size - w_size) / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  NULL
print(p)

p <- ggplot(hgs, aes(b_bp_distance / 1e3, (b_size - w_size) / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  NULL
print(p)

p <- ggplot(hgs, aes(pmin(b_bp_distance, bo_bp_distance) / 1e3, (b_size - w_size) / 1e3, color = paste(genes_DE > 0, genes_DE1.5 > 0))) +
  geom_point() +
  NULL
print(p)

p <- ggplot(hg, aes(w_bp_distance / 1e3, log2FoldChange, color = signf)) +
  geom_point() +
  NULL
print(p)

p <- ggplot(hg, aes(wo_bp_distance / 1e3, log2FoldChange, color = signf)) +
  geom_point() +
  NULL
print(p)

dev.off()

#
#  Manhattan plots around TAD breakpoints
#

dth <- copy(halves)
dth[, id_trimmed := factor(sub("_left|_right", "", id), levels_id_trimmed)]
dth[, allele := "i"]
# dth <- rbind(dth, with(dth[w_start < b_ref_start], data.table(w_start = w_start, w_end = b_ref_start - 1L, id_trimmed = id_trimmed, bp_mean = bp_mean, allele = "w")), fill = T)
# dth <- rbind(dth, with(dth[b_ref_end < w_end], data.table(w_start = b_ref_end + 1L, w_end = w_end, id_trimmed = id_trimmed, bp_mean = bp_mean, allele = "w")), fill = T)
# dth <- rbind(dth, with(dth[b_ref_start < w_start], data.table(w_start = b_ref_start, w_end = w_start - 1L, id_trimmed = id_trimmed, bp_mean = bp_mean, allele = "b")), fill = T)
# dth <- rbind(dth, with(dth[w_end < b_ref_end], data.table(w_start = w_end + 1L, w_end = b_ref_end, id_trimmed = id_trimmed, bp_mean = bp_mean, allele = "b")), fill = T)

dt_bp_mean <- rbind(bp_mean[id %in% dth$id_trimmed, ],
  data.table(id = "chr3_Rend", breakpoint = dm6_chrom_sizes$chrom_size[match("chr3R", dm6_chrom_sizes$chrom)]),
  data.table(id = "chr3_1a", breakpoint = round(mean(INV_selected_dt$start))),
  data.table(id = "chr3_1b", breakpoint = round(mean(INV_selected_dt$end))),
  fill = T)
dth[, bp_mean := dt_bp_mean$breakpoint[match(id_trimmed, dt_bp_mean$id)]]

dt_bp_unique <- rbind(bp_unique[id %in% dth$id_trimmed, ],
  data.table(id = "chr3_1a", breakpoint = INV_selected_dt$start),
  data.table(id = "chr3_1b", breakpoint = INV_selected_dt$end),
  fill = T)
dt_bp_unique[, id_trimmed := factor(id, levels_id_trimmed)]
dt_bp_unique[, bp_mean := dt_bp_mean$breakpoint[match(id_trimmed, dt_bp_mean$id)]]

ymax <- 7.5
dt <- copy(hg)
dt[, id_trimmed := factor(sub("_left|_right", "", id), levels_id_trimmed)]
dt[, bp_mean := dt_bp_mean$breakpoint[match(id_trimmed, dt_bp_mean$id)]]
dt[, y := pmin(-log10(padj), ymax) * sign(log2FoldChange)]


pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_manhattan.pdf"), width = 5.5, height = 8)

p <- ggplot(dt) +
  geom_rect(data = dth, aes(xmin = (w_start - 1L - bp_mean) / 1e3, xmax = (w_end - bp_mean) / 1e3, fill = allele), ymin = -Inf, ymax = Inf, alpha = 0.2) +
  # geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "#E17C05") +
  # geom_hline(yintercept = +log10(0.05), linetype = "dashed", col = "#E17C05") +
  geom_vline(data = dt_bp_unique, aes(xintercept = (breakpoint - bp_mean) / 1e3), color = "#984ea3", alpha = 1, lty = "dashed") +
  # geom_segment(aes(y = 0, yend = y, x = (tss - bp_mean) / 1e3, xend = (tss - bp_mean) / 1e3, col = signf), size = 0.2) +
  # geom_point(aes(x = (tss - bp_mean) / 1e3, y = y, shape = signf, col = signf, alpha = signf), size = 1.1) + 
  geom_segment(aes(y = 0, yend = log2FoldChange, x = (tss - bp_mean) / 1e3, xend = (tss - bp_mean) / 1e3, col = signf), size = 0.2) +
  geom_point(aes(x = (tss - bp_mean) / 1e3, y = log2FoldChange, shape = signf, col = signf, alpha = signf), size = 1.1) + 
  scale_alpha_manual(values = gene.alpha, labels = gene.labels, name = NULL) +
  scale_color_manual(values = ggbio.gene.colors, labels = gene.labels, name = NULL) +
  scale_fill_manual(values = gene.colors, labels = c(i = "disrupted TADs"), name = NULL) +
  scale_shape_manual(values = c(s = 15, i = 16, n = 16), labels = gene.labels, name = NULL) +
  facet_grid(id_trimmed ~ .) +
  # scale_x_continuous(breaks = scales::pretty_breaks(3), labels = format_Mb, expand = c(0,10e3), limits = c(0, NA)) +
  # scale_y_continuous(breaks = c(-10,-5,0,5,10), labels = function(x) {minuslog10(abs(x))}) +
  scale_y_continuous(breaks = -1:1 * ymax, limits = c(-1, 1) * ymax) +
  background_grid(major = "xy", minor = "x") +
  xlab("Genomic position around the breakpoint (kb)") +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
  theme(strip.text.y = element_text(size = 6, angle = 90)) +
  # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  NULL
print(p)

dev.off()

if ("signfs" %in% colnames(genes))
{
  pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_manhattan_separation.pdf"), width = 5.5, height = 8)

  p <- ggplot(dt) +
    geom_rect(data = dth, aes(xmin = (w_start - 1L - bp_mean) / 1e3, xmax = (w_end - bp_mean) / 1e3, fill = allele), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    # geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "#E17C05") +
    # geom_hline(yintercept = +log10(0.05), linetype = "dashed", col = "#E17C05") +
    geom_vline(data = dt_bp_unique, aes(xintercept = (breakpoint - bp_mean) / 1e3), color = "#984ea3", lty = "dashed") +
    # geom_segment(aes(y = 0, yend = y, x = (tss - bp_mean) / 1e3, xend = (tss - bp_mean) / 1e3, col = signfs), size = 0.2) +
    # geom_point(aes(x = (tss - bp_mean) / 1e3, y = y, shape = signfs, col = signfs), size = 1.1) + 
    geom_segment(aes(y = 0, yend = log2FoldChange, x = (tss - bp_mean) / 1e3, xend = (tss - bp_mean) / 1e3, col = signfs), size = 0.2) +
    geom_point(aes(x = (tss - bp_mean) / 1e3, y = log2FoldChange, shape = signfs, col = signfs), size = 1.1) + 
    scale_color_manual(values = ggbio.gene.colors, labels = gene.labels, name = NULL) +
    scale_fill_manual(values = gene.colors, labels = c(i = "disrupted TADs"), name = NULL, guide = F) +
    scale_shape_manual(values = c(s = 15, i = 16, ns = 16, nu = 16), labels = gene.labels, name = NULL) +
    facet_grid(id_trimmed ~ .) +
    # scale_x_continuous(breaks = scales::pretty_breaks(3), labels = format_Mb, expand = c(0,10e3), limits = c(0, NA)) +
    # scale_y_continuous(breaks = c(-10,-5,0,5,10), labels = function(x) {minuslog10(abs(x))}) +
    scale_y_continuous(breaks = -1:1 * ymax, limits = c(-1, 1) * ymax) +
    background_grid(major = "xy", minor = "x") +
    xlab("Genomic position around the breakpoint (kb)") +
    ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
    theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
    theme(strip.text.y = element_text(size = 6, angle = 90)) +
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    NULL
  print(p)

  dev.off()

  pdf(paste0("analysis/balancer/diffTAD", suffix, "_breakpoints_manhattan_separation_jitter.pdf"), width = 5.5, height = 8)

  p <- ggplot(dt) +
    geom_rect(data = dth, aes(xmin = (w_start - 1L - bp_mean) / 1e3, xmax = (w_end - bp_mean) / 1e3, fill = allele), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    # geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "#E17C05") +
    # geom_hline(yintercept = +log10(0.05), linetype = "dashed", col = "#E17C05") +
    geom_vline(data = dt_bp_unique, aes(xintercept = (breakpoint - bp_mean) / 1e3), color = "#984ea3", lty = "dashed") +
    # geom_segment(aes(y = 0, yend = y, x = (tss - bp_mean) / 1e3, xend = (tss - bp_mean) / 1e3, col = signfs), size = 0.2) +
    # geom_point(aes(x = (tss - bp_mean) / 1e3, y = y, shape = signfs, col = signfs), size = 1.1) + 
    geom_point(aes(x = (tss - bp_mean) / 1e3, y = 0, shape = signfs, col = signfs), position = "jitter", size = 1.1) + 
    scale_color_manual(values = ggbio.gene.colors, labels = gene.labels, name = NULL) +
    scale_fill_manual(values = gene.colors, labels = c(i = "disrupted TADs"), name = NULL, guide = F) +
    scale_shape_manual(values = c(s = 15, i = 16, ns = 16, nu = 16), labels = gene.labels, name = NULL) +
    facet_grid(id_trimmed ~ .) +
    # scale_x_continuous(breaks = scales::pretty_breaks(3), labels = format_Mb, expand = c(0,10e3), limits = c(0, NA)) +
    # scale_y_continuous(breaks = c(-10,-5,0,5,10), labels = function(x) {minuslog10(abs(x))}) +
    scale_y_continuous(breaks = NULL, limits = c(-1, 1)) +
    background_grid(major = "xy", minor = "x") +
    xlab("Genomic position around the breakpoint (kb)") +
    ylab(NULL) +
    theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
    theme(strip.text.y = element_text(size = 6, angle = 90)) +
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    NULL
  print(p)

  dev.off()
}
