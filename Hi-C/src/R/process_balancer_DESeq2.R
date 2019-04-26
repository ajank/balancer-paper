args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
{
  message("Usage:  Rscript process_balancer_DESeq2.R <dataset> {dm6|dm6bal3} <bin size> <output rda file>")
  q()
}

dataset <- args[1]
genome <- args[2]
bin_size <- as.integer(args[3])
fout <- args[4]

# dataset <- "HiC_DB_6-8h"
# genome <- "dm6"
# bin_size <- as.integer(50e3)
# fout <- "data/hicexplorer/rda/dm6_HiC_DB_6-8h_combined_BAL_vs_VRG_filtered_50000_corrected.rda"

# dataset <- "HiC_DB_6-8h"
# genome <- "dm6"
# bin_size <- as.integer(5e3)
# fout <- "data/hicexplorer/rda/dm6_HiC_DB_6-8h_combined_BAL_vs_VRG_dist_100kb_filtered_5000_corrected.rda"

# require(BiocParallel)
# register(MulticoreParam(16))
require(DESeq2)
require(stringr)

source("src/R/functions_HiCExplorer.R")

across_breakpoint <- !grepl("_not_across_breakpoint", fout)
diagonal_dist_str <- sub("^_dist_", "", sub("kb$", "e3", regmatches(fout, regexpr("_dist_[^_]*", fout))))
diagonal_dist <- if (length(diagonal_dist_str) > 0) as.double(diagonal_dist_str) else Inf
message("across_breakpoint: ", across_breakpoint, ", diagonal_dist: ", diagonal_dist)

allele_VRG <- "VRG"
allele_BAL <- "BAL"
alleles <- c(allele_VRG, allele_BAL)
replicates <- c("R1", "R2")
genome_VRG <- "dm6"
genome_BAL <- "dm6bal3"

genome_other <- switch(genome, dm6 = "dm6bal3", dm6bal3 = "dm6", stop("unknown genome assembly"))

#
# put together all the input data
# in particular, fill in countData and arr.ind, having the same number of rows
#

prefix <- sub(".*/", "", sub(".rda$", "", fout))
message(prefix)

countData <- NULL
colData <- NULL

for (allele in alleles)
{
  message("\n")
  fin <- paste0("data/hicexplorer/txt/", genome, "_", dataset, "_combined_", allele, "_filtered_", bin_size, "_corrected.txt.gz")
  fm <- fread_map(fin)
  chrom_map <- fm$chrom_map
  bins <- fm$bins
  params <- list(genome = genome, bin_size = bin_size, bins = nrow(bins))
  rm(fm)

  if (is.null(countData))
  {
    d <- max(chrom_map$map_end) # dimension of all contact matrices

    countData <- data.frame(fill = 1:(d * (d + 1) / 2))
    countData$fill <- NULL # tricky way to construct data.frame of appropriate number of rows
    colData <- data.frame(allele = factor(NULL, alleles), replicate = factor())

    # bins_masked <- rep(F, d)
  }

  # message(allele, ": out of ", length(bins$masked), " bins, removing ", sum(bins$masked), " and keeping ", sum(!bins$masked))
  # bins_masked <- bins_masked | bins$masked

  for (replicate in replicates)
  {
    id <- paste0(genome, "_", dataset, "_", replicate, "_", allele, "_filtered_", bin_size)
    colData <- rbind(colData, data.frame(allele = allele, replicate = replicate))

    f_id <- paste0("data/hicexplorer/txt/", id, ".txt.gz")
    t_id <- fread_map(f_id)$map

    # account for the fact that values on the diagonal are NOT double-counted
    # stopifnot(diag(t_id) %% 2L == 0L)
    # diag(t_id) <- diag(t_id) %/% 2L

    countData[, id] <- t_id[!lower.tri(t_id)]
  }
}

# message("\nIn total: out of ", length(bins_masked), " bins, removing ", sum(bins_masked), " and keeping ", sum(!bins_masked))

arr.ind <- which(!lower.tri(t_id), arr.ind = T)
colnames(arr.ind) <- c("bin1", "bin2")
stopifnot(nrow(arr.ind) == nrow(countData))
rm(t_id)

countData <- as.matrix(countData)

#
# calculate bin positions along wild-type and balancer genomes
#

bin_bed <- data.table(
  chrom = rep(chrom_map$chrom, chrom_map$map_end - chrom_map$map_start + 1L),
  chrom_size = rep(chrom_map$length, chrom_map$map_end - chrom_map$map_start + 1L),
  start = as.integer(unlist(lapply(chrom_map$map_end - chrom_map$map_start + 1, function(len) seq(from = 0L, length.out = len))) * bin_size)
)
bin_bed[, end := pmin(bin_bed$start + bin_size, bin_bed$chrom_size)]
bin_bed[, chrom_size := NULL]
bin_bed[, id := seq_len(nrow(bin_bed))]

# take only balancer chromosomes
bin_bed <- bin_bed[grepl("^chr[23]", chrom), ]

ftmp1 <- tempfile()
ftmp2 <- tempfile()
ftmp3 <- tempfile()
write.table(bin_bed, file = ftmp1, sep = "\t", quote = F, row.names = F, col.names = F)
system(paste("Rscript src/R/balancer_split_and_liftOver_bed.R", ftmp1, genome, genome_other, ftmp2, ftmp3))
cat("\nShould be only deletions:\n")
system(paste("cat", ftmp3))
cat("\n")

bin_bed_other <- fread(ftmp2, sep = "\t", header = F)
bin_bed_other <- bin_bed_other[, 1:4, with = F] # keep only 4 first columns
setnames(bin_bed_other, paste0("V", 1:4), c("chrom", "start", "end", "id_interval"))
bin_bed_other[, id := as.integer(stringr::str_split_fixed(id_interval, "_", 2)[, 1])]
unlink(c(ftmp1, ftmp2, ftmp3))

bin_bed_VRG <- if (genome == genome_VRG) bin_bed else bin_bed_other
bin_bed_BAL <- if (genome == genome_BAL) bin_bed else bin_bed_other

#
# calculate genomic distance between bins, in wild-type and balancer genomes
# bin_fit should have the same number of rows as countData and arr.ind
#

# http://stackoverflow.com/questions/25888706/r-data-table-cross-join-not-working/27347397#27347397
# There is no cross join functionality available in data.table out of the box.
CJ.dt = function(X,Y)
{
  stopifnot(is.data.table(X),is.data.table(Y))
  k = NULL
  X = X[, c(k=1, .SD)]
  setkey(X, k)
  Y = Y[, c(k=1, .SD)]
  setkey(Y, NULL)
  X[Y, allow.cartesian=TRUE][, k := NULL][]
}

bin_bed_to_fit <- function(bin_bed, genome, replicate, allele)
{
  dt <- data.table(
    bin1 = bin_bed$id,
    chrom1 = bin_bed$chrom,
    midpoint1 = (bin_bed$start + bin_bed$end) / 2,
    weight1 = (bin_bed$end - bin_bed$start) / bin_size
  )

  # consider all pairs of genomic bins
  jdt <- CJ.dt(dt, dt)
  setnames(jdt, c("i.bin1", "i.chrom1", "i.midpoint1", "i.weight1"), c("bin2", "chrom2", "midpoint2", "weight2"))
  jdt <- jdt[bin1 <= bin2, ]

  # load distance-decay fit
  fname <- paste0("data/distance_decay/decay_", genome, "_", dataset, "_", replicate, "_", allele, "_filtered_", bin_size, ".rda")
  cat(paste0("Loading ", fname, "\n"))
  load(fname)

  # take the distance fit only for chr2[LR], chr3[LR]
  require(locfit, quietly = T)
  # distance_fit <- fit_raw$chr23$fit.glm.poisson
  distance_fit <- fit_raw$chr23$locfit_nolog
  # very nasty trick to make locfit predict anything
  attr(distance_fit$terms, ".Environment") <- globalenv()
  mean_trans <- fit_raw$chr23$mean_trans
  # but normalize by the number of reads genome-wide
  count <- fit_raw$genome$count

  # distance == Inf for trans-chromosomal contacts, distance == NA for bin pairs that are not considered
  jdt[, distance := ifelse(chrom1 == chrom2, abs(midpoint2 - midpoint1), Inf)]
  # expected distance for same-bin interactions
  jdt$distance[jdt$distance == 0] <- bin_size / 3
  # remove unused columns
  jdt[, chrom1 := NULL]
  jdt[, chrom2 := NULL]
  jdt[, midpoint1 := NULL]
  jdt[, midpoint2 := NULL]

  # estimate the number of contacts, make sure the estimate is not smaller than trans-interaction frequency
  sel <- is.finite(jdt$distance)
  jdt$value <- Inf
  jdt$value[sel] <- pmax(predict(distance_fit, newdata = jdt[sel, ]), mean_trans)
  jdt$value[!sel] <- mean_trans
  jdt[, value := value * weight1 * weight2]
  stopifnot(jdt$value > 0)
  # remove unused columns
  jdt[, weight1 := NULL]
  jdt[, weight2 := NULL]
  # normalize the estimate by the total number of contacts
  jdt$value <- jdt$value / count

  # append "null" data table to make sure that all bin pairs are covered
  dt_null <- data.table(bin1 = seq_len(d))
  jdt_null <- CJ.dt(dt_null, dt_null)
  setnames(jdt_null, "i.bin1", "bin2")
  jdt_null <- jdt_null[bin1 <= bin2, ]
  jdt_null[, distance := Inf]
  jdt_null[, value := 0]

  jdt_combined <- rbind(jdt, jdt_null)
  jdt_combined <- jdt_combined[, list(distance = min(distance), value = sum(value)), by = c("bin1", "bin2")]
  # marks bin pairs that are not considered with distance == NA, value == NA
  jdt_combined$distance[jdt_combined$value == 0] <- NA
  jdt_combined$value[jdt_combined$value == 0] <- NA

  # put rows in the same order as in arr.ind
  setkey(jdt_combined, bin2, bin1)
  return(jdt_combined)
}

# note that here we use the distance-decay fit from each replicate

bin_fit <- NULL
normMatrix_raw <- NULL

for (allele in alleles)
  for (replicate in replicates)
  {
    # here use the distance-decay fit from the matching genome assembly
    this_bin_fit <- bin_bed_to_fit(
      if (allele == allele_VRG) bin_bed_VRG else bin_bed_BAL,
      if (allele == allele_VRG) genome_VRG else genome_BAL,
      replicate,
      allele
    )

    stopifnot(nrow(this_bin_fit) == nrow(arr.ind))
    stopifnot(this_bin_fit$bin1 == arr.ind[, "bin1"])
    stopifnot(this_bin_fit$bin2 == arr.ind[, "bin2"])

    if (is.null(bin_fit))
    {
      bin_fit <- this_bin_fit[, c("bin1", "bin2")]
      normMatrix_raw <- data.table(V1 = this_bin_fit$value)
    }
    else
    {
      stopifnot(nrow(this_bin_fit) == nrow(bin_fit))
      stopifnot(this_bin_fit$bin1 == bin_fit$bin1)
      stopifnot(this_bin_fit$bin2 == bin_fit$bin2)
      normMatrix_raw <- cbind(normMatrix_raw, this_bin_fit$value)
    }

    if (allele == allele_VRG)
      bin_fit$distance_VRG <- this_bin_fit$distance
    if (allele == allele_BAL)
      bin_fit$distance_BAL <- this_bin_fit$distance
  }

#
# limit the analysis to bin pairs within a certain distance (either in wild-type or in balancer genome)
#

# distance == Inf for trans-chromosomal contacts, distance == NA for bin pairs that are not considered
sel <- with(bin_fit, !is.na(distance_VRG) & !is.na(distance_BAL) & (distance_VRG <= diagonal_dist | distance_BAL <= diagonal_dist))
if (!across_breakpoint)
  sel <- sel & with(bin_fit, is.finite(distance_VRG) & is.finite(distance_BAL) & distance_VRG == distance_BAL)
sel[is.na(sel)] <- F

countData <- countData[sel, ]
arr.ind <- arr.ind[sel, ]
bin_fit <- bin_fit[sel, ]
normMatrix_raw <- normMatrix_raw[sel, ]

#
# now do the real DESeq2 job
#

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ allele)
# rm(countData)
mcols(dds) <- DataFrame(mcols(dds), arr.ind)
# rm(arr.ind)

normMatrix <- normMatrix_raw / exp(rowMeans(log(normMatrix_raw)))

dds <- estimateSizeFactors(dds, normMatrix = normMatrix)
dds <- DESeq(dds) # , parallel = T
res <- results(dds) # , parallel = T
print(summary(res))

#
# save the outcome in a format suitable for plotting
#

map <- list(genome = list())
ind <- as.matrix(mcols(dds)[, c("bin1", "bin2")])

# log2FoldChange
map$genome$log2FoldChange <- matrix(NA_real_, d, d)
map$genome$log2FoldChange[ind] <- res$log2FoldChange
print(summary(as.vector(map$genome$log2FoldChange)))

# minusLog10Padj
map$genome$minusLog10Padj <- matrix(NA_real_, d, d)
map$genome$minusLog10Padj[ind] <- -log10(res$padj)
print(summary(as.vector(map$genome$minusLog10Padj)))

# # distance_VRG
# map$genome$distance_VRG <- matrix(NA_integer_, d, d)
# map$genome$distance_VRG[arr.ind] <- bin_fit$distance_VRG
# print(summary(as.vector(map$genome$distance_VRG)))

# # distance_BAL
# map$genome$distance_BAL <- matrix(NA_integer_, d, d)
# map$genome$distance_BAL[arr.ind] <- bin_fit$distance_BAL
# print(summary(as.vector(map$genome$distance_BAL)))

# # value
# map$genome$value <- matrix(NA_integer_, d, d)
# map$genome$value[arr.ind] <- bin_fit$value
# print(summary(as.vector(map$genome$value)))

# # value_other
# map$genome$value_other <- matrix(NA_integer_, d, d)
# map$genome$value_other[arr.ind] <- bin_fit$value_other
# print(summary(as.vector(map$genome$value_other)))

for (i in seq_len(nrow(chrom_map)))
{
  chrom <- as.character(chrom_map$chrom[i])
  int <- chrom_map$map_start[i]:chrom_map$map_end[i]
  map[[chrom]] <- list(log2FoldChange = map$genome$log2FoldChange[int, int],
    minusLog10Padj = map$genome$minusLog10Padj[int, int])
}

message("Saving ", fout)
save(params, chrom_map, bins, map, file = fout)
message("saved\n")

#
# diagnostics plot to assess the goodness-of-fit
#

dt <- NULL
countDataSums <- apply(countData, 2, sum)
vdist <- 1:200 * 100e3
vdist <- vdist[vdist <= diagonal_dist]
for (dist in vdist)
{
  sel <- dist - bin_size / 2 <= bin_fit$distance_VRG & bin_fit$distance_VRG < dist + bin_size / 2
  sel[is.na(sel)] <- F
  if (sum(sel) > 0)
  {
    expected <- mean(as.matrix(normMatrix_raw[sel, 1:2]))
    v <- apply(countData[sel, , drop = F], 2, mean) / countDataSums
    dt <- rbind(dt, data.frame(distance = dist, expected = expected, observed = mean(v[1:2]), source = "VRG"))
  }

  sel <- dist - bin_size / 2 <= bin_fit$distance_BAL & bin_fit$distance_BAL < dist + bin_size / 2
  sel[is.na(sel)] <- F
  if (sum(sel) > 0)
  {
    expected <- mean(as.matrix(normMatrix_raw[sel, 3:4]))
    v <- apply(countData[sel, , drop = F], 2, mean) / countDataSums
    dt <- rbind(dt, data.frame(distance = dist, expected = expected, observed = mean(v[3:4]), source = "BAL"))
  }
}

if (length(vdist) > 0)
{
  require(ggplot2)
  pdf(paste0("analysis/DESeq2/bin_fit_", prefix, ".pdf"))
  print(ggplot(dt, aes(distance, observed, color = source)) + geom_line() + scale_x_log10() + scale_y_log10() + theme(legend.position="top"))
  print(ggplot(dt, aes(distance, expected, color = source)) + geom_line() + scale_x_log10() + scale_y_log10() + theme(legend.position="top"))
  print(ggplot(dt, aes(expected, observed, color = source)) + geom_line() + scale_x_log10() + scale_y_log10() + theme(legend.position="top"))
  dev.off()
}

#
# save the outcome in a format for debugging
#

save(chrom_map, bin_fit, dds, res, file = paste0("data/DESeq2/", prefix, ".rda"))

#
# make a plot
#

pdf(paste0("analysis/DESeq2/", prefix, ".pdf"))
plotMA(res, main = paste("DESeq2", prefix), ylim = c(-2, 2))
dev.off()
