options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 1)
  stop("Usage:  Rscript fragment_distance_balancer_cap2.R designDir")

designDir <- options[1]

require(data.table)
require(stringr)

genome <- "dm6"
genome_other <- switch(genome, dm6 = "dm6bal3", dm6bal3 = "dm6", stop("unknown genome assembly"))
diagonal_dist <- 500e3

#
# read bait map
#

ftmp <- tempfile()
system(paste0("Rscript ../33_Hi-C/src/R/balancer_split_and_liftOver_bed.R ", designDir, "/", genome, "_DpnII.baitmap ", genome, " ", genome_other, " ", designDir, "/", genome_other, "_DpnII.baitmap.liftedOver_from_", genome_other, " ", ftmp))
cat("\nShould be only deletions:\n")
system(paste("cat", ftmp))
unlink(ftmp)
cat("\n")

baitmap <- fread(paste0(designDir, "/", genome, "_DpnII.baitmap"), sep = "\t", header = F)
setnames(baitmap, paste0("V", 1:5), c("chrom", "start", "end", "id", "name"))

baitmap_other <- fread(paste0(designDir, "/", genome_other, "_DpnII.baitmap.liftedOver_from_", genome_other), sep = "\t", header = F)
baitmap_other <- baitmap_other[, 1:5, with = F] # keep only 5 first columns
setnames(baitmap_other, paste0("V", 1:5), c("chrom", "start", "end", "id_interval", "name"))
baitmap_other[, id := as.integer(str_split_fixed(id_interval, "_", 2)[, 1])]
# baitmap_other$chrom <- factor(baitmap_other$chrom, levels(baitmap$chrom))

#
# read restriction site map
#

ftmp <- tempfile()
system(paste0("Rscript ../33_Hi-C/src/R/balancer_split_and_liftOver_bed.R ", designDir, "/", genome, "_DpnII.rmap ", genome, " ", genome_other, " ", designDir, "/", genome_other, "_DpnII.rmap.liftedOver_from_", genome_other, " ", ftmp))
cat("\nShould be only deletions:\n")
system(paste("cat", ftmp))
unlink(ftmp)
cat("\n")

rmap <- fread(paste0(designDir, "/", genome, "_DpnII.rmap"), sep = "\t", header = F)
setnames(rmap, paste0("V", 1:4), c("chrom", "start", "end", "id"))

rmap_other <- fread(paste0(designDir, "/", genome_other, "_DpnII.rmap.liftedOver_from_", genome_other), sep = "\t", header = F)
rmap_other <- rmap_other[, 1:4, with = F] # keep only 4 first columns
setnames(rmap_other, paste0("V", 1:4), c("chrom", "start", "end", "id_interval"))
rmap_other[, id := as.integer(str_split_fixed(id_interval, "_", 2)[, 1])]
# rmap_other$chrom <- factor(rmap_other$chrom, levels(rmap$chrom))

#
# calculate genomic distance between restriction fragments, in wild-type and balancer genomes
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

rmap_to_distance <- function(baitmap, rmap)
{
  dt.baitmap <- data.table(
    baitID = baitmap$id,
    chrom1 = baitmap$chrom,
    midpoint1 = (baitmap$start + baitmap$end) / 2,
    weight1 = baitmap$end - baitmap$start
  )

  dt.rmap <- data.table(
    otherEndID = rmap$id,
    chrom2 = rmap$chrom,
    midpoint2 = (rmap$start + rmap$end) / 2,
    weight2 = rmap$end - rmap$start
  )

  jdt <- CJ.dt(dt.rmap, dt.baitmap)

  jdt[, distance := ifelse(chrom1 == chrom2, abs(midpoint2 - midpoint1), NA)]
  # remove unused columns
  jdt[, chrom1 := NULL]
  jdt[, chrom2 := NULL]
  jdt[, midpoint1 := NULL]
  jdt[, midpoint2 := NULL]

  # # estimate the number of contacts, make sure the estimate is not smaller than trans-interaction frequency
  # jdt$value[!is.na(jdt$distance)] <- pmax(predict(distance_fit, newdata = jdt), mean_trans)
  # jdt$value[is.na(jdt$distance)] <- mean_trans
  # jdt[, value := value * weight1 * weight2]
  # remove unused columns
  jdt[, weight1 := NULL]
  jdt[, weight2 := NULL]
  # # normalize the estimate by the total number of contacts
  # jdt$value <- jdt$value / count

  # # append "null" data table to make sure that all bin pairs are covered
  # dt_null <- data.table(baitID = seq_len(d))
  # jdt_null <- CJ.dt(dt_null, dt_null)
  # setnames(jdt_null, "i.baitID", "otherEndID")
  # jdt_null <- jdt_null[baitID <= otherEndID, ]
  # jdt_null[, distance := Inf]
  # jdt_null[, value := 0]

  # jdt_combined <- rbind(jdt, jdt_null)
  return(jdt[, list(distance = min(distance), value = 1), by = c("baitID", "otherEndID")])
  # FIXME add distance-based decay as normalization factors
  # FIXME calculate distance as min instead of distance between midpoints?
}

dist_fit <- rmap_to_distance(baitmap, rmap)
dist_fit_other <- rmap_to_distance(baitmap_other, rmap_other)

setnames(dist_fit_other, "distance", "distance_other")
setnames(dist_fit_other, "value", "value_other")
dist_fit <- merge(dist_fit, dist_fit_other, by = c("baitID", "otherEndID"), all = T)

#
# limit the analysis to bin pairs within a certain distance (either in wild-type or in balancer genome)
#

sel <- dist_fit$distance <= diagonal_dist | dist_fit$distance_other <= diagonal_dist | diagonal_dist == Inf
sel[is.na(sel)] <- F
sel[dist_fit$value == 0 | dist_fit$value_other == 0] <- F
dist_fit <- dist_fit[sel, ]

# FIXME consider here also fragments crossing the breakpoint boundaries
dist_fit <- dist_fit[!is.na(distance) & !is.na(distance_other) & distance == distance_other, ]

save(diagonal_dist, dist_fit, file = paste0(designDir, "/", genome, "_", genome_other, "_DpnII_fragment_distance.Rdata"))
