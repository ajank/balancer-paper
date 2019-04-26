require(data.table)
require(GenomicRanges)

bin_size <- 5000L
genome <- "dm6"

chrom_map <- read.table(paste0("data/hic_maps/chrom_map_", genome, "_", bin_size, ".txt"), header = T, sep = "\t")

boundary_gr <- GRanges()
for (i in 1:nrow(chrom_map))
{
  len <- chrom_map$map_end[i] - chrom_map$map_start[i] + 1
  interval_start <- 0:(len - 1) * bin_size + 1L
  interval_end <- pmin(1:len * bin_size, chrom_map$length[i])
  midpoint <- 0:(len - 1) * bin_size + 1L
  boundary_gr <- c(boundary_gr, GRanges(chrom_map$chrom[i], IRanges(midpoint, midpoint), midpoint = midpoint, interval_start = interval_start, interval_end = interval_end))
}

liftOver_to_dm6 <- function(d)
{
  ftmp1 <- tempfile()
  ftmp2 <- tempfile()
  ftmp3 <- tempfile()
  write.table(d, file = ftmp1, sep = "\t", quote = F, row.names = F, col.names = F)
  system(paste("Rscript src/R/balancer_split_and_liftOver_bed.R", ftmp1, "dm6bal3", "dm6", ftmp2, ftmp3))
  cat("\nShould be only deletions:\n")
  system(paste("cat", ftmp3))
  cat("\n")

  dl <- fread(ftmp2, sep = "\t", header = F)
  unlink(c(ftmp1, ftmp2, ftmp3))
  dl[, V4 := stringr::str_split_fixed(V4, "_", 2)[, 1]]
  if (length(names(d)) == 3)
    dl[, V4 := NULL]
  else if (typeof(d[[4]]) == "integer")
    dl[, V4 := as.integer(V4)]
  else if (typeof(d[[4]]) == "double")
    dl[, V4 := as.double(V4)]
  names(dl) <- names(d)[seq_along(names(dl))]
  return(dl)
}

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
  dl[, V4 := stringr::str_split_fixed(V4, "_", 2)[, 1]]
  if (length(names(d)) == 3)
    dl[, V4 := NULL]
  else if (typeof(d[[4]]) == "integer")
    dl[, V4 := as.integer(V4)]
  else if (typeof(d[[4]]) == "double")
    dl[, V4 := as.double(V4)]
  names(dl) <- names(d)[seq_along(names(dl))]
  return(dl)
}

read_IS_raw <- function(dataset, genome = "dm6", target_genome = "dm6", source = "filtered_5000", colnum = NA)
{
  fname <- paste0("data/hicexplorer/h5/", genome, "_", dataset, "_", source, "_corrected_tad_score.bm")
  message("Reading ", fname)
  d <- fread(fname, header = F, sep = "\t", skip = 1)
  names(d)[1:3] <- c("chrom", "start", "end")

  # take the values from a given column -- or take the mean if the column is not specified
  if (is.na(colnum))
    d$value <- apply(d[, 4:ncol(d)], 1, mean)
  else
    d$value <- apply(d[, 4:ncol(d)], 1, "[", colnum)

  d <- d[, c("chrom", "start", "end", "value")]

  # convert to dm6 genome assembly if needed
  if (genome == "dm6bal3" && target_genome == "dm6")
    d <- liftOver_to_dm6(d)

  # convert to 1-based coordinates
  d[, start := start + 1L]

  d_gr <- GRanges(d)
  return(d_gr)
}

read_IS <- function(dataset, suffix, genome = "dm6", source = "filtered_5000", colnum = NA)
{
  d_gr <- read_IS_raw(dataset = dataset, genome = genome, source = source, colnum = colnum)
  ov <- findOverlaps(boundary_gr, d_gr)

  stopifnot(!anyDuplicated(queryHits(ov)))
  elementMetadata(boundary_gr)[, paste0("IS", suffix)] <<- NA
  elementMetadata(boundary_gr)[, paste0("IS", suffix)][queryHits(ov)] <<- d_gr$value[subjectHits(ov)]
}

read_IS_all_columns <- function(dataset, genome = "dm6", source = "filtered_5000", bin_size = 5000L)
{
  fname <- paste0("data/hicexplorer/h5/", genome, "_", dataset, "_", source, "_corrected_tad_score.bm")
  message("Reading ", fname)
  d <- fread(fname, header = F, sep = "\t", skip = 1)
  names(d)[1:3] <- c("chrom", "start", "end")

  # split the regions larger than bin_size
  sel <- d$end - d$start > bin_size
  ds <- d[sel, ]
  d$end[sel] <- d$end[sel] - bin_size
  ds[, start := end - bin_size]
  d <- rbind(d, ds)

  # take only the coordinate where the IS value is calculated, +/- 1 bp
  d[, end := start + 1L]
  d[, start := start - 1L]

  # convert to dm6 genome assembly if needed
  if (genome == "dm6bal3")
    d <- liftOver_to_dm6(d)

  # convert to 1-based coordinates
  d[, start := start + 1L]

  dt <- NULL
  for (i in 4:ncol(d))
      dt <- rbind(dt, data.table(chrom = d$chrom, start = d$start, end = d$end, midpoint = d$start, value = d[[i]], group = i - 3L))
  setkey(dt, chrom, start)
  return(dt)
}

read_PC1 <- function(dataset, suffix, genome = "dm6", source = "filtered_5000", colnum = NA)
{
  fname <- paste0("data/hicexplorer/h5/", genome, "_", dataset, "_", source, "_corrected_PCA1.bedgraph")
  message("Reading ", fname)
  d <- fread(fname, header = F, sep = "\t")
  names(d)[1:4] <- c("chrom", "start", "end", "value")

  # convert to dm6 genome assembly if needed
  if (genome == "dm6bal3")
    d <- liftOver_to_dm6(d)

  # convert to 1-based coordinates
  d[, start := start + 1L]

  d_gr <- GRanges(d)
  ov <- findOverlaps(boundary_gr, d_gr)

  stopifnot(!anyDuplicated(queryHits(ov)))
  elementMetadata(boundary_gr)[, paste0("PC1", suffix)] <<- NA
  elementMetadata(boundary_gr)[, paste0("PC1", suffix)][queryHits(ov)] <<- d_gr$value[subjectHits(ov)]
}

extract_max_value_scaled <- function(is_gr, ind, which = "left")
{
  is_dt <- data.table(chrom = as.character(seqnames(is_gr)), value_scaled = is_gr$value_scaled)
  if (which == "left")
    is_dt[, boundary := cut(seq_along(is_gr), c(0, ind), labels = F)]
  else
    is_dt[, boundary := cut(seq_along(is_gr), c(ind, length(is_gr)), labels = F)]
  is_dt <- is_dt[, list(max_value_scaled = max(value_scaled)), by = c("chrom", "boundary")]

  b_dt <- data.table(chrom = as.character(seqnames(is_gr[ind])), boundary = seq_along(ind))
  b_dt <- merge(b_dt, is_dt, by = c("chrom", "boundary"), all.x = T)
  return(b_dt$max_value_scaled)
}

calculate_delta <- function(is_gr, ind)
{
  chrom <- as.integer(as.factor(as.character(seqnames(is_gr))))
  value_scaled <- is_gr$value_scaled
  v <- rep(NA, length(ind))

  for (i in seq_along(ind))
  {
    ii <- ind[i]
    iv <- ii + c(-10:-1, 1:10)
    iv <- iv[chrom[iv] == chrom[ii]]
    v[i] <- mean(value_scaled[iv], na.rm = T) - value_scaled[ii]
  }

  return(v)
}

annotate_boundaries <- function(d, dataset, genome, domain_source)
{
  # annotate using complete IS profiles
  is_gr <- read_IS_raw(dataset = dataset, genome = genome, target_genome = genome, source = domain_source, colnum = NA)
  d_gr <- GRanges(d$chrom, IRanges(d$end, d$end))

  load("analysis/balancer/differential_IS.Rdata") # lm.fit
  b.diffIS <- b

  if (grepl("VRG", dataset))
  {
    is_gr$value_scaled <- predict(lm.fit, data.table(IS.VRG = is_gr$value))
    d$tad_sep_scaled <- predict(lm.fit, data.table(IS.VRG = d$tad_sep))
  }
  else
  {
    is_gr$value_scaled <- is_gr$value
    d$tad_sep_scaled <- d$tad_sep
  }

  ov <- findOverlaps(d_gr, is_gr)
  stopifnot(queryHits(ov) == seq_len(nrow(d)))
  # print(summary(abs(d$tad_sep - is_gr$value[subjectHits(ov)])))
  # stopifnot(abs(d$tad_sep - is_gr$value[subjectHits(ov)])[!is.na(d$tad_sep)] < 1e-6)
  d[, value_scaled := is_gr$value_scaled[subjectHits(ov)]]
  d[, max_value_left := extract_max_value_scaled(is_gr, subjectHits(ov), which = "left")]
  d[, max_value_right := extract_max_value_scaled(is_gr, subjectHits(ov), which = "right")]
  d[, prominence := pmin(max_value_left, max_value_right) - value_scaled]
  d[, mydelta := calculate_delta(is_gr, subjectHits(ov))]

  if (genome == "dm6bal3")
    d <- liftOver_to_dm6(d)
  setkey(d, chrom, start)

  d_gr <- GRanges(d$chrom, IRanges(d$end, d$end))
  ov <- findOverlaps(d_gr, GRanges(b.diffIS))
  stopifnot(queryHits(ov) == seq_along(ov))

  IS.VRG <- rep(NA, nrow(d))
  IS.VRG[queryHits(ov)] <- b.diffIS$IS.VRG[subjectHits(ov)]
  d[, IS.VRG := IS.VRG]

  IS.BAL <- rep(NA, nrow(d))
  IS.BAL[queryHits(ov)] <- b.diffIS$IS.BAL[subjectHits(ov)]
  d[, IS.BAL := IS.BAL]

  IS.BAL.predicted.from.VRG <- rep(NA, nrow(d))
  IS.BAL.predicted.from.VRG[queryHits(ov)] <- b.diffIS$IS.BAL.predicted.from.VRG[subjectHits(ov)]
  d[, IS.BAL.predicted.from.VRG := IS.BAL.predicted.from.VRG]

  diff.IS <- rep(NA, nrow(d))
  diff.IS[queryHits(ov)] <- b.diffIS$diff.IS[subjectHits(ov)]
  d[, diff.IS := diff.IS]

  d_width <- d[,
    list(
      end = end,
      width_left = end - c(1L, head(end, -1)),
      width_right = c(tail(end, -1), chrom_map$length[match(chrom, chrom_map$chrom)] + 1L) - end
    ),
    by = "chrom"]
  d <- merge(d, d_width, by = c("chrom", "end"), all.x = T)
  d[, width_smaller := pmin(width_left, width_right)]
  d[, width_bigger := pmax(width_left, width_right)]

  # get rid of +/- 1 bp
  d[, start := start + 1L]
  d[, end := end - 1L]
  stopifnot(d$end - d$start == 0)
  # convert to 1-based coordinates
  d[, start := start + 1L]

  d[, chrom_length := chrom_map$length[match(chrom, chrom_map$chrom)]]
  d <- d[d$start <= d$chrom_length, ]
  d[, chrom_length := NULL]

  d[, midpoint := start]
  # d[, allele := NA]

  return(d)
}

read_boundaries <- function(dataset, genome = "dm6", domain_source = "filtered_5000", only23 = F)
{
  fname <- paste0("data/hicexplorer/h5/", genome, "_", dataset, "_", domain_source, "_corrected_domains.bed")
  message("Reading ", fname)
  d <- fread(fname, header = F, sep = "\t")
  names(d)[1:3] <- c("chrom", "start", "end")

  # take only TAD boundaries, +/- 1 bp
  d <- with(d, unique(rbind(
      data.table(chrom, start = start - 1L, end = start + 1L),
      data.table(chrom, start = end - 1L, end = end + 1L)
  )))
  setkey(d, chrom, start)
  if (only23)
    d <- d[grepl("^chr[23]", chrom), ]

  # read TAD boundary annotations separately
  fname <- paste0("data/hicexplorer/h5/", genome, "_", dataset, "_", domain_source, "_corrected_boundaries.gff")
  message("Reading ", fname)
  ba <- rtracklayer::import(fname)
  if (only23)
    ba <- ba[grepl("^chr[23]", as.character(seqnames(ba)))]

  stopifnot(nrow(d) == length(ba))
  stopifnot(d$chrom == as.character(seqnames(ba)))
  stopifnot(start(ba) <= d$start + 1)
  stopifnot(d$end - 1 <= end(ba))
  stopifnot(abs(ba$score - as.numeric(ba$tad_sep)) < 1e-12)

  d$delta <- as.numeric(ba$delta)
  d$pvalue <- as.numeric(ba$pvalue)
  d$tad_sep <- as.numeric(ba$tad_sep)

  d <- annotate_boundaries(d, dataset = dataset, genome = genome, domain_source = domain_source)

  message("Read ", nrow(d), " boundaries")
  return(d)
}
