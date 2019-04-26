args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
{
  message("Usage:  Rscript balancer_split_and_liftOver_bed.R <input BED file> <input assembly, e.g. dm6> <target assembly, e.g. dm6bal3> <output BED file> <unmapped BED file>")
  q()
}

fin <- args[1]
ain <- args[2]
aout <- args[3]

ftemp <- tempfile()
fout_temp <- tempfile()

fout <- args[4]
funmapped <- args[5]

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(GenomicRanges))

# read input BED file
bed_dt <- fread(fin, header = F, sep = "\t")
bed_gr <- GRanges(bed_dt$V1, IRanges(bed_dt$V2 + 1L, bed_dt$V3))

# read balancer breakpoints
balancer.dir <- "/g/furlong/project/39_Balancer_Hi-C"
if (ain == "dm6")
{
  br <- fread(paste0(balancer.dir, "/analysis/breakpoints_dm6_v", sub("^dm6bal", "", aout), ".tab"), header = T)
} else if (aout == "dm6")
{
  br <- fread(paste0(balancer.dir, "/analysis/breakpoints_", ain, ".tab"), header = T)
} else stop("neither input or output assembly is dm6")

# add non-balancer chromosomes
dm6_chrom_sizes <- fread(paste0(balancer.dir, "/liftOver/dm6.chrom.sizes"), header = F, col.names = c("chrom", "chrom_size"))
nonbal_dt <- dm6_chrom_sizes[!(chrom %in% br$chrom), list(chrom = chrom, start = 0L, end = chrom_size)]

# convert breakpoints to intervals
br_unique <- unique(br[, c("id", "chrom", "breakpoint"), with = F])
int_dt <- br_unique[, list(start = head(breakpoint, -1), end = tail(breakpoint, -1)), by = "chrom"]
int_dt <- rbind(int_dt, nonbal_dt)
int_gr <- GRanges(int_dt$chrom, IRanges(int_dt$start + 1L, int_dt$end))

# for each feature, determine which intervals it overlaps
ov <- findOverlaps(bed_gr, int_gr)
missing <- setdiff(seq_along(bed_gr), queryHits(ov))
if (length(missing) > 0)
  stop("missing ", length(missing), " rows among overlap hits, are they zero-length features?")

# split the features according to inter-breakpoint intervals
bed_dt_split <- bed_dt[queryHits(ov), 1:3, with = F]
bed_dt_split$V2 <- pmax(bed_dt_split$V2, start(int_gr)[subjectHits(ov)] - 1L)
bed_dt_split$V3 <- pmin(bed_dt_split$V3, end(int_gr)[subjectHits(ov)])
# store ID in 4th column
bed_dt_split$V4 <- seq_len(nrow(bed_dt_split))

# save the split BED file and liftOver it
write.table(bed_dt_split, file = ftemp, sep = "\t", quote = F, row.names = F, col.names = F)
system(paste0("liftOver -multiple -bedPlus=1 -tab ", ftemp, " ",
  balancer.dir, "/liftOver/", ain, "To", sub("^d", "D", aout), ".over.chain.gz ",
  fout_temp, " ",
  funmapped))
unlink(ftemp)

# add back missing features from columns 4+
restore_features <- function(f_temp, f)
{
  f_dt <- fread(f_temp, header = F, sep = "\t")
  unlink(f_temp)

  sel <- queryHits(ov)[f_dt$V4]
  f_dt$V4 <- paste0(bed_dt$V4[sel], "_", subjectHits(ov)[f_dt$V4])
  f_dt[, V5 := NULL]
  if (ncol(bed_dt) > 4)
    f_dt <- cbind(f_dt, bed_dt[sel, 5:ncol(bed_dt), with = F])

  write.table(f_dt, file = f, sep = "\t", quote = F, row.names = F, col.names = F)
}

restore_features(fout_temp, fout)
