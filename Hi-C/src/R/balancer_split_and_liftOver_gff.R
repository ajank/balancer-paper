args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
{
  message("Usage:  Rscript balancer_split_and_liftOver_gff.R <input GFF file> <input assembly, e.g. dm6> <target assembly, e.g. dm6bal3> <output GFF file> <unmapped GFF file>")
  q()
}

fin <- args[1]
ain <- args[2]
aout <- args[3]

ftemp <- tempfile()

fout <- args[4]
funmapped <- args[5]

suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(GenomicRanges))

# read input BED file
gff <- rtracklayer::import(fin)

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
ov <- findOverlaps(gff, int_gr)
missing <- setdiff(seq_along(gff), queryHits(ov))
if (length(missing) > 0)
  stop("missing ", length(missing), " rows among the overlaps, are they zero-length features?")

# split the features according to inter-breakpoint intervals
gff_split <- gff[queryHits(ov)]
gff_split$ID <- ifelse(is.na(gff_split$ID), NA, paste0(gff_split$ID, "_", subjectHits(ov)))
start(gff_split) <- pmax(start(gff_split), start(int_gr)[subjectHits(ov)])
end(gff_split) <- pmin(end(gff_split), end(int_gr)[subjectHits(ov)])

# save the split BED file and liftOver it
rtracklayer::export.gff3(gff_split, ftemp)
system(paste0("liftOver -multiple -gff ", ftemp, " ",
  balancer.dir, "/liftOver/", ain, "To", sub("^d", "D", aout), ".over.chain.gz ",
  fout, " ",
  funmapped))
unlink(ftemp)
