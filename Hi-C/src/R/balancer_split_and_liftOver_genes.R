suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(GenomicRanges))

Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"
Alek.dir <- "/g/furlong/project/33_Hi-C"
gff.file <- "/g/furlong/genome/D.melanogaster/Dm6/6.05/gff/dmel-all-filtered-r6.05.UCSC_names.genes.gff.gz"
gff <- rtracklayer::import(gff.file)

system(paste0("ln -f -s ", gff.file, " ", Alek.dir, "/data/gff/"))

balancer.dir <- "/g/furlong/project/39_Balancer_Hi-C"
assembly <- "3"
br <- fread(paste0(balancer.dir, "/analysis/breakpoints_dm6_v", assembly, ".tab"), header = T)

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
gff_split$Parent <- CharacterList(mapply(function(v, ...) if (length(v) == 0) character(0) else paste0(v, ...), gff_split$Parent, "_", subjectHits(ov)))
start(gff_split) <- pmax(start(gff_split), start(int_gr)[subjectHits(ov)])
end(gff_split) <- pmin(end(gff_split), end(int_gr)[subjectHits(ov)])

rtracklayer::export.gff3(gff_split, paste0(Alek.dir, "/data/gff/dmel-all-filtered-r6.05.UCSC_names.split.genes.gff"))

system(paste0("liftOver -multiple -gff ", Alek.dir, "/data/gff/dmel-all-filtered-r6.05.UCSC_names.split.genes.gff ",
  balancer.dir, "/liftOver/dm6ToDm6bal", assembly, ".over.chain.gz ",
  Alek.dir, "/data/gff/dmel-all-filtered-r6.05-dm6bal", assembly, ".UCSC_names.split.genes.gff ",
  Alek.dir, "/data/gff/dmel-all-filtered-r6.05-unmapped-dm6bal", assembly, ".UCSC_names.split.genes.gff"))

unlink(paste0(Alek.dir, "/data/gff/dmel-all-filtered-r6.05.UCSC_names.split.genes.gff"))

system(paste0("gzip -f ", Alek.dir, "/data/gff/dmel-all-filtered-r6.05-dm6bal", assembly, ".UCSC_names.split.genes.gff"))
system(paste0("gzip -f ", Alek.dir, "/data/gff/dmel-all-filtered-r6.05-unmapped-dm6bal", assembly, ".UCSC_names.split.genes.gff"))
