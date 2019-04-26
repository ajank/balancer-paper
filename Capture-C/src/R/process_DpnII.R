options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 1)
  stop("Usage:  Rscript process_DpnII.R genome")

genome <- options[1]

require(data.table)

message("genome: ", genome)
chrom_lst <- switch(genome,
  dm3 = c("chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet", "chr3R", "chr3RHet", "chr4", "chrX", "chrXHet", "chrYHet"),
  dm6 = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"),
  stop("unknown genome assembly")
)
chrom_sizes_file <- switch(genome,
  dm3 = "/g/furlong/genome/D.melanogaster/Dm3/chrsizes/Dmelanogaster_UCSC_dm3_chrsizes.txt",
  dm6 = "/g/furlong/genome/D.melanogaster/Dm6/chrsizes/dm6.ucsc.chrom.sizes"
)

txt <- fread(paste0("analysis/digest_DpnII_", genome, ".txt"))
txt <- txt[Chromosome %in% chrom_lst, ]
# reorder fragments according to chromosome order in chrom_lst
txt <- txt[order(match(Chromosome, chrom_lst)), ]

frag <- data.table(
  frag_id = 1:nrow(txt),
  chrom = factor(txt$Chromosome, chrom_lst),
  start = txt$Fragment_Start_Position,
  end = txt$Fragment_End_Position
)

chrom_sizes <- fread(chrom_sizes_file, col.names = c("chrom", "size"))
chrom_sizes <- chrom_sizes[chrom %in% chrom_lst, ]
# reorder chromosomes according to their order in chrom_lst
chrom_sizes <- chrom_sizes[order(match(chrom, chrom_lst)), ]

write.table(frag, file = paste0("analysis/digest_DpnII_", genome, ".tab"), sep = "\t", quote = F, row.names = F)

bed <- data.table(chrom = frag$chrom, start = frag$start - 1L, end = frag$end, name = frag$frag_id, score = 0L, strand = ".")
write.table(bed, file = paste0("analysis/digest_DpnII_", genome, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)

save(frag, chrom_lst, chrom_sizes, file = paste0("analysis/digest_DpnII_", genome, ".Rdata"))
