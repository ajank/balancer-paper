options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 1)
  stop("Usage:  Rscript extract_statistics.R collection")

collection <- options[1]

# collection <- "TS_Capture"

require(data.table)

ss <- NULL

datasets <- gsub("\\.nodups\\.stats$", "", list.files(path = paste0("data/", collection, "/bam/filtered_reads"), pattern = ".*\\.nodups\\.stats"))

dummy <- data.table(data.frame(key = character(0), value = character(0)))

for (dataset in datasets)
{
  print(dataset)

  s1 <- dummy
  f <- paste0("data/", collection, "/bam/filtered_reads/", dataset, ".stats")
  if (file.exists(f))
  {
    s1 <- fread(f, col.names = c("key", "value"))
    s1 <- s1[!grepl("^chrom_freq/|^dist_freq/", key), ]

    lev <- c("total", paste0("pair_types/", c("NN", "NM", "NC", "NL", "MM", "MC", "ML", "CC", "CL", "CX", "LL", "DD")), "total_mapped")
    s1 <- s1[order(factor(s1$key, lev)), ]
    s1[, key := paste0("fastq/", key)]
    s1[, value := as.character(value)]
  }

  s2 <- dummy
  f <- paste0("data/", collection, "/bam/filtered_reads/", dataset, ".dedup.stats")
  if (file.exists(f))
  {
    s2 <- fread(f, col.names = c("key", "value"))
    s2[, key := paste0("total_mapped/", key)]
  }

  s3 <- dummy
  f <- paste0("data/", collection, "/bam/filtered_reads/", dataset, ".nodups.stats")
  if (file.exists(f))
  {
    s3 <- fread(f, col.names = c("key", "value"))
    s3 <- s3[!grepl("^chrom_freq/|^dist_freq/|^total_", key), ]

    s3 <- s3[order(factor(s3$key, lev)), ]
    s3[, key := paste0("n_nodups/", key)]
  }

  s4p <- fread(paste0("data/hicexplorer/qc/", dataset, "_filtered_rs/QC_table.txt"), header = T)
  s4p <- s4p[, !grepl("File", names(s4p)), with = F]
  key <- names(s4p)
  s4 <- data.table(key, value = as.integer(s4p[1, ]))

  lev <- c("Pairs considered", "One mate low quality", "One mate not unique", "One mate unmapped", "Pairs mappable, unique and high quality",
    "One mate not close to rest site", "dangling end", "duplicated pairs", "same fragment", "self circle", "self ligation (removed)", "Pairs used",
    "inter chromosomal", "short range < 20kb", "long range", "inward pairs", "outward pairs", "left pairs", "right pairs")
  s4 <- s4[order(factor(s4$key, lev)), ]
  s4[, key := paste0("n_nodups/", key)]

  s5 <- dummy
  f <- paste0("data/", collection, "/fastq/", gsub("^dm[^_]*_", "", dataset), ".adapter.stats")
  if (file.exists(f))
  {
    s5 <- fread(f, col.names = c("key", "value"))
    s5[, key := paste0("adapter/", key)]
  }

  s <- rbind(s1, s2, s3, s4, s5)
  setnames(s, "value", dataset)

  if (is.null(ss))
    ss <- s
  else
    ss <- merge(ss, s, by = "key", all = T, sort = F)
}

write.table(ss, paste0("analysis/", collection, "/combined_stats.tsv"), sep = "\t", quote = F, row.names = F)
