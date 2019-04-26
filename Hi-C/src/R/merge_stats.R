args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
{
  message("Usage:  Rscript merge_stats.R <outfile.stats> <infile1.stats> ... <infileN.stats>")
  q()
}

outfile <- args[1]
infiles <- tail(args, -1)

require(data.table)

dt <- data.table()

for (f in infiles)
{
  message(f)
  newdt <- fread(f, header = F, col.names = c("key", "value"))
  dt <- rbind(dt, newdt, fill = T)
}

dt <- dt[, list(value = sum(value, na.rm = T)), by = "key"]

write.table(dt, file = outfile, sep = "\t", quote = F, row.names = F, col.names = F)
