library(data.table)

inv = fread(snakemake@input[["inv"]])
colnames(inv) = c("chrom", "start", "end", "id", "gt")
inv

reads = rbind(cbind(fread(snakemake@input[["inner"]]), where = "inner"),
              cbind(fread(snakemake@input[["left"]]),  where = "left"),
              cbind(fread(snakemake@input[["right"]]),  where = "right"))
colnames(reads) = c("chrom","start","end","inv","gt","rd.chrom","rd.start","rd.end","rd.strand", "where")
reads

summary = reads[, .(plus = sum(rd.strand=="+"), minus = sum(rd.strand=="-")),
                by = .(where,inv)]
summary = merge(summary, inv, by.x="inv", by.y="id")
summary

write.table(summary, file = snakemake@output[[1]], quote=F, sep="\t", col.names=T, row.names=F)
