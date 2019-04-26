library(data.table)

# get MEI candidates from a list of Yad (email from 22/07/2016)
d = fread("potential_meis.csv")
d_ = d[distance!="" | grepl("upstream|downstream|fused",Note, ignore.case = T)]
message("\nThere are ", nrow(d_), " candidate genes with a weird TSS (based on an early analysis by Yad)")

# get DE genes
e = fread("../../deseq/DESeq.N1_6-8h.standardFormat.txt")
d_e = d_[gene_id %in% e[padj < 0.05, gene_id]]
message("  - ", nrow(d_e), " of them are among the DE genes (", nrow(e[padj<0.05]), ")")

# which are considered MEI?
f = fread("../MEI.gene_ids.txt", header=F)
d_f = d_[gene_id %in% f$V1]
message("  --- ", nrow(d_f), " of those are verified MEI cases")

# which might be fused
d_fused = d_[grepl("fused", Note, ignore.case = T) & !(gene_id %in% f$V1)]
message("  --- another ", nrow(d_fused), " genes are fused to some gene and the DE signal might be caused by that*")
message("")
message("  * I don't know whether these are fused to MEI genes or not.")


