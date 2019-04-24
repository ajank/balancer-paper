library(data.table)
D = NULL
for (file in snakemake@input) {
  d = fread(file)
  D = rbind(D, data.table(file = file, 
                          mapped    = d[V1 == "total read pairs",V2],
                          balancer  = d[V1 == "balancer read pairs",V2],
                          wildtype  = d[V1 == "wild type read pairs",V2],
                          ambiguous = d[V1 == "ambiguous read pairs",V2],
                          error     = d[V1 == "errorneous read pairs",V2]))
}
D$perc_bal    = D$balancer / D$mapped
D$perc_wt     = D$wildtype / D$mapped
D$perc_phased = (D$balancer + D$wildtype) / D$mapped
D$perc_err    = D$error / D$mapped
write.table(D, file = snakemake@output[[1]], row.names = F, quote=F, sep = "\t", col.names = T)
