library(ggplot2)
library(data.table)
library(scales)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")
source("src/R/functions_DE_colors.R")

# Combine SVs
d = rbind(cbind(SNV_dt, `Variant type` = "SNVs"),
          cbind(DEL_dt, `Variant type` = "Deletions"),
          cbind(DUP_dt, `Variant type` = "Duplications"),
          # cbind(MEI_dt, `Variant type` = "Mobile Element Insertions"),
          cbind(DHS_deleted_dt, `Variant type` = "DHS peaks >=5% deleted"), fill = T)
stopifnot(grepl('^chr[23]', d$chrom))

# Size classes:
d$addon = ""
d[`Variant type` == "Deletions" & end-start+1 < 50,]$addon = " (<50 bp)"
d[`Variant type` == "Deletions" & end-start+1 >= 50 & end-start+1 < 160,]$addon = " (50-159 bp)"
d[`Variant type` == "Deletions" & end-start+1 >= 160,]$addon = " (>=160 bp)"
d <- d[, `Variant type`:= paste0(`Variant type`, addon)]
d <- d[, `Variant type`:= factor(`Variant type`, levels = c("DHS peaks >=5% deleted", "Mobile Element Insertions", "Duplications", "Deletions (>=160 bp)", "Deletions (50-159 bp)", "Deletions (<50 bp)", "SNVs"), ordered = T)][]
d[, type := factor(type, c("vrg", "bal", "common"))]

# Summary table
print(table(d$type, d$`Variant type`))


d <- d[type != "common", ]
e = merge(d[, .(count = nrow(.SD)), by = .(`Variant type`, type)], 
          d[, .(total = nrow(.SD)), by = .(`Variant type`)], 
          by = c("Variant type"))


pdf("analysis/balancer/SV_ratio.pdf", width = 5.5, height = 2.7)

plt <- ggplot(e) + 
    aes(x = `Variant type`, y = count/total, fill=type) + 
    geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=11, colour = "black"),
      axis.text.y = element_text(size=11, colour = "black"),
      legend.text = element_text(size=11, colour = "black")) +
    theme(legend.position = "right") +
    theme(panel.grid.major.y = element_blank()) +
    scale_fill_manual(NULL, values = c(`vrg` = as.character(my_colors["wildtype"]),
              `bal`  = as.character(my_colors["balancer"]),
              `common`    = as.character(my_colors["common"])),
              labels = c(`vrg` = "wild\u00adtype", `bal` = "balancer")) +
    coord_flip() +
    xlab(NULL) +
    scale_y_continuous(labels = percent, expand = c(0, 0)) +
    ylab("Fraction (common variation excluded)") +
    ggtitle(NULL) +
    geom_text(data = e[, .(count = sum(count)), by = `Variant type`], 
              aes(x = `Variant type`, label = paste("n =",format(count,big.mark=",",trim=T))), 
              y = 0.97, hjust = 1, col = "white", size = 3, inherit.aes = F)
#+
#    scale_x_log10(breaks = c(10,100,1000,10000,100000), labels=c("10bp","100bp","1kb","10kb","100kb")) +
#    scale_y_log10(breaks = c(1,2,11,101,1001), labels=c("None",1,10,100,1000),
#                  minor_breaks = c(6,51,501)) +
#    scale_fill_manual(values=c(DEL="dodgerblue1", DUP="red1")) +
#    ggtitle("CNV size distribution") +
#    theme(legend.position = c(0.85,0.75), legend.background = element_rect(fill="white", size=0.2)) +
#    ylab("Count") + xlab("SV size")
print(plt)

dev.off()
