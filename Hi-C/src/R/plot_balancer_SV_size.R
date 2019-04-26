library(ggplot2)
library(data.table)
library(scales)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")
source("src/R/functions_DE_colors.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))


d = rbind(cbind(DEL_dt, `SV type` = "Deletions"),
          cbind(DUP_dt, `SV type` = "Duplications"))
d <- d[type != "common", ] # only allele-specific SVs
d <- merge(d, data.table(d[, list(label = paste0(`SV type`, " (", format(.N, big.mark = ","), ")")), by = `SV type`]))


pdf("analysis/balancer/SV_size.pdf", width = 4, height = 3)

plt <- ggplot(d) + 
    aes(x=(end-start+1) /1e3, y = ..count.., fill=`SV type`) + 
    facet_wrap(~ label, ncol = 1, scales = "free_y") +
    geom_histogram(binwidth=0.1, boundary = 0, position=position_dodge(), alpha=1) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=11, colour = "black"),
      axis.text.y = element_text(size=11, colour = "black"),
      strip.text = element_text(size=11, colour = "black"),
      legend.text = element_text(size=11, colour = "black")) +
    scale_x_log10(breaks = 10^(-2:2), labels = c("0.01","0.1","1","10","100")) +
    scale_y_continuous() +
    scale_fill_manual(NULL, values=c(Deletions = as.character(my_colors["del"]), 
                               Duplications = as.character(my_colors["dup"]))) +
    ggtitle(NULL) +
    theme(legend.position="none") +
    # theme(legend.position = c(1,1), legend.justification = c(1, 1), legend.background = element_rect(fill="white", color = NA)) +
    ylab("Count") + xlab("Structural variant size (kb)")
    #annotation_logticks(size = 0.2, short = unit(0.5,"mm"), mid = unit(1,"mm"), long = unit(2,"mm")) #+ 
    #coord_cartesian(ylim = c(1.55,3000), xlim = c(6, 290000), expand = F)
print(plt)

dev.off()


DEL_with_chrX_dt <- fread(evalq(FN.bed.sv.del, env.dm6), col.names = c("chrom", "start", "end", "name"))
DEL_with_chrX_dt[, start := start + 1L]

DUP_with_chrX_dt <- fread(evalq(FN.bed.sv.dup, env.dm6), col.names = c("chrom", "start", "end", "name"))
DUP_with_chrX_dt[, start := start + 1L]

d = rbind(cbind(DEL_with_chrX_dt, `SV type` = "Deletions"),
          cbind(DUP_with_chrX_dt, `SV type` = "Duplications"))
d <- d[!grepl("1/1_1/1_", name), ] # only allele-specific SVs
d <- merge(d, data.table(d[, list(label = paste0(`SV type`, " (", format(.N, big.mark = ","), ")")), by = c("SV type", "chrom")]))


pdf("analysis/balancer/SV_size_by_chrom.pdf", width = 4, height = 5.25)

plt <- ggplot(d) + 
    aes(x=(end-start+1) /1e3, y = ..count.., fill=`SV type`) + 
    facet_wrap(~ paste(label, substr(chrom, 1, 4)), ncol = 1, scales = "free_y") +
    geom_histogram(binwidth=0.1, boundary = 0, position=position_dodge(), alpha=1) +
    theme_minimal() +
    theme(axis.text.x = element_text(size=11, colour = "black"),
      axis.text.y = element_text(size=11, colour = "black"),
      strip.text = element_text(size=11, colour = "black"),
      legend.text = element_text(size=11, colour = "black")) +
    scale_x_log10(breaks = 10^(-2:2), labels = c("0.01","0.1","1","10","100")) +
    scale_y_continuous() +
    scale_fill_manual(NULL, values=c(Deletions = as.character(my_colors["del"]), 
                               Duplications = as.character(my_colors["dup"]))) +
    ggtitle(NULL) +
    theme(legend.position="none") +
    # theme(legend.position = c(1,1), legend.justification = c(1, 1), legend.background = element_rect(fill="white", color = NA)) +
    ylab("Count") + xlab("Structural variant size (kb)")
    #annotation_logticks(size = 0.2, short = unit(0.5,"mm"), mid = unit(1,"mm"), long = unit(2,"mm")) #+ 
    #coord_cartesian(ylim = c(1.55,3000), xlim = c(6, 290000), expand = F)
print(plt)

dev.off()
