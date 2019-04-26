library(data.table)
library(ggplot2)
library(scales)
source("scripts/colors.R")

args = commandArgs(trailingOnly=T)
d = as.data.table(fread(args[1]))
colnames(d) = c("chrom", "class", "count")

d[class == "1/1__1/1",]$class    = "common"
d[class == "0/1__1/1",]$class    = "wild type-specific"
d[class == "0/1__0/0",]$class    = "balancer-specific"
d[grepl("__0/1$",class),]$class  = "heterozygous"
d[grepl('[.01/]', class),]$class = "errorneous"
d = d[, .(count = sum(count)), by=.(chrom, class)]
d$class = factor(d$class, 
                 levels = c("balancer-specific", "wild type-specific", "common", "heterozygous", "errorneous"),
                 ordered = T)
d$chrom = factor(d$chrom, levels=c("chrX","chr3","chr2"), ordered = T)


#Chrom_lengths:
chrom_lengths = data.table(chrom = c("chr2", "chr3", "chrX"),
                           chrom_length = c(23513712+25286936, 28110227+32079331, 23542271))
d = merge(d, chrom_lengths, by = "chrom", all.x=T)


# Labels    
d <- d[, label := as.character(format(count,format="d", big.mark=","))][]
d <- d[, label_width := strwidth(label, 'inches')][]

# print numbers 
d[, .(chrom, class, count, chrom_length, perMb = count/chrom_length*1e6)]

# New
p <- ggplot(d) + 
    aes(class, count/chrom_length*1e6, fill=class) + 
    geom_bar(stat="identity", position = position_dodge(1)) + 
    coord_flip(expand = F, ylim = c(0, 1.2*max(d$count/d$chrom_length*1e6))) +
    ylab("SNVs per Megabase") +
    xlab("Chromosome") +
    facet_grid(chrom ~ ., switch = "y") +
    geom_text(aes(y = count/chrom_length*1e6 + -450 + 1350 * label_width,  # Hand-tuned to place the labels nicely!
                  label = label), 
               hjust = "right", size=3) +
    theme_minimal() +
    theme(axis.ticks.x = element_blank(), axis.text.y = element_blank()) +
    scale_fill_manual(values = c(`wild type-specific` = as.character(my_colors["wildtype"]),
                                 `balancer-specific`  = as.character(my_colors["balancer"]),
                                 `common`             = as.character(my_colors["common"]),
                                 `heterozygous`       = as.character(my_colors["hetero"]), 
                                 `errorneous`         = as.character(my_colors["error"]))) +
    guides(fill = guide_legend(reverse=T, title = NULL)) +
    scale_y_continuous(label=comma, limits=c(0,max(d$count/d$chrom_length*1e6) * 1.25))

f_out = sub('\\.pdf', '', args[2])
ggsave(p + theme(legend.position = "bottom"), 
       filename = args[2], width = 6, height = 4)
ggsave(p + theme(legend.position = "right"), 
       filename = paste0(f_out, '_2.pdf'), width = 7, height = 3)
