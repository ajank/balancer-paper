assert_that(exists("hack_map_factors_to_other_factor"))
assert_that(is.function(hack_map_factors_to_other_factor))

assert_that(exists("FN.bed.chip.ctcf"))
assert_that(exists("FN.bed.chip.suhw")) 
assert_that(exists("FN.bed.chip.zw5"))
assert_that(exists("FN.bed.chip.beaf"))
assert_that(exists("FN.bed.chip.cp190"))

assert_that(exists("Alek.dir"))
source(paste0(Alek.dir, "/src/R/plot_map_binned_functions.R"))


get_TADs <- function(files)
{
    d <- NULL
    for (i in seq_along(files))
    {
        f <- files[[i]]
        source <- names(files)[i]
        if (is.null(source) || source == "")
            source <- sub('^(.*/)?([^/]+).bed$', '\\2',f)

        gr <- rtracklayer::import(f)
        gr <- gr[grepl('^chr[234XY][LR]?$',seqnames(gr)), ]
        d <- data.frame(chrom = seqnames(gr), 
                        start = start(gr), 
                        end = end(gr), 
                        source = source,
                        color = if (grepl("_colors_", f) || grepl("Sexton_2012",f)) sapply(strsplit(gr$name, "_"), nth, 1) else NA,
                        even = factor(rep(c("odd", "even"), length(gr))[1:length(gr)])) %>%
            arrange(chrom, start) %>%
            rbind(d, .)
    }
    d %>% arrange(chrom, start) %>% as.data.table
}


#' Read ChIP chip data from bed files
#' and return as data.table
read_chip_peaks <- function(f, sample = "unknown") 
{
    bed <- fread(f, sep = "\t", header = F)
    data.table(chrom = bed[[1]],
               pos   = round((bed[[2]] + bed[[3]])/2),
               sample = sample)
}


### TADs
message("[global] loading TAD calls...")
all_TADs <- get_TADs(list(
    "Chromatin colors" = paste0(Alek.tad.dir, "/Filion_2010_colors_", genome, ".bed.gz"),
    "Hug...Vaquerizas 2017" = paste0(Alek.tad.dir, "/TADs_Hug_2017_embryo_3-4h_", genome, ".bed.gz"),
    "Ramírez...Manke 2017" = paste0(Alek.tad.dir, "/TADs_Ramirez_2017_Kc167_", genome, ".bed.gz"),
    "Eagen...Kornberg 2017" = paste0(Alek.tad.dir, "/TADs_Eagen_2017_Kc167_", genome, ".bed.gz"),
    "Eagen...Kornberg 2015" = paste0(Alek.tad.dir, "/TADs_Eagen_2015_salivary_glands_", genome, ".bed.gz"),
    "Hou...Corces 2012" = paste0(Alek.tad.dir, "/TADs_Hou_2012_Kc167_", genome, ".bed.gz"),
    "Sexton...Cavalli 2012"= paste0(Alek.tad.dir, "/TADs_Sexton_2012_embryo_16-18h_", genome, ".bed.gz"),
    # "IS 100 kb, Hi-C 2-4h"= paste0(Alek.tad.dir, "/TADs_HiC_WE_2-4h_IS100k_", genome, ".bed.gz"),
    # "IS 100 kb, Hi-C All"= paste0(Alek.tad.dir, "/TADs_HiC_DB_6-8h_All_IS100k_", genome, ".bed.gz")
    "HiCExplorer, Hi-C balancer"= paste0(Alek.tad.dir, "/TADs_HiC_DB_6-8h_BAL_HiCExplorer_", genome, ".bed.gz"),
    "HiCExplorer, Hi-C wild\u00adtype"= paste0(Alek.tad.dir, "/TADs_HiC_DB_6-8h_VRG_HiCExplorer_", genome, ".bed.gz"),
    "HiCExplorer, Hi-C All"= paste0(Alek.tad.dir, "/TADs_HiC_DB_6-8h_All_HiCExplorer_", genome, ".bed.gz")
))


### architectural proteins
message("[global] loading ChIP data...")
chip <-             read_chip_peaks(FN.bed.chip.ctcf , "CTCF")
chip <- rbind(chip, read_chip_peaks(FN.bed.chip.suhw , "su(Hw)"))
chip <- rbind(chip, read_chip_peaks(FN.bed.chip.zw5  , "ZW5"))
chip <- rbind(chip, read_chip_peaks(FN.bed.chip.beaf , "BEAF-32"))
chip <- rbind(chip, read_chip_peaks(FN.bed.chip.cp190, "CP190"))


plt_TADcalls <- function(wh)
{
    tad.wh <- all_TADs[  all_TADs$chrom == as.character(seqnames(wh))[1]
                         & all_TADs$start < end(wh)
                         & all_TADs$end > start(wh),]
    chip.wh <- chip[chrom == as.character(seqnames(wh))[1] & pos < end(wh) & pos > start(wh)]
    chip.wh[, new_fac := hack_map_factors_to_other_factor(sample, tad.wh$source, reverse_levels = T)]
    tad.wh[, source := factor(tad.wh$source)]
    p20 <- ggplot(tad.wh[is.na(color)]) +
        aes(x = start, xend = end, y = source, yend = source, alpha=even) + 
        geom_segment(size=3.5) +
        guides(alpha = 'none') + 
        ylab(NULL) +
        scale_y_discrete(drop = F, position = "right") +
        scale_alpha_manual(values = c(0.7, 0.4))

    col <- list(BLACK = "#000000", BLUE = "#0072B2", GREEN = "#009E73", RED = "#D55E00", YELLOW = "#F0E442", 
                Active = "darkolivegreen2", Null = "gray80", HP1_centromeric = "firebrick3", PcG = "deepskyblue3")
    for (name in names(col))
        p20 <- p20 + geom_segment(data = tad.wh[color == name], size = 3.5, alpha = 0.7, color = col[[name]])

    # p21 <- p20 + #geom_vline(data=chip.wh, aes(xintercept=pos, col = sample), size = 0.2) + 
    #     geom_point(data = chip.wh, 
    #                aes(x=pos, col=sample, y = new_fac), 
    #                inherit.aes = F, alpha = 0.7, shape=20,
    #                position = position_jitter(height=0.7)) +
    #     theme(legend.key.size = unit(10, "pt")) +
    #     ylab(NULL) + scale_color_brewer(type='qual', palette = 6) +
    #     # guides(colour = guide_legend(override.aes = list(size=3), title = NULL))
    #     guides(colour = NULL)
    return(p20)
}
