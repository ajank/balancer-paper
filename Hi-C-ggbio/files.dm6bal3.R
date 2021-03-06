assert_that(exists("Sascha.dir"))
assert_that(exists("Alek.dir"))

### GFF
# the file below is generated by /g/furlong/project/33_Hi-C/src/R/balancer_split_and_liftOver_gff.R
FN.gff.FlyBase     <- paste0(Alek.dir, "/data/gff/dmel-all-filtered-r6.05-dm6bal3.UCSC_names.split.genes.gff.gz")
FN.sql.reducedGFF  <- paste0(Alek.dir, "/data/gff/dmel-all-filtered-r6.05-dm6bal3.UCSC_names.split.reduced.sqlite")

### Bed files for SVs and enhancers:
Sascha.enhancer.dir <- paste0(Sascha.dir, "/analyses/tracks_dm6bal3/enhancer/")
Sascha.tracks.dir <- paste0(Sascha.dir, "/analyses/tracks_dm6bal3/")
FN.bed.sv.dup <- paste0(Sascha.dir, "/analyses/SV_filter/FINAL/DUP.final.bed")
FN.bed.sv.del <- paste0(Sascha.dir, "/analyses/SV_filter/FINAL/DEL.final.bed")
FN.bed.enh.macs <- paste0(Sascha.enhancer.dir, "macs_summits_p50_slop50_merge25_clean.htseq.dm6.onlyChr2and3.bed")
FN.bed.enh.cad4 <- paste0(Sascha.enhancer.dir, "CAD4_plus_vienna_dm6.core80percent.onlyChr2and3.bed")

### TAD calls (external)
Alek.tad.dir <- paste0(Alek.dir, "/data/hic_annotations/dm6bal3")

### Arc. protein ChIP
Sascha.chip.dir <- paste0(Sascha.dir, "/analyses/tracks_dm6bal3/TADs/chipData/")
FN.bed.chip.ctcf  <- paste0(Sascha.chip.dir, "ChipSeq.CTCF.embryo_stage2.Boyle_2014.bed")
FN.bed.chip.suhw  <- paste0(Sascha.chip.dir, "ChipChip.SuHw.embryo_0-12h.White_K_27.bed")
FN.bed.chip.zw5   <- paste0(Sascha.chip.dir, "ChipChip.ZW5.embryo_2-4h.Karpen_G_5265.bed")
FN.bed.chip.beaf  <- paste0(Sascha.chip.dir, "ChipChip.BEAF-32.embryo_0-12h.White_K_21.bed")
FN.bed.chip.cp190 <- paste0(Sascha.chip.dir, "ChipChip.CP190.embryo_0-12h.White_K_22.bed")

### DE Genes
Sascha.de.dir   <- paste0(Sascha.dir, "/analyses/ase/deseq/")
Sascha.expr.dir <- paste0(Sascha.dir, "/analyses/gene_expression/counts/")
FN.txt.DESeq.embryos <- paste0(Sascha.de.dir, "DESeq.N1_6-8h.standardFormat.txt")
FN.txt.DESeq.adults  <- paste0(Sascha.de.dir, "DESeq.F1_heads.standardFormat.txt")
FN.txt.expr.embryos  <- paste0(Sascha.expr.dir, 
                               c("N1_pool_6-8h_rep1.htseq-count-rev.txt",
                                 "N1_pool_6-8h_rep2.htseq-count-rev.txt",
                                 "N1sex_pool_6-8h_rep1.htseq-count-rev.txt",
                                 "N1sex_pool_6-8h_rep2.htseq-count-rev.txt"))
FN.txt.expr.adults   <- paste0(Sascha.expr.dir, 
                               c("F1_CyOTM3_heads_rep1.htseq-count-rev.txt",
                                 "F1_CyOTM3_heads_rep2.htseq-count-rev.txt"))


### DHS data
FN.bedG.DHS <- paste0(Sascha.dir, "/analyses/dhs/tracks/6-8_WE.bedGraph.gz")


### RNA tracks
Sascha.rna.dir <- paste0(Sascha.dir, "/analyses/tracks_dm6bal3/rna/")
FN.bed.RNA.embryos.all.fwd <- paste0(Sascha.rna.dir, "N1_pool_6-8h_rep1.ALL.fwd.bedGraph.gz")
FN.bed.RNA.embryos.all.rev <- paste0(Sascha.rna.dir, "N1_pool_6-8h_rep1.ALL.rev.bedGraph.gz")
FN.bed.RNA.embryos.bal.fwd <- paste0(Sascha.rna.dir, "N1_pool_6-8h_rep1.ALT.fwd.bedGraph.gz")
FN.bed.RNA.embryos.bal.rev <- paste0(Sascha.rna.dir, "N1_pool_6-8h_rep1.ALT.rev.bedGraph.gz")
FN.bed.RNA.embryos.vrg.fwd <- paste0(Sascha.rna.dir, "N1_pool_6-8h_rep1.REF.fwd.bedGraph.gz")
FN.bed.RNA.embryos.vrg.rev <- paste0(Sascha.rna.dir, "N1_pool_6-8h_rep1.REF.rev.bedGraph.gz")
FN.bed.RNA.adults.all.fwd  <- paste0(Sascha.rna.dir, "F1_CyOTM3_heads_rep1.ALL.fwd.bedGraph.gz")
FN.bed.RNA.adults.all.rev  <- paste0(Sascha.rna.dir, "F1_CyOTM3_heads_rep1.ALL.rev.bedGraph.gz")
FN.bed.RNA.adults.bal.fwd  <- paste0(Sascha.rna.dir, "F1_CyOTM3_heads_rep1.ALT.fwd.bedGraph.gz")
FN.bed.RNA.adults.bal.rev  <- paste0(Sascha.rna.dir, "F1_CyOTM3_heads_rep1.ALT.rev.bedGraph.gz")
FN.bed.RNA.adults.vrg.fwd  <- paste0(Sascha.rna.dir, "F1_CyOTM3_heads_rep1.REF.fwd.bedGraph.gz")
FN.bed.RNA.adults.vrg.rev  <- paste0(Sascha.rna.dir, "F1_CyOTM3_heads_rep1.REF.rev.bedGraph.gz")



