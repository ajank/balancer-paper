# Hack for Sascha to have the "good" packages  
.libPaths(c(
    "/g/furlong/jankowsk/R-lib/3.4.0-foss-2016b",
    "/g/easybuild/x86_64/CentOS/7/nehalem/software/R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0",
    "/g/easybuild/x86_64/CentOS/7/nehalem/software/R/3.4.0-foss-2016b/lib64/R/library"
))
# END

options(warn = 1)
suppressPackageStartupMessages(library(data.table, quietly = T))
suppressPackageStartupMessages(library(GenomicRanges, quietly = T))
suppressWarnings(suppressPackageStartupMessages(library(ggbio, quietly = T)))
suppressPackageStartupMessages(library(grid, quietly = T))
suppressPackageStartupMessages(library(lattice, quietly = T))
suppressPackageStartupMessages(library(scales, quietly = T)) # label=comma
suppressPackageStartupMessages(library(dplyr, quietly = T))
suppressPackageStartupMessages(library(tidyr, quietly = T))
suppressPackageStartupMessages(library(assertthat, quietly = T))
suppressPackageStartupMessages(library(AnnotationDbi, quietly = T))
suppressPackageStartupMessages(library(GenomicFeatures, quietly = T))
suppressPackageStartupMessages(library(org.Dm.eg.db, quietly = T))
suppressPackageStartupMessages(library(ggplot2, quietly = T))

Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"
Alek.dir <- "/g/furlong/project/33_Hi-C"
source(paste0(Alek.dir, "/src/R/plot_map_binned_functions.R"))
liftover.file <- "/g/furlong/project/39_Balancer_Hi-C/liftOver/dm6ToDm6bal3.over.chain.gz"
env.libs <- environment()
genomeInfo <- Seqinfo(c("chr2L", "chr2R" ,"chr3L", "chr3R"), 
                      genome="dm6", 
                      seqlengths = c(23513712, 25286936, 28110227, 32079331))






# Read regions from file (first argument):
args <- commandArgs(trailingOnly=T)
assert_that(length(args)==2)
regions_raw <- fread(args[1], select=1:4, col.names = c("chrom","start","end","sample"))
regions <- with(regions_raw, GRanges(chrom, IRanges(start, end), sample = sample, seqinfo = genomeInfo))

out_prefix = args[2]




# The following files load global data when they are sourced.
# Upon loading, most of the files report a short message to give 
# a quick impression of progress. Once all global data is loaded, 
# the main function subsets this data to the region of interest.
# The main function also heavily relies on plt_ functions that are
# defined in the sourced files, grouped by subject (e.g. TAD calls
# or genes).
# Note that there are cross-dependencies between the sourced files
# thus the order is important.
# files.dm6.R is a central file that stores all file paths, so
# it can be easily replaced with the dm6bal version later on.

### dm6 ENVIRONMENT
env.dm6 <- new.env(parent = env.libs)
env.dm6$genome <- "dm6"
source("files.dm6.R", local=env.dm6)
source("genes.R",     local=env.dm6)
source("tads.R",      local=env.dm6)
source("vlines.R",    local=env.dm6)
source("enhancers.R", local=env.dm6) 
source("rna.R",       local=env.dm6)
source("overview.R",  local=env.dm6)
source("main.R",      local=env.dm6)


### dm6bal3 ENVIRONMENT
env.dm6bal3 <- new.env(parent = env.libs)
env.dm6bal3$genome <- "dm6bal3"
source("files.dm6bal3.R", local=env.dm6bal3)
source("genes.R",     local=env.dm6bal3)
source("tads.R",      local=env.dm6bal3)
source("vlines.R",    local=env.dm6bal3)
source("enhancers.R", local=env.dm6bal3) 
source("rna.R",       local=env.dm6bal3)
source("overview.R",  local=env.dm6bal3)
source("main.R",      local=env.dm6bal3)


for (i in seq_along(regions))
{
    wh = regions[i]
    my.file = paste0(out_prefix,
                     as.character(seqnames(wh)), "_", 
                     round(start(wh)/1000000,2), "-",
                     round(end(wh)/1000000,2), "Mb_",
                     wh$sample, ".pdf")
    message("---", my.file, "---")
    
    # Add vlines around original coordinates, then increase window by 50%
    my.vlines = c(start(wh), end(wh))
    wh.width  = width(wh)
    if (wh.width < 50e3) wh.width = 50e3
    start(wh) = max(start(wh) - wh.width/2, 0)
    end(wh)   = min(end(wh) + wh.width/2, seqlengths(wh)[as.character(seqnames(wh))])

    cairo_pdf(my.file, width = 10*2, height = 7*2, onefile = F)
    tryCatch(
      print(evalq(main_plot(wh, wh$sample, "dm6", add_lines = my.vlines), env.dm6)),
      error = function(msg) {print(msg); print ("An assertError: dim(d)[1] not greater than 0 means that there are no genes in the region. Plot is ignored")})
    dev.off()
}
