# Hack for Sascha to have the "good" packages  
.libPaths(c(
    "/g/furlong/jankowsk/R-lib/3.4.0-foss-2016b",
    "/g/easybuild/x86_64/CentOS/7/nehalem/software/R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0",
    "/g/easybuild/x86_64/CentOS/7/nehalem/software/R/3.4.0-foss-2016b/lib64/R/library"
))
# END

options(warn = 1)
library(data.table, quietly = T)
library(GenomicRanges, quietly = T)
library(ggbio, quietly = T)
library(grid, quietly = T)
library(lattice, quietly = T)
library(scales, quietly = T) # label=comma
library(dplyr, quietly = T)
library(tidyr, quietly = T)
library(assertthat, quietly = T)
library(AnnotationDbi, quietly = T)
library(GenomicFeatures, quietly = T)
library(org.Dm.eg.db, quietly = T)

Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"
Alek.dir <- "/g/furlong/project/33_Hi-C"
source(paste0(Alek.dir, "/src/R/plot_map_binned_functions.R"))
liftover.file <- "/g/furlong/project/39_Balancer_Hi-C/liftOver/dm6ToDm6bal3.over.chain.gz"
env.libs <- environment()

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
if(!exists("env.dm6")) env.dm6 <- new.env(parent = env.libs)
env.dm6$genome <- "dm6"
source("files.dm6.R", local=env.dm6)


if (!exists("gff", env.dm6) ||
    !exists("genes.embryo", env.dm6) ||
    !exists("genes.adult", env.dm6) ||
    !exists("expr.embryo", env.dm6) ||
    !exists("expr.adult", env.dm6)) {
    source("genes.R",     local=env.dm6)
} else {
    message("[global] genes.R already loaded")
}

if (!exists("chip", env.dm6) ||
    !exists("all_TADs", env.dm6)) {
    source("tads.R",      local=env.dm6)
} else {
    message("[global] tads.R already loaded")
}

source("vlines.R",    local=env.dm6)

if (!exists("SV.dt", env.dm6) ||
    !exists("DHS.dt", env.dm6)) {
    source("enhancers.R", local=env.dm6)
} else {
    message("[global] enhancers.R already loaded")
}

if (!exists("rna.fwd", env.dm6) ||
    !exists("rna.rev", env.dm6) ||
    !exists("rna.bal.fwd", env.dm6) ||
    !exists("rna.bal.rev", env.dm6) ||
    !exists("rna.vrg.fwd", env.dm6) ||
    !exists("rna.vrg.rev", env.dm6) ||
    !exists("rna2.fwd", env.dm6) ||
    !exists("rna2.rev", env.dm6) ||
    !exists("rna2.bal.fwd", env.dm6) ||
    !exists("rna2.bal.rev", env.dm6) ||
    !exists("rna2.vrg.fwd", env.dm6) ||
    !exists("rna2.vrg.rev", env.dm6)) {
    source("rna.R",       local=env.dm6)
} else {
    message("[global] rna.R already loaded")
}

source("insulation.R",  local=env.dm6)
source("overview.R",  local=env.dm6)
source("main_balancer_loops.R",      local=env.dm6)




### dm6bal3 ENVIRONMENT
if(!exists("env.dm6bal3")) env.dm6bal3 <- new.env(parent = env.libs)
env.dm6bal3$genome <- "dm6bal3"
source("files.dm6bal3.R", local=env.dm6bal3)
if (!exists("gff", env.dm6bal3) ||
    !exists("genes.embryo", env.dm6bal3) ||
    !exists("genes.adult", env.dm6bal3) ||
    !exists("expr.embryo", env.dm6bal3) ||
    !exists("expr.adult", env.dm6bal3)) {
    source("genes.R",     local=env.dm6bal3)
} else {
    message("[global] genes.R already loaded")
}

if (!exists("chip", env.dm6bal3) ||
    !exists("all_TADs", env.dm6bal3)) {
    source("tads.R",      local=env.dm6bal3)
} else {
    message("[global] tads.R already loaded")
}

source("vlines.R",    local=env.dm6bal3)
if (!exists("SV.dt", env.dm6bal3) ||
    !exists("DHS.dt", env.dm6bal3)) {
    source("enhancers.R", local=env.dm6bal3)
} else {
    message("[global] enhancers.R already loaded")
}

if (!exists("rna.fwd", env.dm6bal3) ||
    !exists("rna.rev", env.dm6bal3) ||
    !exists("rna.bal.fwd", env.dm6bal3) ||
    !exists("rna.bal.rev", env.dm6bal3) ||
    !exists("rna.vrg.fwd", env.dm6bal3) ||
    !exists("rna.vrg.rev", env.dm6bal3) ||
    !exists("rna2.fwd", env.dm6bal3) ||
    !exists("rna2.rev", env.dm6bal3) ||
    !exists("rna2.bal.fwd", env.dm6bal3) ||
    !exists("rna2.bal.rev", env.dm6bal3) ||
    !exists("rna2.vrg.fwd", env.dm6bal3) ||
    !exists("rna2.vrg.rev", env.dm6bal3)) {
    source("rna.R",       local=env.dm6bal3)
} else {
    message("[global] rna.R already loaded")
}

source("insulation.R",  local=env.dm6bal3)
source("overview.R",  local=env.dm6bal3)
source("main_balancer_loops.R",      local=env.dm6bal3)


### Files independent of environment 
source("regions_balancer.R") # regions, tad.lines.dm6, tad.lines.dm6bal3
source("fix_ggbio_tracks.R")


prepare_regions_of_interest_as_list_of_lists <- function()
{
    # Get all interesting SV loci
    sqlvls <- c("chr2L", "chr2R" ,"chr3L", "chr3R")
    dt <- fread("../analysis/balancer/HiC_DB_6-8h_VRG_dm6_annotations.tsv")
    out <- list()

    # for (i in seq_len(nrow(dt)))
    for (i in c(1,4,7:10))
    {
        regions <- GRanges()
        regions <- c(regions, with(dt, c(
            GRanges(chr1[i], IRanges(x1[i] + 150e3 + 1, x2[i] - 150e3), seqinfo = Seqinfo(sqlvls, genome="dm6")),
            GRanges(chr2[i], IRanges(y1[i] + 150e3 + 1, y2[i] - 150e3), seqinfo = Seqinfo(sqlvls, genome="dm6"))
        )))

        # Liftover coordinates
        regions_lifted <- liftover(as.data.table(regions))
        # Make sure all coordinates could be lifted over
        assert_that(length(regions_lifted) == length(regions))

        out[[length(out)+1]] <- list(
            # grs = list(regions[i]),
            grs = list(regions[1], regions[2]),
            prefix = paste0("balancer_loop_", i),
            span = 500e3
        )
    }
    return(out)
}

regions <- prepare_regions_of_interest_as_list_of_lists()


tmp.file <- tempfile()

par <- list(span = 500e3, width = 6, height = 2.36)

for (i in seq_along(regions))
{
    grs <- regions[[i]]$grs
    prefix <- regions[[i]]$prefix
    # span <- regions[[i]]$span
    
    if(!is.list(grs) || length(grs)<1) { 
        message(paste("Skip entry",i)); 
    } else {
        # Create a single PDF for set of GR
        my.file = paste0(prefix, "_", as.character(seqnames(grs[[1]])), "_", 
                          round(start(grs[[1]])/1000000,2), "Mb_span",
                          round(par$span/1000), "kb.pdf")
        message("---", my.file, "---")
        plts <- list()

        cairo_pdf(tmp.file, width = par$width, height = par$height, onefile = T)
        for (wh in grs) {
            start(wh) = start(wh)
            end(wh) = end(wh)
            my.genome <- genome(wh)[1]
            for (sample in c("embryo"))
            {
                if (my.genome == "dm6") {
                    try(plts <- c(plts, evalq(main_plot(wh, sample, my.genome), env.dm6)))
                } else {
                    try(plts <- c(plts, evalq(main_plot(wh, sample, my.genome), env.dm6bal3)))
                }
            }
        }
        dev.off()

        cairo_pdf(my.file, width = par$width, height = par$height, onefile = T)
        for (p in plts)
        {
            grid.newpage()
            grid.draw(fix_ggbio_tracks(p))
        }
        dev.off()
    }
}

unlink(tmp.file)
