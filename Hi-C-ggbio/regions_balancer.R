sqlvls <- c("chr2L", "chr2R" ,"chr3L", "chr3R")

liftover <- function(dt)
{
    ftemp <- tempfile(c("in", "out", "tmp"))
    write.table(dt, file = ftemp[1], sep = "\t", quote = F, row.names = F, col.names = F)
    system(paste("/g/furlong/jankowsk/ucsc/liftOver", ftemp[1], liftover.file, ftemp[2], ftemp[3]))
    dd <- fread(ftemp[2])
    unlink(ftemp)
    return(GRanges(dd$V1, IRanges(dd$V2, dd$V3), seqinfo = Seqinfo(sqlvls, genome="dm6bal3")))
}


prepare_regions_of_interest_as_list_of_lists <- function()
{
    # Get all interesting SV loci
    sqlvls <- c("chr2L", "chr2R" ,"chr3L", "chr3R")
    r <- fread("/g/furlong/project/37_Capture-C/analysis/balancer/interesting_SV.tab", header = T)
    r <- r[allele != "chr3_Rend", ]
    regions = GRanges(r, seqinfo = Seqinfo(sqlvls, genome="dm6"))

    # Liftover coordinates left and right of the breakpoint
    dist = 10e3
    left <- as.data.table(regions)[,1:3, with=F]
    left[, start := start - 2 * dist]
    left[, end := end - dist]
    right <- as.data.table(regions)[,1:3, with=F]
    right[, start := start + dist]
    right[, end := end + 2 * dist]

    lleft <- liftover(left)
    lright <- liftover(right)

    # Make sure all coordinates could be lifted over
    assert_that(length(lleft) == length(regions))
    assert_that(length(lright) == length(regions))

    # Subtract the earlier shift by dist
    target <- 2L * start(lleft) - end(lleft)
    start(lleft) <- 0L
    start(lleft) <- end(lleft) <- target
    target <- 2L * start(lright) - end(lright)
    start(lright) <- 0L
    start(lright) <- end(lright) <- target
    
    out <- list()
    for (i in seq_along(regions)) {
        out[[length(out)+1]] <- list(
            # grs = list(regions[i]),
            grs = list(regions[i], lleft[i], lright[i]),
            prefix = paste0(regions$class[i], "_", regions$allele[i]),
            span = 300e3
        )
    }
    return(out)
}

regions <- prepare_regions_of_interest_as_list_of_lists()
