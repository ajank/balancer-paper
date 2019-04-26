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
    # Get all breakpoints
    regions = GRanges()
    sqlvls <- c("chr2L", "chr2R" ,"chr3L", "chr3R")
    for (chr in sqlvls)
    {
        x <- extract_breakpoints("dm6", chr, bin_size = 5e3, summarize="mean")
        x <- x[2:(length(x)-1)]
        regions = c(regions, GRanges(seqnames=rep(chr, length(x)), ranges = IRanges(x, x), 
                                     seqinfo = Seqinfo(sqlvls, genome="dm6")))
    }

	# Manually add extra regions (here: big DUPs)
	regions = c(regions, GRanges("chr2L", IRanges(9e6,9e6+1), seqinfo = Seqinfo(sqlvls, genome="dm6")))
	regions = c(regions, GRanges("chr2R", IRanges(3500000,3500001), seqinfo = Seqinfo(sqlvls, genome="dm6")))
    

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
        out[[length(out)+1]] <- list(regions[i], lleft[i], lright[i])
    }
    return(out)
}

bed <- fread("TADs_manual.dm6.bed")
dd <- with(bed, data.table(chrom = V1, pos = c(V2, V3)))
dd <- dd[!is.na(pos), ]
setorder(dd, chrom, pos)
dd <- unique(dd)
tad.lines.dm6 <- GRanges(dd$chrom, IRanges(dd$pos, dd$pos + 1L), seqinfo = Seqinfo(sqlvls, genome="dm6"))
tad.lines.dm6bal3 <- liftover(as.data.table(tad.lines.dm6)[,1:3, with=F])
regions <- prepare_regions_of_interest_as_list_of_lists()
