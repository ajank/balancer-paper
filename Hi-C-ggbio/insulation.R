bin_size <- 5000L

liftOver_bw <- function(d, assembly_from, assembly_to)
{
    if (assembly_from == assembly_to)
        return(d)

    ftmp1 <- tempfile()
    ftmp2 <- tempfile()
    ftmp3 <- tempfile()
    write.table(d, file = ftmp1, sep = "\t", quote = F, row.names = F, col.names = F)
    system(paste0("Rscript ", Alek.dir, "/src/R/balancer_split_and_liftOver_bed.R ", ftmp1, " ", assembly_from, " ", assembly_to, " ", ftmp2, " ", ftmp3))
    cat("\nShould be only deletions:\n")
    system(paste("cat", ftmp3))
    cat("\n")

    dl <- fread(ftmp2, sep = "\t", header = F)
    dl[, V4 := as.numeric(str_split_fixed(V4, "_", 2)[, 1])]
    unlink(c(ftmp1, ftmp2, ftmp3))
    names(dl) <- names(d)[seq_along(names(dl))]
    return(dl)
}

if (!exists("IS_separate"))
{
    IS_separate <- list()
    IS_combined <- list()
}

read_IS <- function(dataset, allele, source = "filtered_5000")
{
    prefix <- paste(dataset, allele)
    if (!is.null(IS_separate[[prefix]]))
        return()

    calculation_genome <- if (allele == "VRG") "dm6" else if (allele == "BAL") "dm6bal3" else NA
    fname <- paste0(Alek.dir, "/data/hicexplorer/h5/", calculation_genome, "_", dataset, "_", allele, "_", source, "_corrected_tad_score.bm")
    message("Reading ", fname)
    d <- fread(fname, header = F, sep = "\t", skip = 1)
    names(d)[1:3] <- c("chrom", "start", "end")

    # split the regions larger than bin_size
    sel <- d$end - d$start > bin_size
    ds <- d[sel, ]
    d$end[sel] <- d$end[sel] - bin_size
    ds[, start := end - bin_size]
    d <- rbind(d, ds)

    # take only the coordinate where the IS value is calculated, +/- 1 bp
    d[, end := start + 1L]
    d[, start := start - 1L]

    # convert to the genome assembly used for plotting
    d <- liftOver_bw(d, calculation_genome, genome)

    # convert to 1-based coordinates
    d[, start := start + 1L]

    dt <- NULL
    for (i in 4:ncol(d))
        dt <- rbind(dt, data.table(chrom = d$chrom, start = d$start, end = d$end, midpoint = d$start, value = d[[i]], group = i - 3L))
    setkey(dt, chrom, start)
    IS_separate[[prefix]] <<- dt
    IS_combined[[prefix]] <<- dt[, list(value = mean(value), group = 0L), by = c("chrom", "start", "end", "midpoint")]
}

if (!exists("PC1"))
    PC1 <- list()

read_PC1 <- function(dataset, allele, source = "filtered_5000")
{
    prefix <- paste(dataset, allele)
    if (!is.null(PC1[[prefix]]))
        return()

    calculation_genome <- if (allele == "VRG") "dm6" else if (allele == "BAL") "dm6bal3" else NA
    fname <- paste0(Alek.dir, "/data/hicexplorer/h5/", calculation_genome, "_", dataset, "_", allele, "_", source, "_corrected_PCA1.bedgraph")
    message("Reading ", fname)
    d <- fread(fname, header = F, sep = "\t", skip = 1)
    names(d)[1:4] <- c("chrom", "start", "end", "value")

    # split the regions larger than bin_size
    sel <- d$end - d$start > bin_size
    ds <- d[sel, ]
    d$end[sel] <- d$end[sel] - bin_size
    ds[, start := end - bin_size]
    d <- rbind(d, ds)

    # convert to the genome assembly used for plotting
    d <- liftOver_bw(d, calculation_genome, genome)

    # convert to 1-based coordinates
    d[, start := start + 1L]

    # add extra data points where the plot crosses 0
    sel <- which(d$chrom == c(tail(d$chrom, -1), NA) & d$value * c(tail(d$value, -1), NA) < 0)
    ds <- d[sel, ]
    next_start <- d$start[sel + 1L]
    next_value <- d$value[sel + 1L]
    ds[, start := (start * next_value - next_start * value) / (next_value - value)]
    ds[, value := 0]
    d <- rbind(d, ds)

    d[, midpoint := start]
    setkey(d, chrom, start)
    PC1[[prefix]] <<- d
}

plot_IS_track <- function(wh, allele)
{
    dataset <- "HiC_DB_6-8h_combined"
    read_IS(dataset, allele)
    prefix <- paste(dataset, allele)

    dt_separate <- IS_separate[[prefix]][chrom == as.character(seqnames(wh)) & start(wh) <= end & start <= end(wh) + bin_size]
    dt_combined <- IS_combined[[prefix]][chrom == as.character(seqnames(wh)) & start(wh) <= end & start <= end(wh) + bin_size]

    ggplot(dt_separate, aes(x = midpoint, y = value, group = group)) +
        geom_line(color = "#6f6f6f", alpha = 0.3) +
        geom_line(data = dt_combined, color = "dodgerblue4") +
        ylab(NULL)
}

plot_PC1_track <- function(wh, allele)
{
    dataset <- "HiC_DB_6-8h_combined"
    read_PC1(dataset, allele)
    prefix <- paste(dataset, allele)

    dt <- PC1[[prefix]][chrom == as.character(seqnames(wh)) & start(wh) <= end & start <= end(wh) + bin_size]

    ggplot(dt[value >= 0, ], aes(x = midpoint, y = value)) +
        geom_area(fill = "#cc3333") +
        geom_area(data = dt[value <= 0, ], fill = "#336699") +
        geom_hline(yintercept = 0, color = "#6f6f6f") +
        scale_y_continuous(position = "right") +
        coord_cartesian(ylim = c(-0.022, 0.022), expand = F) +
        ylab(NULL)
}

read_boundaries <- function(dataset, allele, subset = "all")
{
    dt <- fread("../analysis/balancer/differential_TADs.tab", header = T)

    if (allele == "VRG")
        dt <- dt[allele == "w", ]
    else if (allele == "BAL")
        dt <- dt[allele == "b", ]
    else
        stop("unknown allele: ", allele)

    if (subset == "matched")
        dt <- dt[class %in% c("matched", "shifted"), ]
    else if (subset == "specific")
        dt <- dt[!(class %in% c("matched", "shifted")), ]
    else if (subset != "all")
        stop("unknown subset: ", subset)

    # take +/- 1 bp in 0-based coordinates
    d <- with(dt, 
        data.table(chrom, start = start - 2L, end = end + 1L),
    )

    # convert to the genome assembly used for plotting
    d <- liftOver_bw(d, "dm6", genome)

    # get rid of +/- 1 bp
    d[, start := start + 1L]
    d[, end := end - 1L]
    stopifnot(d$end - d$start == 0)
    # convert to 1-based coordinates
    d[, start := start + 1L]

    setkey(d, chrom, start)
    return(d)
}

if (!exists("tad.lines.vrg"))
{
    tad.lines.vrg <- GRanges(read_boundaries("HiC_DB_6-8h_combined", "VRG", subset = "specific"))
    tad.lines.vrg.relaxed <- GRanges(read_boundaries("HiC_DB_6-8h_combined", "VRG", subset = "matched"))
}
if (!exists("tad.lines.bal"))
{
    tad.lines.bal <- GRanges(read_boundaries("HiC_DB_6-8h_combined", "BAL", subset = "specific"))
    tad.lines.bal.relaxed <- GRanges(read_boundaries("HiC_DB_6-8h_combined", "BAL", subset = "matched"))
}
