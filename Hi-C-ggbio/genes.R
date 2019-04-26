# Need global variables:
assert_that(exists("FN.gff.FlyBase"))
assert_that(exists("FN.sql.reducedGFF"))


#'  input: GFF dataset (GenomicRanges object)
#'  output: GFF dataset with all the transcripts renamed from FBtrNNNNNNN to trFBgnNNNNNNN
gff_rename_transcripts_to_gene_name <- function(gff)
{
  sel <- grepl("^FBtr", gff$ID)
  tr_dt <- data.table(old_tr = as.vector(gff$ID[sel]), new_tr = paste0("tr", as.vector(gff$Parent[sel])))
  map_tr <- function(l)
  {
    v <- unlist(l)
    m <- match(v, tr_dt$old_tr)
    v[!is.na(m)] <- tr_dt$new_tr[na.omit(m)]
    return(CharacterList(lapply(relist(v, l), unique)))
  }
  gff$ID <- unlist(map_tr(gff$ID))
  gff$Parent <- map_tr(gff$Parent)
  # gff$Name[sel] <- unlist(gff$Parent[sel])
  return(gff)
}

#'  input: GFF dataset (GenomicRanges object), desired feature type (character) and desired Parent ID (character)
#'  output: data.table containing the union of all such features, with very basic metadata
gff_extract_reduced_features <- function(this_gff, this_type, this_Parent)
{
  gff_id_reduced <- reduce(this_gff[this_gff$type == this_type])
  if (length(gff_id_reduced) > 0)
  {
    gff_id_reduced$type <- this_type
    gff_id_reduced$Parent <- list(this_Parent)
    return(as.data.table(gff_id_reduced))
  }
  else
    return(NULL)
}

#'  input: GFF dataset (GenomicRanges object)
#'  output: GFF dataset where all the transcripts of the same ID are collapsed to a single transcript
#'  exons and CDSes in the reduced transcript are also collapsed by taking their union
gff_reduce_transcripts_of_identical_name <- function(gff)
{
  # start with entries without a parent (i.e. genes)
  reduced_dt <- as.data.table(gff[sapply(gff$Parent, function(v) length(v) == 0)])

  # iterate over all transcript IDs
  for (id in unique(na.omit(gff$ID[gff$type != "gene"])))
  {
    sel <- !is.na(gff$ID) & gff$ID == id
    gff_sel <- gff[sel]

    # extract the transcripts of a given ID, collapse them all
    if (sum(sel) == 1)
      reduced_dt <- rbind(reduced_dt, as.data.table(gff_sel), fill = T)
    else if (sum(sel) > 1)
    {
      gff_id_reduced <- reduce(gff_sel)
      stopifnot(length(unique(gff_sel$type)) == 1)
      gff_id_reduced$type <- gff_sel$type[1]
      gff_id_reduced$Parent <- list(unique(unlist(gff_sel$Parent)))
      gff_id_reduced$ID <- id
      reduced_dt <- rbind(reduced_dt, as.data.table(gff_id_reduced), fill = T)
    }

    # collapse seaparately each kind of children features
    sel <- sapply(gff$Parent, function(v) id %in% v)
    gff_sel <- gff[sel]
    reduced_dt_children <- rbind(
      gff_extract_reduced_features(gff_sel, "exon", id),
      # gff_extract_reduced_features(gff_sel, "five_prime_UTR", id),
      # gff_extract_reduced_features(gff_sel, "three_prime_UTR", id),
      gff_extract_reduced_features(gff_sel, "CDS", id)
    )
    reduced_dt <- rbind(reduced_dt, reduced_dt_children, fill = T)
  }

  reduced_gr <- with(reduced_dt, GRanges(seqnames, IRanges(start, end), strand))
  # for (col in setdiff(names(reduced_dt), c("seqnames", "start", "end", "width", "strand")))
  #   elementMetadata(reduced_gr)[, col] <- reduced_dt[[col]]
  reduced_gr$type <- reduced_dt$type
  reduced_gr$Name <- reduced_dt$Name
  reduced_gr$ID <- reduced_dt$ID
  reduced_gr$Parent <- CharacterList(reduced_dt$Parent)
  return(sort(reduced_gr))
}

#'  Get columns chrom, start, end, strand, type and gene_id
#'  from a gff object (imported via rtracklayer::import).
#'  Also, some weird entries are removed and chromosomes are
#'  filtered down to the main ones.
#'  Moreover, the column interval stores the identifier of
#'  between-breakpoint interval in the balancer assembly.
#'  In the case of reference assembly, the interval is NA.
extract_gff_into_data.table <- function(gff) 
{
    t1 <- data.table(chrom = as.vector(seqnames(gff)),
                     start = start(gff),
                     end   = end(gff),
                     strand = as.vector(strand(gff)))
    t2 <- as.data.table(mcols(gff))[, c("ID", "Name", "Parent", "type", "gene_id"), with=F]
    t <- cbind(t1, t2)
    # Remove bad categories (types)
    regex_bad <- 'point_mutation|breakpoint|polyA_site|insertion_site|rescue_fragment|sequence_variant|complex_substitution|deletion|DNA_motif'
    t <- t[!grepl(regex_bad, t$type), ]
    t$type <- factor(t$type)
    # Select only main chromosomes
    t <- t[grepl('^chr[234XY][LR]?$',t$chrom), ]

    # Determine the between-breakpoint interval for each feature
    t[, interval_from_ID := sapply(strsplit(ID, "_"), nth, 2)]
    f <- function(v)
    {
        r <- unique(sapply(strsplit(v, "_"), nth, 2, default = NA))
        stopifnot(length(r) < 2)
        if (length(r) == 0) NA else r
    }
    t[, interval_from_Parent := sapply(Parent, f)]
    stopifnot(t$interval_from_ID == t$interval_from_Parent | is.na(t$interval_from_ID) | is.na(t$interval_from_Parent))
    t[, interval := ifelse(is.na(interval_from_ID), interval_from_Parent, interval_from_ID)]
    stopifnot(all(is.na(t$interval)) || all(!is.na(t$interval)))
    t[, ID := NULL]
    t[, Parent := NULL]
    t[, gene_symbol := ifelse(type == "gene", Name, NA)]
    t[, Name := NULL]
    t[, interval_from_ID := NULL]
    t[, interval_from_Parent := NULL]

    t
}

#' I want to assign each gene to exactly one category.
#' To that end I order the possible labels by decreasing
#' significance and I only choose the first one for each 
#' gene_id.
#' Disadvantage is, that this will fail on non-hardcoded
#' types.
select_type_of_gene <- function(all_types) 
{
    lvls <- c("tRNA", "ncRNA", "snoRNA", "pre_miRNA", "snRNA", "rRNA",
              "pseudogene", "mRNA","gene", "five_prime_UTR", 
              "three_prime_UR", "intron", "exon", "CDS")
    ts <- factor(as.character(all_types), levels= lvls, ordered = T)
    as.character(dplyr::first(sort(ts)))
}


#' Add DE significance from a DESeq table
#' per gene and return as data.table
#' Give it a sample prefix to name the columns
add_DE_to_gene_table <- function(gff.tbl, de_file, alpha=0.05)
{
    # Read significant genes
    de <- read.table(de_file, header=T) %>%
        dplyr::select(gene_id, log2FoldChange, baseMean, padj) %>%
        as.data.table
    gff.sum <- merge(gff.tbl, de, all.x = T, by="gene_id")
    
    # factor c("significant", "insignificant", "not tested")
    na_rows <- is.na(gff.sum$padj)
    gff.sum$signf = "not tested"
    gff.sum$signf[!na_rows] = if_else(gff.sum$padj[!na_rows] < alpha, "significant", "insignificant")
    gff.sum$signf = factor(gff.sum$signf)
    
    # Replace NA values
    gff.sum$log2FoldChange[na_rows] = 0
    gff.sum$padj[na_rows] = 0.99
    gff.sum
}


#' Use a special GFF file to get the gene length (summed length
#' of allnon-overlapping exons)
get_gene_sizes <- function(file = FN.gff.FlyBase)
{
    # Reduce partly overlapping intervals and sum up their length
    get_total_exon_length <- function(start, end) {
        ir <- IRanges(start, end+1)
        return(sum(width(reduce(ir))))
    }
    
    d <- fread(paste("zcat -f", file), sep = "\t", header = F) %>%
        #dplyr::filter(V3 == "exon") %>% 
        dplyr::filter(V3 %in% c("CDS", "five_prime_UTR", "three_prime_UTR", "exon")) %>%
        dplyr::mutate(gene_id = sub('.*(FBgn[0-9]{7}).*', '\\1', V9)) %>%
        dplyr::select(chrom = V1, start = V4, end=V5, strand = V7, gene_id) %>%
        dplyr::group_by(gene_id) %>%
        dplyr::summarize(length = get_total_exon_length(start, end), start = min(start), end = max(end)) %>%
        as.data.table
    return(d)
}



#' Read HTSeq-count file for general gene expression
read_htseqcount <- function(f, gene_size = get_gene_sizes())
{
    expr_lvl <- fread(f, sep = "\t", header = F) %>%
        dplyr::rename(gene_id = V1, count = V2) %>%
        dplyr::filter(!grepl("^__", gene_id))
    assert_that(noNA(expr_lvl))
    total_read_count = sum(expr_lvl$count)
    expr_lvl <- suppressWarnings(dplyr::left_join(expr_lvl, gene_size, by="gene_id")) %>%
        mutate(FPKM = count*1000000/total_read_count*1000/length) %>%
        dplyr::select(gene_id, count, FPKM) %>% 
        filter(!is.na(FPKM))
    assert_that(noNA(expr_lvl))
    expr_lvl
}

#' Read multiple HTSeq-count files and return average gene expression
read_htseqcount_multi <- function(..., gene_size = get_gene_sizes())
{
    d <- NULL
    for (f in unlist(list(...)))
    {
        d <- rbind(cbind(read_htseqcount(f, gene_size = gene_size),
                         source = as.character(f)),
                   d)
    }
    d
}

#' Given multiple file names of HTSeq-count files,
#' plot a box plot (median) expression of the gene
plot_gene_expression <- function(df, ylims)
{
    assert_that(is.data.frame(df))
    assert_that('start_' %in% colnames(df))
    assert_that('end_' %in% colnames(df))
    assert_that('gene_id' %in% colnames(df))
    assert_that('interval' %in% colnames(df))
    assert_that('FPKM' %in% colnames(df))
    
    df$FPKM[df$FPKM < ylims[1]] = ylims[1]
    df$FPKM[df$FPKM > ylims[2]] = ylims[2]
    plt <- plot_background_polygons(df, ylims) +
        geom_boxplot(aes(x = (start_ + end_)/2, y = FPKM, group = paste(gene_id, interval)), data = df, inherit.aes = F, outlier.shape = NA) +
        geom_point(aes(x = (start_ + end_)/2, y = FPKM, group = paste(gene_id, interval)), data = df, inherit.aes = F) +
        scale_y_log10() + 
        guides(fill = F) + 
        coord_cartesian(ylim = ylims, expand=F ) +
        geom_hline(yintercept = .1, linetype = "dotted") + 
        geom_hline(yintercept = 10, linetype = "dotted") +
        ylab(NULL)

    return(plt)
}


#' Add FPKM values from HTSEQ-count files to the gene table.
#' If multiple files are given, average is used.
add_FPKM_to_gene_table <- function(gff.tbl, ...)
{
    files <- list(...)
    assert_that(is.list(files))
    assert_that(length(files)>0)
    
    d <- NULL
    for (f in unlist(files))
    {
        if (is.null(d))
        {
            d <- fread(f, sep = "\t", header = F)
        } else {
            d <- fread(f, sep = "\t", header = F) %>%
                merge(d, ., by = "V1", all = T)   
        }
    }
    
    
    
    de <- read.table(de_file, row.names = 1, header=T) %>%
        mutate(gene_id = rownames(.)) %>%
        dplyr::select(gene_id, log2FoldChange, padj) %>%
        as.data.table
    gff.sum <- merge(gff.tbl, de, all.x = T, by="gene_id")
}


#' Mapping of "gene space" (start, end), to 
#' evenly spaced columns (start_, end_)
distribute_coords <- function(d, start, end) 
{
    # d is expected to have columns "start" and "end"
    assert_that(is.data.frame(d))
    assert_that("start" %in% colnames(d))
    assert_that("end" %in% colnames(d))
    assert_that(dim(d)[1]>0)
    
    # output will go into "start_" and "end_"
    assert_that(! "start_" %in% colnames(d))
    assert_that(! "end_" %in% colnames(d))
    
    d$start_ = start + seq(0,dim(d)[1]-1)*(end - start)/(dim(d)[1])
    d$end_   = d$start_ + (end - start)/(dim(d)[1])*0.95
    d
}


#' Map an arbitrary position within [start, end]
#' to the gene space (given by distribute_coords() )
map_pos_to_gene_space <- function(df, bp) 
{
    assert_that(is.data.frame(df))
    assert_that(nrow(df)>0)
    assert_that('start' %in% colnames(df))
    assert_that('end' %in% colnames(df))
    assert_that('start_' %in% colnames(df))
    assert_that('end_' %in% colnames(df))
    
    assert_that(is.numeric(bp))
    assert_that(bp > min(df$start))
    assert_that(bp < max(df$end))

    df.ir <- IRanges(start = df$start, end = df$end)
    bp.ir <- IRanges(bp, bp)
    # first take the end coordinate of the preceding gene
    idx <- follow(bp.ir, df.ir, select = "last")
    coordinate <- df$end_[idx]
    # change to the start of the following gene, if there is any
    idx <- precede(bp.ir, df.ir, select = "first")
    id <- idx[!is.na(idx)]
    coordinate[!is.na(idx)] <- df$start_[id]
    # change it to the position within the gene, if we overlap any
    idx <- findOverlaps(bp.ir, df.ir, select = "first")
    id <- idx[!is.na(idx)]
    coordinate[!is.na(idx)] = df$start_[id] + (df$end_[id] - df$start_[id]) * (bp - df$start[id]) / (df$end[id] - df$start[id])
    stopifnot(length(bp) == length(coordinate))
    coordinate
    # print(coordinate)
}



#' I assing slightly different colors to neighboring genes,
#' so they are easier to distinguish in the plot. Significant
#' genes overrule this colour scheme
assign_bg_colors <- function(df) 
{
    assert_that("start" %in% colnames(df))
    assert_that("start_" %in% colnames(df))
    assert_that("signf" %in% colnames(df))
    assert_that( all( levels(df$signf) == c("insignificant", "not tested", "significant")) )
    assert_that(!is.unsorted(df$start))
    assert_that(!is.unsorted(df$start_))
    
    df$bg_col = factor(as.character(rep(1:5, nrow(df))[1:nrow(df)]),
                       levels = c("1","2","3","4","5","insignificant", "significant"))
    df$bg_col[df$signf != "not tested"] = as.character(df$signf[df$signf != "not tested"])
    df
}



#' geom_polygons to highlight columns belonging to genes
#' LFC or expression bar plots can be plotted above
plot_background_polygons <- function (df, ylims = c(-2,2), alpha = 0.3)
{
    assert_that("start_" %in% colnames(df))
    assert_that("end_" %in% colnames(df))
    assert_that("gene_id" %in% colnames(df))
    assert_that("interval" %in% colnames(df))
    assert_that("signf" %in% colnames(df))
    assert_that("bg_col" %in% colnames(df))
    polydata <- rbind(
        data.frame(x = df$start_, y = ylims[1], group = paste(df$gene_id, df$interval), bg_col = df$bg_col),
        data.frame(x = df$start_, y = ylims[2], group = paste(df$gene_id, df$interval), bg_col = df$bg_col),
        data.frame(x = df$end_,   y = ylims[2], group = paste(df$gene_id, df$interval), bg_col = df$bg_col),
        data.frame(x = df$end_,   y = ylims[1], group = paste(df$gene_id, df$interval), bg_col = df$bg_col) )
    ggplot(polydata) + aes(x, y, group = group, fill = bg_col) + geom_polygon(alpha = alpha)
}


#' Plot polygons that indicate the positional mapping from
#' genome coordinates to evenly spaced columns
plot_gene2column_mapping <- function(df, alpha = 0.3)
{
    assert_that("start_" %in% colnames(df))
    assert_that("end_" %in% colnames(df))
    assert_that("gene_id" %in% colnames(df))
    assert_that("interval" %in% colnames(df))
    assert_that("signf" %in% colnames(df))
    # Polygons
    e_poly <- df %>% 
        dplyr::select(gene_id, interval, p1 = start, p2 = end, p4 = end_,
                      p6 = start_, signf, bg_col) %>% 
        mutate(p3 = p2, p5 = p4, p7 = p6, p8 = p1) %>% 
        group_by(gene_id, interval) %>%
        summarize_all(funs(dplyr::first)) # Unique entry per gene_id ! plots are wrong otherwise
    e_poly$cc <- rep(c(1,2,3,4,5),nrow(df))[1:(nrow(e_poly))]
    e_poly <- e_poly %>%
        gather(source, x, -c(gene_id, interval, signf,cc, bg_col)) %>% 
        mutate(y=0) %>%
        arrange(gene_id, source)
    
    e_poly$y[e_poly$source == "p1"] = 1
    e_poly$y[e_poly$source == "p2"] = 1
    e_poly$y[e_poly$source == "p3"] = 0.96 - 0.06*e_poly$cc[e_poly$source == "p3"]
    e_poly$y[e_poly$source == "p4"] = 0.46 -0.06*e_poly$cc[e_poly$source == "p4"]
    e_poly$y[e_poly$source == "p5"] = 0
    e_poly$y[e_poly$source == "p6"] = 0
    e_poly$y[e_poly$source == "p7"] = 0.46 -0.06*e_poly$cc[e_poly$source == "p7"]
    e_poly$y[e_poly$source == "p8"] = 0.96 - 0.06*e_poly$cc[e_poly$source == "p8"]
    e_poly
    ggplot(e_poly) + 
        aes(x=x, y=y, group = paste(gene_id, interval), fill=bg_col) + 
        geom_polygon(alpha=alpha) + 
        theme_bw_no_border() + 
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) + 
        ylab(NULL) + guides(fill='none')
}


#'  geom_bar for LFC. outliers are trimmed of and replaced 
#'  by a symbol
plot_lfc <- function(df, ylims)
{
    assert_that("log2FoldChange" %in% colnames(df))
    assert_that("start_" %in% colnames(df))
    assert_that("end_" %in% colnames(df))
    assert_that("gene_id" %in% colnames(df))
    assert_that("interval" %in% colnames(df))
    assert_that("bg_col" %in% colnames(df))
    assert_that("signf" %in% colnames(df))
    
    df <- df %>%
        dplyr::select(gene_id, interval, start_, end_, log2FoldChange, bg_col, signf, gene_symbol) %>%
        group_by(gene_id, interval) %>%
        summarize_all(funs(dplyr::first))
    
    # plot point/character if out of range
    outliers <- df %>%
        filter(log2FoldChange < ylims[1] | log2FoldChange > ylims[2]) %>%
        mutate(log2FoldChange = ifelse(log2FoldChange < 0, ylims[1]*.95, ylims[2]*.95))
    
    # shorten bars to allowed range
    df$log2FoldChange[df$log2FoldChange < ylims[1]] = ylims[1]
    df$log2FoldChange[df$log2FoldChange > ylims[2]] = ylims[2]

    plt <- plot_background_polygons(df, ylims) +
        geom_bar(aes((start_+end_)/2, log2FoldChange, fill=bg_col), data=df, inherit.aes=F, stat="identity") + 
        coord_cartesian(ylim = ylims, expand=F) +
        geom_hline(yintercept=0, size=0.3) +
        geom_point(aes((start_+end_)/2, log2FoldChange, fill=signf), data=outliers, inherit.aes=F, shape=5, size=0.5) +
        scale_y_continuous(expand = c(0, 0)) +
        ylab(NULL) + 
        geom_text(aes(x=(start_+end_)/2, label=gene_symbol), data=df, inherit.aes=F, y=ylims[1]*.95,
                               hjust = 0, angle=90, size = 2) #  y=ylims[1]*.95, 
    
    return(plt)
}

hack_map_factors_to_other_factor <- function(fac, out_levels, reverse_levels = F) 
{
    if (is.factor(out_levels)) out_levels = levels(out_levels)
    out <- as.character(out_levels)
    if (reverse_levels)
        i <- -as.integer(factor(fac)) %% length(out) + 1L
    else
        i <- (as.integer(factor(fac)) - 1L) %% length(out) + 1L
    factor(out[i], levels = out)
}



#'  internal function to draw vertical lines as grid system annotation
drawDetails.vline_grob <- function(x, recording = TRUE)
{
  limits <- x$params$limits
  ymax <- x$params$ymax

  grid.polyline(
    x = unit((x$params$coords$x - limits[1]) / diff(limits), "npc"),
    y = unit((2 * x$params$coords$y - 1) / (2 * ymax - 1) / 2 + 0.5, "npc"),
    id = x$params$coords$id,
    gp = x$params$gp,
    vp = viewport(height = unit(2 * ymax - 1, "npc"), clip = "off")
  )
}

#'  Draw breakpoints as vertical lines across the given plot and above.
#'  xintercepts: position of lines (e.g. breakpoints)
#'  df: data.frame that contains the mapping between genome space and columns
#'     <- if df is given, the plot produces a zig-zag shape. Only applicable 
#'        to the genes-to-columns mapping panel
#'  ymax: how far up to draw (in multiplicities of current plot height)
#'  gp: gpar() graphical parameters to use while drawing
ggplot_arbitrary_vertical_lines <- function(which, xintercepts, df = NULL, ymax = 1, gp = NULL)
{
    stopifnot(class(which) == "GRanges")
    stopifnot(length(which) == 1)
    
    x <- xintercepts
    assert_that(is.numeric(x))
    
    if (!is.null(df))
    {
        assert_that(is.data.frame(df))
        assert_that(nrow(df)>0)
        assert_that('start' %in% colnames(df))
        assert_that('end' %in% colnames(df))
        assert_that('start_' %in% colnames(df))
        assert_that('end_' %in% colnames(df))
        assert_that(all(x > min(df$start)))
        assert_that(all(x < max(df$end)))
    }
    
    if (length(x)<1) return(NULL)
    else 
    {
        
        coords <- NULL
        for (i in seq_along(x))
        {
            # coordinates in df should match the ones calculated by plot_gene2column_mapping
            if (!is.null(df))
            {
                xbase <- map_pos_to_gene_space(df, x[i])
                coords <- rbind(coords, data.table(x = c(x[i], x[i], xbase, xbase), 
                    y = c(ymax, 0.96 - 0.06*3, 0.46 - 0.06*3, 0), id = i))
            } else {
                coords <- rbind(coords, data.table(x = c(x[i], x[i]), y=c(ymax, 0), id = i))
            }
        }
        
        limits <- c(start(which), end(which))
        grob <- grob(cl = "vline_grob")
        grob$params <- list(limits = limits, ymax = ymax, gp = gp, coords = coords)
        
        return(annotation_custom(grob = grob, xmin = limits[1], xmax = limits[2]))
    }
}

get_gene_sizes

subset_gene_data <- function(wh, dataset)
{
    assert_that(is.character(dataset))
    dataset = dataset[1]
    assert_that(dataset == "embryo" || dataset == "adult")
    
    if (dataset == "embryo") {
        genes.embryo.wh <- filter(genes.embryo,
                                  chrom == as.character(seqnames(wh)), 
                                  end >= start(wh), 
                                  start <= end(wh)) %>%
            arrange(start) %>%
            distribute_coords(., start(wh), end(wh)) %>%
            assign_bg_colors(.)
        expr.embryo.wh <- merge(genes.embryo.wh, 
                                expr.embryo, 
                                by = "gene_id", all.x = T)
        return(expr.embryo.wh)
    
    } else {
    
        genes.adult.wh <- filter(genes.adult,
                                 chrom == as.character(seqnames(wh)), 
                                 end >= start(wh), 
                                 start <= end(wh)) %>%
            arrange(start) %>%
            distribute_coords(., start(wh), end(wh)) %>%
            assign_bg_colors(.)
        expr.adult.wh <- merge(genes.adult.wh, 
                               expr.adult, 
                               by = "gene_id", all.x = T)
        return(expr.adult.wh)
    }    
}


# global data


### gene models
message("[global] loading reduced gff file...")
gff <- rtracklayer::import(FN.gff.FlyBase)

### reduced GFF (txdb)
result <- try(txdb <- AnnotationDbi::loadDb(FN.sql.reducedGFF), silent = T)
if (class(result) == "try-error")
{
    message("[global] preparing reduced gene models. This might take a bit...")
    gff.reduced <- gff %>%
        gff_rename_transcripts_to_gene_name %>%
        gff_reduce_transcripts_of_identical_name
    txdb <- GenomicFeatures::makeTxDbFromGRanges(gff.reduced, taxonomyId = 7227)
    AnnotationDbi::saveDb(txdb, file = FN.sql.reducedGFF)
}


select_nonNA_gene_symbol <- function(sym)
{
    candidates <- unique(sym[!is.na(sym)])
    if (length(candidates) == 0)
        NA
    else if (length(candidates) == 1)
        candidates
    else
        stop("non-unique gene symbol: ", candidates)
}

### Gene overview & DE (embryonic)
assert_that(exists("FN.txt.DESeq.embryos"))
assert_that(exists("FN.txt.DESeq.adults"))
assert_that(exists("FN.txt.expr.embryos"))
assert_that(exists("FN.txt.expr.adults"))
message("[global] loading DE gene data (embryo)...")
genes.embryo <- extract_gff_into_data.table(gff) %>%
    group_by(gene_id, interval) %>% 
    summarize(chrom = dplyr::first(chrom),
              start = min(start),
              end = max(end),
              strand = dplyr::first(strand),
              type = select_type_of_gene(type),
              gene_symbol = select_nonNA_gene_symbol(gene_symbol)) %>%
    as.data.table %>%
    add_DE_to_gene_table(., FN.txt.DESeq.embryos)


### Gene expression data (embryonic)
message("[global] loading expression gene data (embryo)...")
expr.embryo <- read_htseqcount_multi(FN.txt.expr.embryos)


### Gene overview & DE (adult)
message("[global] loading DE gene data (adult)...")
genes.adult <- genes.embryo %>%
    dplyr::select(-padj, -log2FoldChange, -signf) %>%
    add_DE_to_gene_table(., FN.txt.DESeq.adults)


### Gene expression data (adult)
message("[global] loading expression gene data (adult)...")
expr.adult <- read_htseqcount_multi(FN.txt.expr.adults)
