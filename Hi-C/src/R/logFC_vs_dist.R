require(cowplot)
require(ggplot2)
require(assertthat)
require(dplyr)

logFC_vs_dist <- function(
  genes,     # a complete list of genes incl. there ASE state
  center,    # regions to center. If attribute "midpoint" is provided, counting starts there, otherwise counting starts at their center.
  forbidden, # distance of gene to feature is not allowed to cross these regions (breakpoints!)
  mode = "closest", # "closest" or "all", see below
  max_dist = 100e3,
  extra_filtering = NULL
)
{
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(GenomicRanges))

  assert_that(is.data.table(genes))
  assert_that(all(c("gene_id","chrom","start","end","strand","log2FoldChange","padj") %in% colnames(genes)))
  assert_that(all(genes$strand %in% c("+", "-")))
  assert_that("GRanges" %in% class(center))
  assert_that("GRanges" %in% class(forbidden))
  assert_that(mode %in% c("closest", "all", "overlap"))
  stopifnot(as.vector(seqnames(center)) %in% unique(genes$chrom))

  tss_dt <- genes[, c("chrom", "start", "end")] # we do not want "strand" here, since it affects how "precede" and "follow" work below
  tss_dt[, start := ifelse(genes$strand == "+", start, end)]
  tss_dt[, end := ifelse(genes$strand == "-", end, start)]
  tss_dt$id <- seq_len(nrow(tss_dt))
  tss <- GRanges(tss_dt)

  if (mode != "overlap")
  {
    if (!("midpoint" %in% colnames(center)))
      center$midpoint <- (start(center) + end(center)) / 2
    start(center) <- center$midpoint
    end(center) <- center$midpoint - 1L
  }

  # Alternative A: assign each gene to its two closest features (one upstream of TSS and one downstream of TSS) plus the ones overlapping TSS
  if (mode == "closest")
  {
    # consider TSSes overlapping regions separately
    ov <- findOverlaps(center, tss, minoverlap = 0)

    # combine preceding, overlapping and following TSSes for all the regions
    # note that function names "precede" and "follow" are counterintuitive (refer to the documentation!)
    dt <- rbind(
      data.table(center_id = precede(tss, center), tss_id = seq_along(tss), ref_coords_direction = -1L),
      data.table(center_id = queryHits(ov), tss_id = subjectHits(ov), ref_coords_direction = rep(0L, length(ov))),
      data.table(center_id = follow(tss, center), tss_id = seq_along(tss), ref_coords_direction = 1L)
    )
    dt <- dt[!is.na(center_id), ]
  }

  # Alternative B: count each gene possibly to multiple features.
  else if (mode == "all")
  {
    ov <- findOverlaps(center, tss + max_dist, minoverlap = 0)
    dt <- data.table(center_id = queryHits(ov), tss_id = subjectHits(ov))
    dt[, ref_coords_direction := ifelse(start(tss)[tss_id] < start(center)[center_id], -1L,
      ifelse(end(center)[center_id] < end(tss)[tss_id], 1L,
      ifelse(start(center)[center_id] <= end(tss)[tss_id] & start(tss)[tss_id] <= end(center)[center_id], 0L, NA_integer_)))]
  }

  # Alternative C: count each gene at most once
  else
  {
    ov <- findOverlaps(center, tss, minoverlap = 0)
    dt <- data.table(center_id = queryHits(ov), tss_id = subjectHits(ov), ref_coords_direction = rep(0L, length(ov)))
  }

  # focus on the gap between the region and the TSS
  pg <- pgap(center[dt$center_id], tss[dt$tss_id])
  dt[, gap := width(pg)]
  # if the gap overlaps a forbidden region (breakpoint), then remove this match
  ov <- findOverlaps(pg, forbidden, minoverlap = 0L)
  dt <- dt[!(seq_len(nrow(dt)) %in% queryHits(ov)), ]

  if (mode == "closest")
    dt <- dt[gap <= max_dist, ]
  else
    assert_that(all(dt$gap <= max_dist))

  dt[, ref_coords_distance := gap * ref_coords_direction]
  dt[, strand := as.factor(genes$strand[tss_id])]
  dt[, directed_distance := ifelse(strand == "+", ref_coords_distance, -ref_coords_distance)]
  setkey(dt, tss_id, center_id)

  # Return: genes table (with possibly duplicate entries) with extra columns:
  # center_id: row number of the matching region
  # ref_coords_distance: distance to TSS; calculated in the direction of increasing reference genome coordinates
  # directed_distance: distance to TSS; negative if the gene is transcribed towards the region, negative otherwise.
  #   Note that overlapping TSSes have ref_coords_distance == directed_distance == 0.
  ret <- genes[dt$tss_id]
  ret[, center_id := dt$center_id]
  ret[, ref_coords_distance := dt$ref_coords_distance]
  ret[, directed_distance := dt$directed_distance]

  if (!is.null(extra_filtering) && extra_filtering == "bait")
  {
    # for Hi-C data: take only the differential contacts originating from the TSS of a given gene
    tss_gr <- with(ret, GRanges(chrom, IRanges(ifelse(strand == "+", start, end), ifelse(strand == "-", end, start))))
    bait_gr <- with(ret, GRanges(center$baitChr[center_id], IRanges(center$baitStart[center_id], center$baitEnd[center_id])))
    ret <- ret[width(pintersect(tss_gr, bait_gr)) > 0, ]
  }

  if (!is.null(extra_filtering) && extra_filtering == "gene_id")
  {
    # for Capture-C data: take only the differential contacts originating from the viewpoint associated with a given gene
    ret <- ret[!is.na(center$gene_id[center_id]) & center$gene_id[center_id] == gene_id, ]
  }

  return(ret)
}

# 1) Run logFC_vs_dist on all genes
# 2) Run logFC_vs_dist on only a subset of genes (i.e. deletions/duplications, breakpoints removed a.s.o.)

# dt1 <- logFC_vs_dist(genes.all, diffIS_upper_5percent, bp_gr, mode = "all", max_dist = 100e3)
# dt2 <- logFC_vs_dist(genes.all, diffIS_upper_5percent, bp_gr, mode = "closest", max_dist = 100e3)
# v1 <- with(dt1, paste(gene_id, center_id, directed_distance))
# v2 <- with(dt2, paste(gene_id, center_id, directed_distance))
# stopifnot(v2 %in% v1)

feature_vs_dist <- function(
  feature,   # GenomicRanges object with features
  center,    # regions to center. If attribute "midpoint" is provided, counting starts there, otherwise counting starts at their center.
  forbidden, # distance of gene to feature is not allowed to cross these regions (breakpoints!)
  mode = "closest", # "closest" or "all", see below
  max_dist = 100e3,
  extra_filtering = NULL
)
{
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(GenomicRanges))

  assert_that("GRanges" %in% class(feature))
  assert_that("GRanges" %in% class(center))
  assert_that("GRanges" %in% class(forbidden))
  assert_that(mode %in% c("closest", "all", "overlap"))
  stopifnot(as.vector(seqnames(center)) %in% as.character(seqnames(feature)))

  if (mode != "overlap")
  {
    if (!("midpoint" %in% colnames(center)))
      center$midpoint <- (start(center) + end(center)) / 2
    start(center) <- center$midpoint
    end(center) <- center$midpoint - 1L
  }

  # Alternative A: assign each gene to its two closest features (one upstream of feature and one downstream of feature) plus the ones overlapping feature
  if (mode == "closest")
  {
    # consider featurees overlapping regions separately
    ov <- findOverlaps(center, feature, minoverlap = 0)

    # combine preceding, overlapping and following featurees for all the regions
    # note that function names "precede" and "follow" are counterintuitive (refer to the documentation!)
    dt <- rbind(
      data.table(center_id = precede(feature, center), feature_id = seq_along(feature), ref_coords_direction = -1L),
      data.table(center_id = queryHits(ov), feature_id = subjectHits(ov), ref_coords_direction = rep(0L, length(ov))),
      data.table(center_id = follow(feature, center), feature_id = seq_along(feature), ref_coords_direction = 1L)
    )
    dt <- dt[!is.na(center_id), ]
  }

  # Alternative B: count each gene possibly to multiple features.
  else if (mode == "all")
  {
    ov <- findOverlaps(center, feature + max_dist, minoverlap = 0)
    dt <- data.table(center_id = queryHits(ov), feature_id = subjectHits(ov))
    dt[, ref_coords_direction := ifelse(start(feature)[feature_id] < start(center)[center_id], -1L,
      ifelse(end(center)[center_id] < end(feature)[feature_id], 1L,
      ifelse(start(center)[center_id] <= end(feature)[feature_id] & start(feature)[feature_id] <= end(center)[center_id], 0L, NA_integer_)))]
  }

  # Alternative C: count each gene at most once
  else
  {
    ov <- findOverlaps(center, feature, minoverlap = 0)
    dt <- data.table(center_id = queryHits(ov), feature_id = subjectHits(ov), ref_coords_direction = rep(0L, length(ov)))
  }

  # focus on the gap between the region and the feature
  pg <- pgap(center[dt$center_id], feature[dt$feature_id])
  dt[, gap := width(pg)]
  # if the gap overlaps a forbidden region (breakpoint), then remove this match
  ov <- findOverlaps(pg, forbidden, minoverlap = 0L)
  dt <- dt[!(seq_len(nrow(dt)) %in% queryHits(ov)), ]

  if (mode == "closest")
    dt <- dt[gap <= max_dist, ]
  else
    assert_that(all(dt$gap <= max_dist))

  dt[, ref_coords_distance := gap * ref_coords_direction]
  dt[, strand := as.factor(feature$strand[feature_id])]
  dt[, directed_distance := ifelse(strand == "+", ref_coords_distance, -ref_coords_distance)]
  setkey(dt, feature_id, center_id)

  # Return: feature table (with possibly duplicate entries) with extra columns:
  # center_id: row number of the matching region
  # ref_coords_distance: distance to feature; calculated in the direction of increasing reference genome coordinates
  # directed_distance: distance to feature; negative if the gene is transcribed towards the region, negative otherwise.
  #   Note that overlapping featurees have ref_coords_distance == directed_distance == 0.
  ret <- as.data.table(feature)[dt$feature_id]
  setnames(ret, "seqnames", "chrom")
  ret[, center_id := dt$center_id]
  ret[, ref_coords_distance := dt$ref_coords_distance]
  ret[, directed_distance := dt$directed_distance]

  if (!is.null(extra_filtering) && extra_filtering == "bait")
  {
    # for Hi-C data: take only the differential contacts originating from the feature of a given gene
    feature_gr <- with(ret, GRanges(chrom, IRanges(start, end)))
    bait_gr <- with(ret, GRanges(center$baitChr[center_id], IRanges(center$baitStart[center_id], center$baitEnd[center_id])))
    ret <- ret[width(pintersect(feature_gr, bait_gr)) > 0, ]
  }

  if (!is.null(extra_filtering) && extra_filtering == "gene_id")
  {
    # for Capture-C data: take only the differential contacts originating from the viewpoint associated with a given gene
    ret <- ret[!is.na(center$gene_id[center_id]) & center$gene_id[center_id] == gene_id, ]
  }

  return(ret)
}

make_feature_plot <- function(genes, gr, forbidden, ggtitle, mode = "closest", max_dist = 100e3, point = T, smooth = T, abs_dist = F, abs_fc = F, color_feature = F)
{
  dt <- logFC_vs_dist(genes, gr, forbidden, mode = mode, max_dist = max_dist)
  dt[, alpha := ifelse(signf == "i", 0.6, 1)]
  dtext <- logFC_vs_dist(genes, gr, forbidden, mode = mode, max_dist = 2 * max_dist)

  dist_transform_function <- ifelse(abs_dist, abs, identity)
  fc_transform_function <- ifelse(abs_fc, abs, identity)

  # default plotting
  if (abs_fc)
    p <- ggplot(dt[signf != "n", ], aes(dist_transform_function(ref_coords_distance / 1e3), fc_transform_function(log2FoldChange), color = signf, alpha = signf))
  else
    p <- ggplot(dt[signf != "n", ], aes(dist_transform_function(ref_coords_distance / 1e3), fc_transform_function(log2FoldChange), color = signf2, alpha = signf2))
  p <- p +
    scale_color_manual(values = gene.colors, labels = gene.labels, name = "Gene expression") +
    scale_alpha_manual(values = gene.alpha, labels = gene.labels, name = "Gene expression")

  if (point)
    p <- p + geom_point()

  # special plotting to color according to the breakpoint
  if (color_feature)
  {
    gr_names <- paste0(seqnames(gr), ": ", format(start(gr)/1e6, digits = 3, trim = T), " Mb")
    dt[, gr_name := factor(gr_names[center_id], gr_names)]
    p <- ggplot(dt[signf != "n", ], aes(dist_transform_function(ref_coords_distance / 1e3), fc_transform_function(log2FoldChange),
      shape = signf, color = gr_name)) +
      geom_point() +
      scale_shape_manual(values = c(16, 19), labels = gene.labels, name = "Gene expression") +
      scale_colour_hue(name = "Breakpoint")
  }

  # general features again
  p <- p + geom_vline(xintercept = 0, linetype = "dashed")

  if (smooth)
  {
    if (abs_fc)
      p <- p + geom_smooth(aes(fill = signf), data = dtext[signf != "n", ], method = "auto", show.legend = !point)
    else
      p <- p + geom_smooth(aes(fill = signf2), data = dtext[signf != "n", ], method = "auto", show.legend = !point)
  }

  p <- p +
    coord_cartesian(xlim = c(ifelse(abs_dist, 0, -100), 100)) +
    # scale_y_continuous(limits = c(ifelse(abs_fc, 0, -6), 6), expand = c(0.005, 0), oob = scales::squish, breaks = scales::pretty_breaks(n = 7)) +
    scale_y_continuous(limits = c(ifelse(abs_fc, 0, -10), 10), expand = c(0.005, 0), oob = scales::squish, breaks = scales::pretty_breaks(n = 7)) +
    scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    background_grid(major = "xy", minor = "xy") +
    # theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
    xlab("Genomic position (kb)") +
    ylab(ifelse(abs_fc, "Absolute expression log2FoldChange", "Expression log2FoldChange")) +
    ggtitle(if (is.null(ggtitle)) NULL else paste0(ggtitle, " (", format(length(gr), big.mark = ","), ")"))

  print(p)

  return(invisible(NULL))
}

extract_map_features <- function(genes, gr, forbidden, ggtitle, mode = "closest", max_dist = 100e3, bin_size = 5000L, abs_dist = F, abs_fc = F, fc_method = "most_significant", extra_filtering = NULL)
{
  dt <- logFC_vs_dist(genes, gr, forbidden, mode = mode, max_dist = max_dist, extra_filtering = extra_filtering)
  dt[, row_id := center_id]
  dt[, distance := abs(ref_coords_distance)]
  dt[, position := round((ref_coords_distance + bin_size / 2) / bin_size) * bin_size - bin_size / 2]

  if (!exists("glm_pvalues"))
    glm_pvalues <<- NULL
  # fit <- glm(signf == "s" ~ distance, data = dt[signf != "n", ], family = "binomial")
  # glm_pvalues <<- rbind(glm_pvalues, data.table(ggtitle = as.character(ggtitle), N = length(gr), formula = format(formula(fit)), t(coef(summary(fit))[2, ])))
  fit <- glm(signf == "s" ~ log10(distance + 1), data = dt[signf != "n", ], family = "binomial")
  glm_pvalues <<- rbind(glm_pvalues, data.table(ggtitle = as.character(ggtitle), N = length(gr), formula = format(formula(fit)), t(coef(summary(fit))[2, ])))

  amap_fill <- as.data.table(expand.grid(
    row_id = seq_along(gr),
    position = (-max_dist %/% bin_size + 1):(max_dist %/% bin_size) * bin_size - bin_size / 2
  ))
  amap_fill[, padj := 1]
  amap_fill[, log2FoldChange := 0]
  amap_merged <- rbind(amap_fill, dt, fill = T)

  if (fc_method == "most_significant")
    amap <- amap_merged[, list(log2FoldChange = log2FoldChange[which.min(padj)],
      count.signf = sum(signf == "s", na.rm = T), count.insignf = sum(signf == "i", na.rm = T)), by = c("row_id", "position")]
  else if (fc_method == "max_abs")
    amap <- amap_merged[, list(log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))],
      count.signf = sum(signf == "s", na.rm = T), count.insignf = sum(signf == "i", na.rm = T)), by = c("row_id", "position")]
  else
    stop("unknown fc_method: ", fc_method)

  if (abs_dist)
    amap[, position := abs(position)]
  if (abs_fc)
    amap[, log2FoldChange := abs(log2FoldChange)]

  return(amap)
}

extract_map_genes <- function(genes, gr, forbidden, ggtitle, mode = "closest", max_dist = 100e3, bin_size = 5000L, abs_dist = F, fc_method = "most_significant", extra_filtering = NULL)
{
  gene_id_by_row_id <- unique(genes$gene_id)
  dt <- logFC_vs_dist(genes, gr, forbidden, mode = mode, max_dist = max_dist, extra_filtering = extra_filtering)
  dt[, distance := abs(ref_coords_distance)]
  # below we multiply by -1 to convert from center-to-TSS distances to TSS-to-center distances
  dt[, position := round((-1 * directed_distance + bin_size / 2) / bin_size) * bin_size - bin_size / 2]

  if (sum(dt$signf == "i") > 0 && sum(dt$signf == "s") > 0)
  {
    suppressWarnings(ks <- with(dt, ks.test(directed_distance[signf == "s"], directed_distance[signf == "i"])))
    message(as.character(ggtitle), ": ", sum(genes$signf == "s"), " DE genes with avg. ", format(sum(dt$signf == "s") / sum(genes$signf == "s")),
      ", ", sum(genes$signf == "i"), " non-DE genes with avg. ", format(sum(dt$signf == "i") / sum(genes$signf == "i")),
      ", Kolmogorov-Smirnov p-value ", format(ks$p.value))
  }
  else
    message("skipped Kolmogorov-Smirnov testing")

  if (!("padj" %in% names(elementMetadata(gr))))
    gr$padj <- 0
  if (!("log2FoldChange" %in% names(elementMetadata(gr))))
    gr$log2FoldChange <- 0

  # collapse multiple TSSes into one heatmap row for a gene
  dt <- unique(dt[, c("gene_id", "signf", "center_id", "position")])
  dt[, row_id := match(gene_id, gene_id_by_row_id)]

  dt[, padj := gr$padj[center_id]]
  dt[, log2FoldChange := gr$log2FoldChange[center_id]]

  amap_fill <- as.data.table(expand.grid(
    row_id = seq_along(gene_id_by_row_id),
    position = (-max_dist %/% bin_size + 1):(max_dist %/% bin_size) * bin_size - bin_size / 2
  ))
  amap_fill[, signf := genes$signf[match(gene_id_by_row_id[row_id], genes$gene_id)]]
  amap_fill[, padj := 1]
  amap_fill[, log2FoldChange := 0]
  amap_merged <- rbind(amap_fill, dt, fill = T)

  if (fc_method == "most_significant")
    amap <- amap_merged[, list(log2FoldChange = log2FoldChange[which.min(padj)], count = .N - 1L),
      by = c("row_id", "signf", "position")]
  else if (fc_method == "max_abs")
    amap <- amap_merged[, list(log2FoldChange = log2FoldChange[which.max(abs(log2FoldChange))], count = .N - 1L),
      by = c("row_id", "signf", "position")]
  else
    stop("unknown fc_method: ", fc_method)

  if (abs_dist)
    amap[, position := abs(position)]

  return(amap)
}

sort_GRanges_by_width <- function(gr, decreasing = T)
{
  wd <- if ("duplication_size" %in% names(elementMetadata(gr))) -(gr$duplication_size) else width(gr)
  return(gr[order(wd, decreasing = decreasing)])
}

sort_GRanges_by_log2FoldChange <- function(gr, decreasing = T)
{
  return(gr[order(gr$log2FoldChange, decreasing = decreasing)])
}

make_heatmap_log2FoldChange_plot <- function(amap, gr, ggtitle, subtitle = NULL, max_dist = 100e3, bin_size = 5000L, ratio = 0.5, label = "Genomic position (kb)\n", legend_name = "Gene expression\nlog2 fold change")
{
  p <- ggplot(amap, aes(position / 1e3, row_id)) +
    geom_tile(aes(fill = pmin(pmax(log2FoldChange, -3), 3))) +
    # scale_fill_gradient2(low = "#336699", high = "#cc3333", name = expression(atop("Gene expression", log[2] * " fold change"))) +
    scale_fill_gradient2(limits = c(-3, 3), low = "#336699", high = "#cc3333", name = legend_name) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = theme_get()$legend.text) +
    xlab(label) +
    ylab(NULL) +
    coord_cartesian(xlim = c(-1, 1) * max_dist / 1e3) +
    scale_y_continuous(expand = c(0.02, 0)) +
    ggtitle(if (is.null(ggtitle)) NULL else paste0(ggtitle, " (", format(length(gr), big.mark = ","), ")"), subtitle = subtitle)

  if (is.finite(ratio))
    p <- p + coord_fixed(ratio = ratio)

  return(p)
}

make_heatmap_count_plot <- function(amap, gr, ggtitle, subtitle = NULL, max_dist = 100e3, bin_size = 5000L, ratio = 0.5, label = "Genomic position (kb)\n", legend_name = "Feature count", zmax = 3)
{
  p <- ggplot(amap, aes(position / 1e3, row_id)) +
    geom_tile(aes(fill = pmin(count, zmax))) +
    scale_fill_gradient(limits = c(0, zmax), low = "white", high = "#cc3333", name = legend_name) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.title = theme_get()$legend.text) +
    xlab(label) +
    ylab(NULL) +
    coord_cartesian(xlim = c(-1, 1) * max_dist / 1e3) +
    scale_y_continuous(expand = c(0.02, 0)) +
    ggtitle(if (is.null(ggtitle)) NULL else paste0(ggtitle, " (", format(length(gr), big.mark = ","), ")"), subtitle = subtitle)

  if (is.finite(ratio))
    p <- p + coord_fixed(ratio = ratio)

  return(p)
}

make_summap_features <- function(genes, gr, amap)
{
  summap.raw <- amap %>% 
    group_by(position)

  summap.do1 <- summap.raw %>%
    do(mean_cl_boot(.$count.signf)) %>%
    dplyr::select(position, count = y, count.lwr = ymin, count.upr = ymax) %>%
    mutate(signf = factor("s", levels(genes$signf)))

  summap.do2 <- summap.raw %>%
    do(mean_cl_boot(.$count.insignf)) %>%
    dplyr::select(position, count = y, count.lwr = ymin, count.upr = ymax) %>%
    mutate(signf = factor("i", levels(genes$signf)))

  return(rbind(summap.do1, summap.do2))
}

make_ribbon_plot_features <- function(genes, gr, summap, ggtitle, subtitle = NULL, max_dist = 100e3, bin_size = 5000L, ymax = 0, xlab = "Genomic position (kb)", ylab = "\nGene density")
{
  p <- ggplot(summap, aes(position / 1e3, count)) +
    # background_grid(major = "y", minor = "y") +
    geom_ribbon(aes(ymin = count.lwr, ymax = count.upr, fill = signf), alpha = 0.3) +
    geom_line(aes(color = signf)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(xlab) +
    ylab(ylab) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    coord_cartesian(ylim = c(0, max(summap$count.upr, ymax, na.rm = T))) +
    background_grid(major = "xy", minor = "x") +
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    ggtitle(if (is.null(ggtitle)) NULL else paste0(ggtitle, " (", format(length(gr), big.mark = ","), ")"), subtitle = subtitle)

  return(p)
}

make_ribbon_plot_features_vs_SV <- function(genes, gr, summap, ggtitle, subtitle = NULL, max_dist = 100e3, bin_size = 5000L, ymax = 0, xlab = "Genomic position (kb)", ylab = "\nGene density")
{
  SV.colors <- c(`SNV` = as.character(gene.colors["s"]),
    `DEL_bal` = as.character(gene.colors["b"]), `DEL_vrg` = as.character(gene.colors["w"]),
    `DUP_bal` = "#d95f02", `DUP_vrg` = "#7570b3")
  SV.labels <- c(`SNV` = "SNV",
    `DEL_bal` = "DEL_bal", `DEL_vrg` = "DEL_vrg",
    `DUP_bal` = "DUP_bal", `DUP_vrg` = "DUP_vrg")

  p <- ggplot(summap, aes(position / 1e3, count)) +
    # background_grid(major = "y", minor = "y") +
    geom_ribbon(aes(ymin = count.lwr, ymax = count.upr, fill = type), alpha = 0.3) +
    geom_line(aes(color = type)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(xlab) +
    ylab(ylab) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    coord_cartesian(ylim = c(0, max(summap$count.upr, ymax, na.rm = T))) +
    background_grid(major = "xy", minor = "x") +
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_manual(values = SV.colors, labels = SV.labels, name = NULL) +
    scale_fill_manual(values = SV.colors, labels = SV.labels, name = NULL) +
    ggtitle(if (is.null(ggtitle)) NULL else paste0(ggtitle, " (", format(length(gr), big.mark = ","), ")"), subtitle = subtitle)

  return(p)
}

make_ratio_plot_features <- function(genes, gr, summap, ggtitle, subtitle = NULL, max_dist = 100e3, bin_size = 5000L, ymax = 0, xlab = "Genomic position (kb)")
{
  summap <- summap %>% 
    group_by(position) %>%
    summarise(
      ratio = count[signf == "s"] / (count[signf == "s"] + count[signf == "i"]),
      ratio.lwr = count.lwr[signf == "s"] / (count.lwr[signf == "s"] + count.upr[signf == "i"]),
      ratio.upr = count.upr[signf == "s"] / (count.upr[signf == "s"] + count.lwr[signf == "i"])
    )

  p <- ggplot(summap, aes(position / 1e3, ratio)) +
    # background_grid(major = "y", minor = "y") +
    # geom_ribbon(aes(ymin = ratio.lwr, ymax = ratio.upr), fill = gene.colors["s"], alpha = 0.3) +
    geom_line(color = gene.colors["s"]) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(xlab) +
    ylab("Fraction of\nDE genes") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    # coord_cartesian(ylim = c(0, max(summap$ratio.upr, ymax, na.rm = T))) +
    coord_cartesian(ylim = c(0, max(summap$ratio, ymax, na.rm = T))) +
    background_grid(major = "xy", minor = "x") +
    # theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    ggtitle(if (is.null(ggtitle)) NULL else paste0(ggtitle, " (", format(length(gr), big.mark = ","), ")"), subtitle = subtitle)

  return(p)
}

make_ribbon_plot_genes <- function(genes, gr, amap, ggtitle, subtitle = NULL, max_dist = 100e3, bin_size = 5000L, ymax = 0, xlab = "Position from the TSS (kb)", ylab = "Average count\n(in a 5 kb bin)", ylabeller = waiver())
{
  amap[, signf := relevel(signf, "s")]

  gene.labels.counted <- gene.labels
  # gene.labels.counted <- list()
  # for (s in names(gene.labels))
  #   gene.labels.counted[[s]] <- paste0(gene.labels[[s]], " (", format(with(amap, length(unique(row_id[signf == s]))), big.mark = ","), ")")

  summap.raw <- amap %>% 
    group_by(position, signf)
 
  summap <- summap.raw %>%
    do(mean_cl_boot(.$count)) %>%
    dplyr::select(position, signf, count = y, count.lwr = ymin, count.upr = ymax)

  p <- ggplot(summap, aes(position / 1e3, count)) +
    # background_grid(major = "y", minor = "y") +
    geom_ribbon(aes(ymin = count.lwr, ymax = count.upr, fill = signf), alpha = 0.3) +
    geom_line(aes(color = signf)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(xlab) +
    ylab(ylab) +
    background_grid(major = "xy", minor = "x") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3), labels = ylabeller) +
    coord_cartesian(xlim = c(-1, 1) * max_dist / 1e3, ylim = c(0, max(summap$count.upr, ymax, na.rm = T))) +
    scale_color_manual(values = gene.colors, labels = gene.labels.counted, name = NULL) +
    scale_fill_manual(values = gene.colors, labels = gene.labels.counted, name = NULL) +
    ggtitle(ggtitle, subtitle = subtitle)

  return(p)
}

make_log2FoldChange_plot <- function(gr, xmax = 2, label = "log2 fold change")
{
  dt <- data.table(
    row_id = if (is.data.table(gr)) seq_len(nrow(gr)) else seq_along(gr),
    signf = if ("signf" %in% names(gr)) gr$signf else "i",
    log2FoldChange = gr$log2FoldChange
  )

  p <- ggplot(dt, aes(x = row_id, y = log2FoldChange)) +
    background_grid(major = "x", minor = "x") +
    # geom_hline(yintercept = 0, linetype = "dashed") +
    geom_col(fill = gene.colors["i"], color = NA, width = 1) +
    # geom_col(aes(fill = signf), color = NA, width = 1) +
    # scale_fill_manual(values = gene.colors, guide = F) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    xlab(NULL) +
    ylab(label) +
    scale_x_continuous(position = "top", expand = c(0.02, 0)) +
    # scale_y_reverse(limits = max(abs(gr$log2FoldChange), xmax) * c(1, -1), breaks = scales::pretty_breaks(n = 3)) +
    scale_y_continuous(limits = max(abs(gr$log2FoldChange), xmax) * c(-1, 1), breaks = scales::pretty_breaks(n = 3)) +
    coord_flip()

  return(p)
}

make_width_plot <- function(gr, xmax = 0, label = "Size (kb)")
{
  wd <- if ("duplication_size" %in% names(elementMetadata(gr))) -(gr$duplication_size) else width(gr)
  dt <- data.table(row_id = seq_along(gr), width = wd)

  p <- ggplot(dt, aes(x = row_id, y = width / 1e3)) +
    background_grid(major = "x", minor = "x") +
    geom_col(fill = gene.colors["i"], color = NA, width = 1) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    xlab(NULL) +
    ylab(label) +
    scale_x_continuous(position = "top", expand = c(0.02, 0)) +
    scale_y_reverse(limits = c(max(wd, xmax), min(wd, 0)) / 1e3, breaks = scales::pretty_breaks(n = 3)) +
    coord_flip()

  return(p)
}

make_composite_plot_features <- function(genes, gr, forbidden, mode = "all", max_dist = 100e3, bin_size = 5000L, abs_dist = F, abs_fc = F, fc_method = "most_significant", ggtitle = NULL, subtitle = NULL, sample_size = 200L, ymax = 0, side_max = NULL, side_label = "Feature size (kb)\n", extra_filtering = NULL)
{
  amap <- extract_map_features(genes, gr, forbidden, ggtitle, mode = mode, fc_method = fc_method, extra_filtering = extra_filtering)
  summap <- make_summap_features(genes, gr, amap)
  p1a <- make_ribbon_plot_features(genes, gr, summap, ggtitle, xlab = NULL, ymax = ymax)
  # p1b <- make_ratio_plot_features(genes, gr, summap, ggtitle = NULL, ymax = 0.2, xlab = NULL)

  size <- min(sample_size, length(gr))
  sel <- sort(sample(length(gr), size))
  gr.sel <- gr[sel]
  elementMetadata(gr.sel) <- elementMetadata(gr)[sel, , drop = F]

  map <- rep(0L, length(gr))
  map[sel] <- seq_len(size)

  amap.sel <- copy(amap)
  amap.sel[, row_id := map[row_id]]
  amap.sel <- amap.sel[row_id > 0, ]

  p2 <- make_heatmap_log2FoldChange_plot(amap.sel, gr, NULL, ratio = NA)
  if ("log2FoldChange" %in% names(elementMetadata(gr)))
    p3 <- make_log2FoldChange_plot(gr.sel, xmax = if (is.null(side_max)) 2.0 else side_max, label = side_label)
  else
    p3 <- make_width_plot(gr.sel, xmax = if (is.null(side_max)) 2000 else side_max, label = side_label)

  p12 <- align_plots(p1a, p2, align = "v")
  rw <- c(0.25, 1, -0.75, 4.25)
  p10 <- plot_grid(NULL, NULL, NULL, p12[[1]], nrow = 1, align = "h", rel_widths = rw) # labels = c(NA, "A")
  p23 <- plot_grid(NULL, p3, NULL, p12[[2]], nrow = 1, align = "h", rel_widths = rw) # labels = c("B", NA),
  p <- plot_grid(p10, p23, ncol = 1, align = "v", rel_heights = c(0.75, 1.75))

  return(p)
}

make_composite_plot_features_with_ratio <- function(genes, gr, forbidden, mode = "all", max_dist = 100e3, bin_size = 5000L, abs_dist = F, abs_fc = F, fc_method = "most_significant", ggtitle = NULL, subtitle = NULL, sample_size = 200L, ymax = 0, side_max = NULL, side_label = "Feature size (kb)\n", extra_filtering = NULL)
{
  amap <- extract_map_features(genes, gr, forbidden, ggtitle, mode = mode, fc_method = fc_method, extra_filtering = extra_filtering)
  summap <- make_summap_features(genes, gr, amap)
  p1a <- make_ribbon_plot_features(genes, gr, summap, ggtitle, xlab = NULL, ymax = ymax)
  p1b <- make_ratio_plot_features(genes, gr, summap, ggtitle = NULL, ymax = 0.2, xlab = NULL)

  size <- min(sample_size, length(gr))
  sel <- sort(sample(length(gr), size))
  gr.sel <- gr[sel]
  elementMetadata(gr.sel) <- elementMetadata(gr)[sel, , drop = F]

  map <- rep(0L, length(gr))
  map[sel] <- seq_len(size)

  amap.sel <- copy(amap)
  amap.sel[, row_id := map[row_id]]
  amap.sel <- amap.sel[row_id > 0, ]

  p2 <- make_heatmap_log2FoldChange_plot(amap.sel, gr, NULL, ratio = NA)
  if ("log2FoldChange" %in% names(elementMetadata(gr)))
    p3 <- make_log2FoldChange_plot(gr.sel, xmax = if (is.null(side_max)) 2.0 else side_max, label = side_label)
  else
    p3 <- make_width_plot(gr.sel, xmax = if (is.null(side_max)) 2000 else side_max, label = side_label)

  p12 <- align_plots(p1a, p1b, p2, align = "v", axis = "lr")
  rw <- c(0.25, 1, -0.75, 4.25)
  p1a0 <- plot_grid(NULL, NULL, NULL, p12[[1]], nrow = 1, align = "h", rel_widths = rw) # labels = c(NA, "A")
  p1b0 <- plot_grid(NULL, NULL, NULL, p12[[2]], nrow = 1, align = "h", rel_widths = rw) # labels = c(NA, "A")
  p23 <- plot_grid(NULL, p3, NULL, p12[[3]], nrow = 1, align = "h", rel_widths = rw) # labels = c("B", NA),
  p <- plot_grid(p1a0, p1b0, p23, ncol = 1, align = "v", rel_heights = c(0.75, 0.5, 1.75))

  return(p)
}

make_composite_plot_genes <- function(genes, gr, forbidden, mode = "all", max_dist = 100e3, bin_size = 5000L, abs_dist = F, abs_fc = F, fc_method = "most_significant", ggtitle = NULL, subtitle = NULL, sample_size = 200L, side_max = NULL, ribbon_label = "Average count\n(in a 5 kb bin)", side_label = "Gene expression\nlog2 fold change", legend_name = "Feature count", extra_filtering = NULL, ylabeller = waiver(), zmax = 3)
{
  gmap <- extract_map_genes(genes, gr, forbidden, ggtitle, mode = mode, fc_method = fc_method, extra_filtering = extra_filtering)
  p1 <- make_ribbon_plot_genes(genes, gr, gmap, ggtitle, max_dist = max_dist, xlab = NULL, ylab = ribbon_label, ylabeller = ylabeller)

  size <- min(sample_size, sum(genes$signf == "s"))
  sel <- sort(sample(seq_len(nrow(genes))[genes$signf == "s"], size))
  genes.sel <- genes[sel, ]

  map <- rep(0L, nrow(genes))
  map[sel] <- seq_len(size)

  gmap.sel <- copy(gmap)
  gmap.sel[, row_id := map[row_id]]
  gmap.sel <- gmap.sel[row_id > 0, ]

  if ("log2FoldChange" %in% names(elementMetadata(gr)))
    p2 <- make_heatmap_log2FoldChange_plot(gmap.sel, genes.sel, NULL, max_dist = max_dist, ratio = NA, label = "Position from the TSS (kb)\n", legend_name = legend_name)
  else
    p2 <- make_heatmap_count_plot(gmap.sel, genes.sel, NULL, max_dist = max_dist, ratio = NA, label = "Position from the TSS (kb)\n", legend_name = legend_name, zmax = zmax)
  p3 <- make_log2FoldChange_plot(genes.sel, xmax = if (is.null(side_max)) max(abs(genes_tested$log2FoldChange)) else side_max, label = side_label)

  p12 <- align_plots(p1, p2, align = "v")
  rw <- c(0.25, 1, -0.75, 4.25)
  p10 <- plot_grid(NULL, NULL, NULL, p12[[1]], nrow = 1, align = "h", rel_widths = rw) # labels = c(NA, "A")
  p23 <- plot_grid(NULL, p3, NULL, p12[[2]], nrow = 1, align = "h", rel_widths = rw) # labels = c("B", NA),
  p <- plot_grid(p10, p23, ncol = 1, align = "v", rel_heights = c(0.75, 2))

  return(p)
}

make_composite_plot_features_vs_SV <- function(SV, gr, forbidden, mode = "all", max_dist = 100e3, bin_size = 5000L, abs_dist = F, abs_fc = F, fc_method = "most_significant", ggtitle = NULL, subtitle = NULL, sample_size = 200L, ymax = 0, side_max = NULL, side_label = "Feature size (kb)\n", extra_filtering = NULL, ylab = "\nSV density")
{
  SV <- as.data.table(SV)
  if (!("chrom" %in% colnames(SV)))
    SV[, chrom := seqnames]
  SV[, gene_id := NA]
  SV[, strand := "+"]
  SV[, log2FoldChange := 1]
  SV[, padj := 0]
  SV[, signf := factor("s")]

  summap <- NULL
  for (t in unique(SV$type))
  {
    amap <- extract_map_features(SV[type == t], gr, forbidden, ggtitle, mode = mode, fc_method = fc_method, extra_filtering = extra_filtering)
    summap <- rbind(summap, data.table(type = t, make_summap_features(SV[type == t], gr, amap))[signf == "s"])
  }
  p1a <- make_ribbon_plot_features_vs_SV(SV, gr, summap, ggtitle, xlab = NULL, ymax = ymax, ylab = ylab)
  # p1b <- make_ratio_plot_features(SV, gr, summap, ggtitle = NULL, ymax = 0.2, xlab = NULL)

  return(p1a)
}
