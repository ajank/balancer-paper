# functionalize y-abstraction

require(data.table, quietly = T)
require(grid, quietly = T)
require(lattice, quietly = T)
require(gplots, quietly = T)
lattice.options(panel.error = NULL)

if (!exists("master.dir"))
  master.dir <- "/g/furlong/project/33_Hi-C"
balancer.dir <- "/g/furlong/project/39_Balancer_Hi-C"

chrom_lst <- c('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX')
chrom_sep_bins <- 10L
breakpoint_summarize <- "all"

loops <- list()


# if (exists("t_hat"))
#   grid.lines((c(0, rep(t_hat, each = 2)) + 0.5) * bin_size1, (c(rep(t_hat, each = 2), max(t_hat)) + 0.5) * bin_size1,
#     default.units = "native", gp = gpar(col = "red"))

help_usage <- function()
{
  cat('

Usage:  make_triangle_plot(dataset, genome, chrom, start, end, bin_size, [options])
        make_plot(dataset1, dataset2, genome, chrom, start, end, bin_size, [options])
        make_plot(dataset1, dataset2, genome, chrom, start, end, bin_size1, bin_size2, [options])
        export_track(dataset, genome, chrom, position, bin_size, [prefix])
        export_tracks_from_bedfile(dataset, genome, bin_size, bedfile)

')
}

help_options <- function()
{
  cat('

Additional options:  pdf.dim = <figure dimensions in inches> [default: 8]
                     what = "observed" [default] / "p-value" / "expected" / "ratio"
                     filter = "none" [default] / "onesignificant" / "bothsignificant"
                     distance_fit = <dataset from which the distance decay fit should be used>
                     colscale = "orangered" [default] / "red" / "rainbow"
                     zmax = <numeric value for maximum value plotted>
                     orientation = "up" [default] / "down"
                     mark_loops = F [default] / T
                     mark_breakpoints = F [default] / T

')
}

mypanel.levelplot.different_bin_sizes <- function(x, y, z, row.values, column.values, secondary, secondary.midpoints, subscripts, panel.internal, background = NA, ...)
{
  stopifnot(current.panel.limits()$xlim == current.panel.limits()$ylim)
  if (!is.na(background))
    panel.fill(background)
  # panel.internal(x, y, z, subscripts, ...)
  panel.internal(x, rev(current.panel.limits()$ylim[2] - y + current.panel.limits()$ylim[1]), z, subscripts, ...)
  sx <- rep(secondary.midpoints, ncol(secondary))
  # sy <- rep(rev(secondary.midpoints), each = nrow(secondary))
  sy <- rep(current.panel.limits()$ylim[2] - secondary.midpoints + current.panel.limits()$ylim[1], each = nrow(secondary))
  panel.internal(sx, sy, t(secondary), seq_along(secondary), ...)
  grid.segments(0, 1, 1, 0)
}

mypanel.levelplot.same_bin_sizes <- function(x, y, z, row.values, column.values, secondary = NULL, secondary.midpoints = NULL, panel.internal, background = NA, ...)
{
  if (!is.na(background))
    panel.fill(background)
  if (length(secondary) > 0)
  {
    panel.internal(x, y, ifelse(lower.tri(secondary), t(secondary), z), ...)
    grid.segments(0, 1, 1, 0)
  }
  else
    panel.internal(x, y, z, ...)
}

mypanel.levelplot.triangle <- function(x, y, z, panel.internal, background = NA, mark_loops = F, mark_breakpoints = F, mark_domains = F, ..., dataset1, genome, this_chrom, start, end, bin_size)
{
  if (!is.na(background))
    grid.polygon(c(0, 1, 1), c(1, 0, 1), gp = gpar(col = NA, fill = background))
  panel.internal(x, y, z, ...)

  if (mark_loops)
    triangle_plot_mark_loops(dataset1, genome = genome, this_chrom = this_chrom, start = start, end = end, bin_size = bin_size, ...)
  if (mark_breakpoints)
    triangle_plot_mark_breakpoints(genome = genome, this_chrom = this_chrom, start = start, end = end, bin_size = bin_size, ...)
  if (identical(mark_domains, T) || is.character(mark_domains))
    triangle_plot_mark_domains(dataset1 = dataset1, domain_source = mark_domains, genome = genome, this_chrom = this_chrom, start = start, end = end, bin_size = bin_size, ...)
}

mypanel.levelplot.compartment <- function(x, y, z, row.values, column.values, secondary, secondary.midpoints, subscripts, panel.internal, background = NA, secondary.zmax = 5, ...)
{
  if (!is.na(background))
    panel.fill(background)
  panel.internal(x, y, z, subscripts, ...)
  sx <- rep(secondary.midpoints, ncol(secondary))
  sy <- rep(rev(secondary.midpoints), each = nrow(secondary))
  at <- -25:25 / 25 * secondary.zmax
  colfun <- colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))
  panel.internal(sx, sy, pmax(pmin(t(secondary), secondary.zmax - 1e-6), -secondary.zmax + 1e-6),
    seq_along(secondary), at = at, col.regions = colfun(length(at) - 1))
  grid.segments(0, 1, 1, 0)
}

ensure_loaded_v3 <- function(dataset, genome, bin_size)
{
  bin_size <- as.integer(bin_size)
  if (!exists("maps"))
    maps <<- list()
  if (is.null(maps[[paste0(dataset, "_", genome)]]))
    maps[[paste0(dataset, "_", genome)]] <<- list()
  if (is.null(maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]]))
  {
    fname <- paste0(master.dir, "/data/rda/pdb_", dataset, "_", genome, "_", bin_size, ".rda", 3)
    cat(paste0("Loading ", fname, "\n"))
    map <- NULL
    load(fname)
    maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]] <<- map
  }
}

ensure_loaded_compartment <- function(dataset, genome, bin_size, type = 'A')
{
  id <- paste0(dataset, ifelse(type == 'A', '', paste0("_", type)))
  if (!exists("pc1s"))
    pc1s <<- list()
  if (is.null(pc1s[[id]]))
    pc1s[[dataset]] <<- list()
  if (is.null(pc1s[[id]][[as.character(bin_size)]]))
  {
    fname <- paste0(master.dir, "/data/rda/pc1_", dataset, "_", genome, "_", bin_size, ".rda")
    cat(paste0("Loading ", fname, "\n"))
    pc1 <- NULL
    load(fname)

    pc1s[[id]][[as.character(bin_size)]] <<- pc1
  }
}

ensure_loaded <- function(dataset, genome, bin_size, type = 'A')
{
  bin_size <- as.integer(bin_size)
  id <- paste0(dataset, "_", genome)

  if (!exists("map_params"))
    map_params <<- list()
  if (is.null(map_params[[id]]))
    map_params[[paste0(dataset, "_", genome)]] <<- list()

  if (!exists("chrom_maps"))
    chrom_maps <<- list()
  if (is.null(chrom_maps[[id]]))
    chrom_maps[[paste0(dataset, "_", genome)]] <<- list()

  if (!exists("map_bins"))
    map_bins <<- list()
  if (is.null(map_bins[[id]]))
    map_bins[[paste0(dataset, "_", genome)]] <<- list()

  if (!exists("maps"))
    maps <<- list()
  if (is.null(maps[[id]]))
    maps[[paste0(dataset, "_", genome)]] <<- list()

  if (is.null(maps[[id]][[as.character(bin_size)]]))
  {
    fname <- paste0(master.dir, "/data/hicexplorer/rda/", genome, "_", dataset, "_filtered_", bin_size, "_corrected.rda")
    if (!file.exists(fname))
      fname <- paste0(master.dir, "/data/rda/pdb_", dataset, "_", genome, "_", bin_size, ".rda")
    if (!file.exists(fname))
      fname <- paste0(master.dir, "/data/hicexplorer/rda/", genome, "_",
        sub("_All_combined$", "_combined_All", sub("_VRG_combined$", "_combined_VRG", sub("_VRGdownBAL_combined$", "_combined_VRGdownBAL",
          sub("_BAL_combined$", "_combined_BAL", sub("$", "_combined", dataset))))),
        "_filtered_", bin_size, "_corrected.rda")
    if (!file.exists(fname))
      fname <- paste0(master.dir, "/data/hicexplorer/rda/", genome, "_", dataset, "_filtered_", bin_size, "_corrected.rda")

    cat(paste0("Loading ", fname, "\n"))
    params <- NULL
    chrom_map <- NULL
    bins <- NULL
    map <- NULL
    load(fname)

    # chrom_map <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size, ".txt"), header = T, sep = "\t")
    # map$genome$expected <- matrix(trans_mean, nrow(map$genome$observed), ncol(map$genome$observed))

    # for (this_chrom in chrom_lst)
    # {
    #   d <- nrow(map[[this_chrom]]$observed)
    #   map[[this_chrom]]$expected <- matrix(NA, d, d)
    #   cis_mean <- stats_chrom$mean[stats_chrom$chrom == this_chrom]
    #   cis_dist <- stats_chrom$dist[stats_chrom$chrom == this_chrom] / bin_size

    #   for (i in seq_along(cis_mean))
    #   {
    #     dist <- cis_dist[i]
    #     stopifnot(dist < d)
    #     map[[this_chrom]]$expected[cbind(1:(d - dist), (1 + dist):d)] <- cis_mean[i]
    #     map[[this_chrom]]$expected[cbind((1 + dist):d, 1:(d - dist))] <- cis_mean[i]
    #   }

    #   m <- match(this_chrom, chrom_map$chrom)
    #   int <- chrom_map$map_start[m]:chrom_map$map_end[m]
    #   map$genome$expected[int, int] <- map[[this_chrom]]$expected
    # }

    map_params[[id]][[as.character(bin_size)]] <<- params
    chrom_maps[[id]][[as.character(bin_size)]] <<- chrom_map
    map_bins[[id]][[as.character(bin_size)]] <<- bins
    maps[[id]][[as.character(bin_size)]] <<- map
  }
}

ensure_loaded_distance_fit <- function(dataset, genome, bin_size)
{
  require(locfit, quietly = T)
  id <- paste0(dataset, "_", genome)
  if (!exists("distance_fits"))
    distance_fits <<- list()
  if (is.null(distance_fits[[id]]))
    distance_fits[[paste0(dataset, "_", genome)]] <<- list()
  if (is.null(distance_fits[[id]][[as.character(bin_size)]]))
  {
    fname <- paste0(master.dir, "/data/distance_decay/decay_", dataset, "_", genome, "_", max(bin_size, 5e3), ".rda")
    cat(paste0("Loading ", fname, "\n"))
    fit_norm <- NULL
    load(fname)

    distance_fits[[id]][[as.character(bin_size)]] <<- fit_norm$genome$locfit
    # very nasty trick to make locfit predict anything
    attr(distance_fits[[id]][[as.character(bin_size)]]$terms, ".Environment") <<- globalenv()
  }
}

ensure_loaded_int <- function(dataset, genome, bin_size)
{
  if (!exists("ints"))
  {
    ints <<- list()
    iuns <<- list()
  }
  if (is.null(ints[[dataset]]))
  {
    ints[[dataset]] <<- list()
    iuns[[dataset]] <<- list()
  }
  if (is.null(ints[[dataset]][[as.character(bin_size)]]))
  {
    fname <- paste0(master.dir, "/data/rda/int_", dataset, "_", genome, "_", bin_size, ".rda")
    cat(paste0("Loading ", fname, "\n"))
    int <- NULL
    load(fname)

    int <- int[int$p_value < 0.001]
    setkey(int, chrom, bin1, bin2)
    int <- unique(int)

    int$cluster <- seq_len(nrow(int))
    done <- F
    while (!done)
    {
      done <- T
      for (i in seq_len(nrow(int)))
      {
        sel <- int$chrom == int$chrom[i] & abs(int$bin1 - int$bin1[i]) <= 1 & abs(int$bin2 - int$bin2[i]) <= 1
        if (max(int$cluster[sel]) != min(int$cluster[sel]))
          done <- F
        int$cluster[sel] <- min(int$cluster[sel])
      }
    }

    # renumber clusters
    int$cluster <- match(int$cluster, sort(unique(int$cluster)))

    print(length(unique(int$cluster)))
    print(table(table(int$cluster)))

    # unique interactions
    iun <- int[, list(bin1 = bin1[which.min(p_value)], bin2 = bin2[which.min(p_value)],
      observed = observed[which.min(p_value)], expected = expected[which.min(p_value)], p_value = min(p_value),
      midbin1 = midbin1[which.min(p_value)], midbin2 = midbin2[which.min(p_value)]), by = c('chrom', 'cluster')]

    ints[[dataset]][[as.character(bin_size)]] <<- int
    iuns[[dataset]][[as.character(bin_size)]] <<- iun
  }
}

ensure_loaded_significant <- function(dataset, genome, bin_size)
{
  if (!exists("ints"))
  {
    ints <<- list()
    iuns <<- list()
  }
  if (is.null(ints[[dataset]]))
  {
    ints[[dataset]] <<- list()
    iuns[[dataset]] <<- list()
  }
  if (is.null(ints[[dataset]][[as.character(bin_size)]]))
  {
    fname <- paste0(master.dir, "/data/rda/significant_", dataset, "_", genome, "_", bin_size, ".rda")
    cat(paste0("Loading ", fname, "\n"))
    int <- NULL
    load(fname)

    ints[[dataset]][[as.character(bin_size)]] <<- int
    iuns[[dataset]][[as.character(bin_size)]] <<- iun
  }
}

fixed_levelplot <- function(matrix, midpointsx, midpointsy, main = NULL, xlab = NULL, colfun = colorRampPalette(c("white", "red")), at = NULL,
  panel.dim = dev.size("in")[1] - 2.3, legend.dim = 1.8, panel = panel.levelplot, ...)
{
  # ensure that the plot size is fixed
  if (is.null(main) || is.na(main))
    main <- " "
  if (is.null(xlab) || is.na(xlab))
    xlab <- " "

  if (is.null(at))
    at <- seq(from = min(matrix, na.rm = T), to = max(matrix, na.rm = T), length.out = 17)

  pl <- levelplot(matrix, main = main, xlab = xlab, ylab = NULL, colorkey = F, col.regions = colfun(length(at) - 1), at = at,
    row.values = midpointsx, column.values = rev(midpointsy), panel = panel, ...,
    xscale.components = function(...) { tmp.xscale <<- xscale.components.default(...); list(num.limit = tmp.xscale$num.limit) },
    yscale.components = function(...) { tmp.yscale <<- yscale.components.default(...); list(num.limit = tmp.yscale$num.limit) }
  )
  print(pl, panel.width = list(panel.dim, "inches"), panel.height = list(panel.dim, "inches"))

  key <- draw.colorkey(list(col = colfun(length(at) - 1), at = at))
  key$vp <- viewport(x = unit(0.5, "npc") + unit(panel.dim / 2, "inches") + unit(2, "char") + 0.5 * grobWidth(key), height = unit(legend.dim, "inches"))
  grid.draw(key)

  trellis.focus("panel", 1, 1, clip.off = TRUE, highlight = FALSE)
  yscaleat <- current.panel.limits()$ylim[2] - tmp.yscale$left$labels$at + current.panel.limits()$ylim[1] # need to flip because of rev(midpointsy) above
  panel.axis("left", yscaleat, labels = prettyNum(as.integer(tmp.yscale$left$labels$labels), big.mark = ","), outside = T)
  panel.axis("right", yscaleat, labels = F, outside = T)
  panel.axis("bottom", tmp.xscale$bottom$labels$at, labels = prettyNum(as.integer(tmp.xscale$bottom$labels$labels), big.mark = ","), rot = 0, outside = T)
  panel.axis("top", tmp.xscale$bottom$labels$at, labels = F, outside = T)
}

axis.fixed_triangle_levelplot <- function(orientation = "up", chrom = chrom, ...)
{
  if (chrom != "genome")
    grid.text(chrom, unit(0, "npc") - unit(0.5, "char"), switch(orientation, up = 0, down = 1), just = c(1, 0.5), gp = do.call(gpar, trellis.par.get("axis.text")))
  panel.axis(switch(orientation, up = "bottom", down = "top"),
    tmp.xscale$bottom$labels$at, labels = prettyNum(as.integer(tmp.xscale$bottom$labels$labels), big.mark = ","), rot = 0, outside = T)
}

fixed_triangle_levelplot <- function(matrix, midpoints, main = NULL, xlab = NULL, chrom = NA, colfun = colorRampPalette(c("white", "red")), at = NULL,
  panel.dim = dev.size("in")[1] - 1, legend.dim = 1.8, orientation = "up", panel = panel.levelplot, axis = axis.fixed_triangle_levelplot, scales = list(), xlim = xlim, ylim = ylim, ...)
{
  stopifnot(!is.null(switch(orientation, up =, down = T)))

  # ensure that the plot size is fixed
  if (is.null(main) || is.na(main))
    main <- " "
  if (is.null(xlab) || is.na(xlab))
    xlab <- " "

  if (is.null(at))
    at <- seq(from = min(matrix, na.rm = T), to = max(matrix, na.rm = T), length.out = 17)

  matrix[!upper.tri(matrix)] <- NA
  basey <- switch(orientation, up = unit(0, "npc") + unit(0.35, "inches"), down = unit(1, "npc") - unit(0.35, "inches"))

  grid.newpage()
  pushViewport(viewport(y = basey, angle = 45, clip = "off"))
  pl <- levelplot(switch(orientation, up = t(matrix), down = matrix), xlab = NULL, ylab = NULL, col.regions = colfun(length(at) - 1), at = at,
    row.values = midpoints, column.values = rev(midpoints), panel = panel, xlim = xlim, ylim = ylim, ..., midpoints = midpoints,
    colorkey = F, scales = list(x = scales, draw = F), par.settings = list(axis.line = list(col = NA)),
    xscale.components = function(...) { tmp.xscale <<- xscale.components.default(...); list(num.limit = tmp.xscale$num.limit) },
    yscale.components = function(...) { tmp.yscale <<- yscale.components.default(...); list(num.limit = tmp.yscale$num.limit) },
  )
  gp <- trellis.par.get("clip")
  trellis.par.set("clip", list(panel = "off"))
  print(pl, panel.width = list(panel.dim / sqrt(2), "inches"), panel.height = list(panel.dim / sqrt(2), "inches"), newpage = F)
  trellis.par.set("clip", gp)
  popViewport()

  grid.text(main, unit(0.1, "inches"), basey + unit(switch(orientation, up = 1, down = -1) * (panel.dim / 2 - 0.05), "inches"),
    just = c(0, switch(orientation, up = 1, down = 0)), gp = do.call(gpar, trellis.par.get("par.main.text")))
  grid.text(xlab, unit(0.1, "inches"), basey + unit(switch(orientation, up = 1, down = -1) * (panel.dim / 2 - 0.3), "inches"),
    just = c(0, switch(orientation, up = 1, down = 0)), gp = do.call(gpar, trellis.par.get("par.xlab.text")))

  key <- draw.colorkey(list(col = colfun(length(at) - 1), at = at))
  key$vp <- viewport(x = unit(panel.dim + 0.9, "inches") - 0.5 * grobWidth(key),
    y = basey + unit(switch(orientation, up = 1, down = -1) * (panel.dim / 2 - legend.dim / 2 - 0.05), "inches"), height = unit(legend.dim, "inches"))
  grid.draw(key)

  pushViewport(viewport(y = basey + unit(switch(orientation, up = 1, down = -1) * panel.dim / 2, "inches"), clip = "off"))
  pl <- levelplot(matrix, xlab = NULL, ylab = NULL, xlim = xlim, ylim = ylim,
    row.values = midpoints, column.values = rev(midpoints),
    panel = function(...) { switch(orientation, up = grid.segments(0, 0, 1, 0), down = grid.segments(0, 1, 1, 1)) },
    colorkey = F, scales = list(draw = F), par.settings = list(axis.line = list(col = NA)),
  )
  print(pl, panel.width = list(panel.dim, "inches"), panel.height = list(panel.dim, "inches"), newpage = F)

  trellis.focus("panel", 1, 1, clip.off = TRUE, highlight = FALSE)
  axis(orientation = orientation, chrom = chrom, ...)
  popViewport()
}

extract_mult <- function(dataset1, dataset2 = NA, genome, chrom, bin_size1, bin_size2 = bin_size1, ...)
{
  bin_size1 <- as.integer(bin_size1)
  bin_size2 <- as.integer(bin_size2)
  compartment <- identical(dataset2, "compartment")
  if (compartment)
    dataset2 <- dataset1

  finfo1 <- paste0(master.dir, "/data/hic_maps/", dataset1, "_", genome, "_", bin_size1, ".info")
  if (file.exists(finfo1))
    s1 <- read.table(finfo1)
  else
    s1 <- as.data.frame(matrix(nrow = 0, ncol = 2))

  finfo2 <- paste0(master.dir, "/data/hic_maps/", dataset2, "_", genome, "_", bin_size2, ".info")
  if (!is.na(dataset2) && file.exists(finfo2))
    s2 <- read.table(finfo2)
  else
    s2 <- as.data.frame(matrix(nrow = 0, ncol = 2))

  null.to.na <- function(x) if (is.null(x)) return(NA) else return(x)
  extract.int <- function(ss, key) strtoi(null.to.na(ss[match(key, ss[, 1]), 2]))

  id1 <- paste0(dataset1, "_", genome)
  id2 <- paste0(dataset2, "_", genome)
  ensure_loaded(dataset1, genome, bin_size1)
  if (!is.na(dataset2))
    ensure_loaded(dataset2, genome, bin_size2)

  if (!is.null(map_params[[id1]][[as.character(bin_size1)]]))
  {
    count1 <- map_params[[id1]][[as.character(bin_size1)]]$count
    count2 <- map_params[[id2]][[as.character(bin_size2)]]$count
    cv <- c(count1, count2)
    if (is.na(dataset2) || (abs(diff(cv)) / mean(cv) < 0.01 && bin_size1 == bin_size2))
    {
      mult1 <- count1 / map_params[[id1]][[as.character(bin_size1)]]$bins
      if (length(mult1) == 0 || !is.finite(mult1))
        mult1 <- 1.0
      mult2 <- count2 / map_params[[id2]][[as.character(bin_size2)]]$bins
      if (length(mult2) == 0 || !is.finite(mult2))
        mult2 <- 1.0
      # message("extract_mult (stored in map_params): ", mult1, ", ", mult2)
    }
    else
    {
      mult1 <- 1.0
      mult2 <- 1.0
    }
  }
  else if (is.na(dataset2) || (!is.na(extract.int(s1, "count_final_raw")) && extract.int(s1, "count_final_raw") == extract.int(s2, "count_final_raw") && bin_size1 == bin_size2))
  {
    mult1 <- extract.int(s1, "count_final_raw") / extract.int(s1, "bins")
    if (!is.finite(mult1))
      mult1 <- 1.0
    mult2 <- extract.int(s2, "count_final_raw") / extract.int(s2, "bins")
    if (!is.finite(mult2))
      mult2 <- 1.0
    # message("extract_mult (stored in .info): ", mult1, ", ", mult2)
  }
  else
  {
    mult1 <- 1.0
    mult2 <- 1.0
  }

  return(c(mult1, mult2))
}

map_distance_fit <- function(map, dataset, genome, bin_size, distance_fit)
{
  if (is.null(map) || is.null(dim(map)) || nrow(map) < 1)
    return(map)
  stopifnot(nrow(map) == ncol(map))
  d <- nrow(map)

  newdata <- data.table(distance = seq_len(d - 1L) * bin_size)
  newdata$value_dataset <- exp(predict(distance_fits[[paste0(dataset, "_", genome)]][[as.character(bin_size)]], newdata = newdata))
  newdata$value_distance_fit <- exp(predict(distance_fits[[paste0(distance_fit, "_", genome)]][[as.character(bin_size)]], newdata = newdata))
  newdata$ratio <- newdata$value_distance_fit / newdata$value_dataset

  map_distance_fit <- matrix(NA, d, d)
  for (i in seq_len(d - 1L))
    map_distance_fit[cbind(c(1:(d - i), (1 + i):d), c((1 + i):d, 1:(d - i)))] <- newdata$ratio[i]

  return(map_distance_fit * map)
}

extract_plot_data <- function(dataset1, dataset2, genome, chrom1, start1 = NA, end1 = NA, bin_size1, chrom2 = chrom1, start2 = start1, end2 = end1, bin_size2 = bin_size1, distance_fit = NA, mult1 = 1, mult2 = 1, what = "observed", filter = "none", zmin = NA, zmax = NA, ...)
{
  bin_size1 <- as.integer(bin_size1)
  bin_size2 <- as.integer(bin_size2)
  dataset1 <- gsub("_iterative", "", dataset1)
  dataset2 <- gsub("_iterative", "", dataset2)

  compartment <- identical(dataset2, "compartment")
  if (compartment)
  {
    dataset2 <- dataset1
    ensure_loaded_compartment(dataset1, genome, bin_size1)
  }

  if (what %in% c("p-value", "observed", "logobserved", "expected", "ratio", "minusLog10Padj", "log2FoldChange"))
  {
    ensure_loaded(dataset1, genome, bin_size1)
    if (!is.na(dataset2))
      ensure_loaded(dataset2, genome, bin_size2)
  }

  if (!is.na(distance_fit))
  {
    ensure_loaded_distance_fit(dataset1, genome, bin_size1)
    ensure_loaded_distance_fit(distance_fit, genome, bin_size1)
    if (!is.na(dataset2))
    {
      ensure_loaded_distance_fit(dataset2, genome, bin_size2)
      ensure_loaded_distance_fit(distance_fit, genome, bin_size2)
    }
  }

  if (chrom1 == "genome")
  {
    chrom1_lst <- chrom_lst
    if (grepl("HiC_DB_", dataset1) || grepl("HiC_DB_", dataset2))
      chrom1_lst <- chrom_lst[1:4]
    print(chrom1_lst)
    coords <- paste(chrom1_lst, collapse = ", ")
    fun_interval <- function(chrom) c((chrom_map$map_start[match(chrom, chrom_map$chrom)]):(chrom_map$map_end[match(chrom, chrom_map$chrom)]),
      rep(NA, ifelse(chrom == tail(chrom1_lst, 1), 0, chrom_sep_bins))) # with chrom_sep_bins bin margin between chromosomes

    chrom_map <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size1, ".txt"), header = T, sep = "\t")
    bins1x <- as.vector(unlist(sapply(chrom1_lst, fun_interval)))
    bins1y <- bins1x

    chrom_map <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size2, ".txt"), header = T, sep = "\t")
    bins2x <- as.vector(unlist(sapply(chrom1_lst, fun_interval)))
    bins2y <- bins2x

    midpointsx <- (seq_along(bins1x) - 0.5) * bin_size1
    midpointsy <- (seq_along(bins1y) - 0.5) * bin_size2
    xlim <- range(0, length(bins1x) * bin_size1, length(bins2x) * bin_size2)
    ylim <- range(0, length(bins1y) * bin_size1, length(bins2y) * bin_size2)
  }
  else
  {
    chrom_map1 <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size1, ".txt"), header = T, sep = "\t")
    if (is.na(start1))
      start1 <- 1
    if (is.na(end1))
      end1 <- chrom_map1$length[match(chrom1, chrom_map1$chrom)]

    chrom_map2 <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size2, ".txt"), header = T, sep = "\t")
    if (is.na(start2))
      start2 <- 1
    if (is.na(end2))
      end2 <- chrom_map2$length[match(chrom2, chrom_map2$chrom)]

    coords <- paste0(chrom1, ":", prettyNum(as.integer(start1), big.mark = ","), "-", prettyNum(as.integer(end1), big.mark = ","))
    if (chrom1 != chrom2 || !identical(start1, start2) || !identical(end1, end2))
      coords <- paste0(chrom2, ":", prettyNum(as.integer(start2), big.mark = ","), "-", prettyNum(as.integer(end2), big.mark = ","), " / ", coords)

    bins1x <- as.integer((start1 - 1L) %/% bin_size1 + 1):as.integer((end1 - 1L) %/% bin_size1 + 1)
    bins1y <- as.integer((start2 - 1L) %/% bin_size1 + 1):as.integer((end2 - 1L) %/% bin_size1 + 1)
    bins2x <- as.integer((start1 - 1L) %/% bin_size2 + 1):as.integer((end1 - 1L) %/% bin_size2 + 1)
    bins2y <- as.integer((start2 - 1L) %/% bin_size2 + 1):as.integer((end2 - 1L) %/% bin_size2 + 1)
    midpointsx <- (bins1x - 0.5) * bin_size1
    midpointsy <- (bins1y - 0.5) * bin_size1
    xlim <- range((bins1x - 1) * bin_size1, bins1x * bin_size1, (bins2x - 1) * bin_size2, bins2x * bin_size2)
    ylim <- range((bins1y - 1) * bin_size1, bins1y * bin_size1, (bins2y - 1) * bin_size2, bins2y * bin_size2)

    if (chrom1 != chrom2)
    {
      bins1x <- bins1x + chrom_map1$map_start[match(chrom1, chrom_map1$chrom)]
      bins1y <- bins1y + chrom_map2$map_start[match(chrom2, chrom_map2$chrom)]
      bins2x <- bins2x + chrom_map1$map_start[match(chrom1, chrom_map1$chrom)]
      bins2y <- bins2y + chrom_map2$map_start[match(chrom2, chrom_map2$chrom)]
      chrom1 <- "genome"
      chrom2 <- "genome"
    }
  }

  if (chrom1 != "genome")
  {
    dim1 <- chrom_map1$map_end[match(chrom1, chrom_map1$chrom)] - chrom_map1$map_start[match(chrom1, chrom_map1$chrom)] + 1L
    if (min(bins1x) < 1 || dim1 < max(bins1x)) warning(paste0("bin coordinates outside range [1, ", dim1, "]"))
    bins1x <- ifelse(bins1x < 1 | dim1 < bins1x, NA, bins1x)
    if (min(bins1y) < 1 || dim1 < max(bins1y)) warning(paste0("bin coordinates outside range [1, ", dim1, "]"))
    bins1y <- ifelse(bins1y < 1 | dim1 < bins1y, NA, bins1y)

    dim2 <- chrom_map2$map_end[match(chrom2, chrom_map2$chrom)] - chrom_map2$map_start[match(chrom2, chrom_map2$chrom)] + 1L
    if (!is.na(dataset2) && (min(bins2x) < 1 || dim2 < max(bins2x)))   warning(paste0("bin coordinates outside range [1, ", dim2, "]"))
    bins2x <- ifelse(bins2x < 1 | dim2 < bins2x, NA, bins2x)
    if (!is.na(dataset2) && (min(bins2y) < 1 || dim2 < max(bins2y)))   warning(paste0("bin coordinates outside range [1, ", dim2, "]"))
    bins2y <- ifelse(bins2y < 1 | dim2 < bins2y, NA, bins2y)
  }

  if (what == "p-value")
  {
    map1 <- -log10(maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$p_value[bins1x, bins1y])
    map2 <- -log10(maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$p_value[bins2x, bins2y])
  }
  else if (what == "observed")
  {
    if (chrom1 == "genome" & !"observed" %in% names(maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]))
      stop("no trans-chromosome contacts saved in .rda file")
    map1 <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$observed[bins1x, bins1y] * mult1

    if (!is.na(dataset2) & chrom2 == "genome" & !"observed" %in% names(maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]))
      stop("no trans-chromosome contacts saved in .rda file")
    map2 <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$observed[bins2x, bins2y] * mult2

    if (!is.na(distance_fit))
    {
      map1 <- map_distance_fit(map1, dataset1, genome, bin_size1, distance_fit)
      map2 <- map_distance_fit(map2, dataset2, genome, bin_size2, distance_fit)
    }
  }
  else if (what == "logobserved")
  {
    map1 <- log10(maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$observed[bins1x, bins1y] * mult1)
    map2 <- log10(maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$observed[bins2x, bins2y] * mult2)
  }
  else if (what == "expected")
  {
    map1 <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$expected[bins1x, bins1y] * mult1
    map2 <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$expected[bins2x, bins2y] * mult2
  }
  else if (what == "ratio")
  {
    map1 <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$observed[bins1x, bins1y] / maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$expected[bins1x, bins1y]
    map2 <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$observed[bins2x, bins2y] / maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$expected[bins2x, bins2y]
  }
  else if (what == "minusLog10Padj")
  {
    map1 <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$minusLog10Padj[bins1x, bins1y]
    map2 <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$minusLog10Padj[bins2x, bins2y]
  }
  else if (what == "log2FoldChange")
  {
    map1 <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$log2FoldChange[bins1x, bins1y]
    map2 <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$log2FoldChange[bins2x, bins2y]
  }
  else
    stop(paste0('unknown value for "what": ', what))

  if (what %in% c("p-value", "observed", "logobserved", "expected", "ratio"))
  {
    if (filter == "onesignificant")
    {
      map1_pval <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$p_value[bins1x, bins1y]
      map2_pval <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$p_value[bins2x, bins2y]
      map1 <- ifelse(map1_pval < 0.01 | map2_pval < 0.01, map1, NA)
      map2 <- ifelse(map1_pval < 0.01 | map2_pval < 0.01, map2, NA)
    }
    else if (filter == "bothsignificant")
    {
      map1_pval <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$p_value[bins1x, bins1y]
      map2_pval <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$p_value[bins2x, bins2y]
      map1 <- ifelse(map1_pval < 0.01, map1, NA)
      map2 <- ifelse(map2_pval < 0.01, map2, NA)
    }
    else if (filter != "none")
      stop(paste0('unknown value for "filter": ', filter))
  }
  if (what %in% c("minusLog10Padj", "log2FoldChange"))
  {
    if (filter == "significant")
    {
      map1_pval <- maps[[paste0(dataset1, "_", genome)]][[as.character(bin_size1)]][[chrom1]]$minusLog10Padj[bins1x, bins1y]
      map2_pval <- maps[[paste0(dataset2, "_", genome)]][[as.character(bin_size2)]][[chrom2]]$minusLog10Padj[bins2x, bins2y]
      map1 <- ifelse(map1_pval > 1, map1, 0)
      map2 <- ifelse(map2_pval > 1, map2, 0)
    }
    else if (filter != "none")
      stop(paste0('unknown value for "filter": ', filter))
  }

  if (compartment)
  {
    vpc <- pc1s[[dataset1]][[as.character(bin_size1)]][bins1x]

    map2 <- matrix(0, length(vpc), length(vpc))
    d <- length(vpc)
    i <- 1
    map2[cbind(c(1:(d - i), (1 + i):d), c((1 + i):d, 1:(d - i)))] <- head(tail(vpc, -1), -1) * -10
    map2[!upper.tri(map2)] <- NA
  }

  if (!is.na(zmax))
  {
    map1 <- pmin(map1, zmax - 1e-6)
    map2 <- pmin(map2, zmax - 1e-6)
  }
  if (!is.na(zmin))
  {
    map1 <- pmax(map1, zmin + 1e-6)
    map2 <- pmax(map2, zmin + 1e-6)
  }
  if (chrom1 == chrom2 && identical(start1, start2) && identical(end1, end2))
  {
    map1[!upper.tri(map1)] <- NA
    map2[!upper.tri(map2)] <- NA
  }

  xlab <- paste0(coords, paste0(" [", genome, "]; resolution ", ifelse(bin_size1 == bin_size2 || is.na(dataset2), bin_size1 / 1e3, paste0(bin_size1 / 1e3, " kb \\ ", bin_size2 / 1e3)), " kb"))

  return(list(map1 = map1, map2 = map2, mult1 = mult1, mult2 = mult2, midpointsx = midpointsx, midpointsy = midpointsy, xlim = xlim, ylim = ylim, xlab = xlab,
    bins1x = bins1x, bins1y = bins1y, bins2x = bins2x, bins2y = bins2y))
}

square_plot_mark_loops <- function(dataset1, dataset2, genome, this_chrom, start, end, bin_size1, bin_size2, radius = 1.5, col = "magenta", ...)
{
  r <- unit(radius * bin_size1, "native") + unit(72.27 / 96 / 2, "points")
  if (!is.null(loops[[dataset1]]))
  {
    l <- loops[[dataset1]]
    l <- l[l$chrom == this_chrom & l$position1 >= start - 2 * bin_size1 & l$position1 <= end + 2 * bin_size1
      & l$position2 >= start - 2 * bin_size1 & l$position2 <= end + 2 * bin_size1, ]
    if (nrow(l) > 0)
    {
      if (is.null(l$col))
        l$col <- col
      l$col[is.na(l$col)] <- col
      grid.circle(l$position1, current.panel.limits()$ylim[2] - l$position2 + current.panel.limits()$ylim[1], r, default.units = "native", gp = gpar(col = l$col, fill = NA))
      grid.text(l$id, l$position1, current.panel.limits()$ylim[2] - (l$position2 - radius * bin_size1) + current.panel.limits()$ylim[1], vjust = -0.5, default.units = "native", gp = gpar(col = l$col, cex = 0.8))
    }
  }
  if (!is.null(loops[[dataset2]]))
  {
      if (is.null(l$col))
        l$col <- col
      l$col[is.na(l$col)] <- col
    l <- loops[[dataset2]]
    l <- l[l$chrom == this_chrom & l$position1 >= start - 2 * bin_size2 & l$position1 <= end + 2 * bin_size2
      & l$position2 >= start - 2 * bin_size2 & l$position2 <= end + 2 * bin_size2, ]
    if (nrow(l) > 0)
    {
      grid.circle(l$position2, current.panel.limits()$ylim[2] - l$position1 + current.panel.limits()$ylim[1], r, default.units = "native", gp = gpar(col = l$col, fill = NA))
      grid.text(l$id, l$position2, current.panel.limits()$ylim[2] - (l$position1 + radius * bin_size1) + current.panel.limits()$ylim[1], vjust = 1.5, default.units = "native", gp = gpar(col = l$col, cex = 0.8))
    }
  }
}

triangle_plot_mark_loops <- function(dataset1, genome, this_chrom, start, end, bin_size1, radius = 1.5, col = "magenta", ...)
{
  r <- unit(radius * bin_size1, "native") + unit(72.27 / 96 / 2, "points")
  if (!is.null(loops[[dataset1]]))
  {
    l <- loops[[dataset1]]
    l <- l[l$chrom == this_chrom & l$position1 >= start - 2 * bin_size1 & l$position1 <= end + 2 * bin_size1
      & l$position2 >= start - 2 * bin_size1 & l$position2 <= end + 2 * bin_size1, ]
    if (nrow(l) > 0)
    {
      # print(l[, 1:7])
      if (is.null(l$col))
        l$col <- col
      l$col[is.na(l$col)] <- col
      grid.circle(l$position2, current.panel.limits()$ylim[2] - l$position1 + current.panel.limits()$ylim[1], r, default.units = "native", gp = gpar(col = l$col, fill = NA))
      grid.text(l$id, l$position2, current.panel.limits()$ylim[2] - (l$position1 + radius * bin_size1) + current.panel.limits()$ylim[1], vjust = -0.5, default.units = "native", gp = gpar(col = l$col, cex = 0.8))
    }
  }
}

extract_breakpoints <- function(genome, this_chrom, start = NA, end = NA, bin_size, summarize = NULL, range_thr = NA, ...)
{
  br <- fread(paste0(balancer.dir, "/analysis/breakpoints_", genome, ".tab"), header = T, sep = "\t")
  if (is.null(br$ref.coords))
    br$ref.coords <- NA

  if (this_chrom == "genome")
  {
    chrom_map <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size, ".txt"), header = T, sep = "\t")
    br$breakpoint <- br$breakpoint + (chrom_map$map_start[match(br$chrom, chrom_map$chrom)] - 1 +
      (match(br$chrom, chrom_map$chrom) - 1) * chrom_sep_bins) * bin_size
  }
  else
    br <- br[br$chrom == this_chrom & ifelse(is.na(start), -Inf, start) <= br$breakpoint & br$breakpoint <= ifelse(is.na(end), Inf, end), ]

  find_match <- function(breakpoint, ref.coords)
  {
    if (is.na(match(summarize, ref.coords)))
    {
      if (length(breakpoint) > 1)
        warning("more than one breakpoint having same id, but none of them has ref.coords == ", summarize)
      return(breakpoint)
    }
    else
      return(breakpoint[summarize == ref.coords])
  }

  mean_if_range_lt <- function(v, range_thr)
  {
    r <- diff(range(v))
    if (r < range_thr)
      return(mean(v))
    else 
      return(as.double(v))
  }

  if (is.null(summarize))
    NULL # if summarize == NULL, do not change anything
  else if (summarize == "mean")
    br <- br[, list(breakpoint = mean(breakpoint)), by = "id"]
  else if (summarize == "mean_if_range_lt")
    br <- br[, list(breakpoint = mean_if_range_lt(breakpoint, range_thr)), by = "id"]
  else
    br <- br[, list(breakpoint = find_match(breakpoint, ref.coords)), by = "id"]

  # print(br$breakpoint)
  return(unique(br$breakpoint))
}

square_plot_mark_breakpoints <- function(genome, chrom1, start1, end1, bin_size1, chrom2 = chrom1, start2 = start1, end2 = end1, bin_size2 = bin_size1, lwd.breakpoints = 1, lty.breakpoints = "dashed", lineend.breakpoints = "butt", col.breakpoints = "magenta", ...)
{
  br1 <- extract_breakpoints(genome, chrom1, start1, end1, bin_size1, summarize = breakpoint_summarize)
  br2 <- extract_breakpoints(genome, chrom2, start2, end2, bin_size2, summarize = breakpoint_summarize)
  if (length(br2) > 0)
    panel.abline(h = current.panel.limits()$ylim[2] - br2 + current.panel.limits()$ylim[1],
    lwd = lwd.breakpoints, lty = lty.breakpoints, lineend = lineend.breakpoints, col = col.breakpoints)
  if (length(br1) > 0)
    panel.abline(v = br1, lwd = lwd.breakpoints, lty = lty.breakpoints, lineend = lineend.breakpoints, col = col.breakpoints)
}

triangle_plot_mark_breakpoints <- function(lwd.breakpoints = 1, lty.breakpoints = "dashed", lineend.breakpoints = "butt", col.breakpoints = "magenta", ...)
{
  br_min <- extract_breakpoints(summarize = "lower", ...)
  br_max <- extract_breakpoints(summarize = "higher", ...)
  if (length(br_min) > 0)
  {
    panel.segments(unit(br_max, "native"), unit(current.panel.limits()$ylim[2] - br_max + current.panel.limits()$ylim[1], "native"),
      unit(1, "npc"), unit(current.panel.limits()$ylim[2] - br_max + current.panel.limits()$ylim[1], "native"),
      lwd = lwd.breakpoints, lty = lty.breakpoints, lineend = lineend.breakpoints, col = col.breakpoints)
    panel.segments(br_min, unit(current.panel.limits()$ylim[2] - br_min + current.panel.limits()$ylim[1], "native"), br_min, unit(1, "npc"),
      lwd = lwd.breakpoints, lty = lty.breakpoints, lineend = lineend.breakpoints, col = col.breakpoints)
  }
}

extract_domains <- function(dataset, domain_source, genome, this_chrom, start = NA, end = NA, bin_size, ...)
{
  if (identical(domain_source, T))
    domain_source <- "5000_IS100000"

  chrom_map <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size, ".txt"), header = T, sep = "\t")
  d <- fread(paste0(master.dir, "/data/domains/", dataset, "_", genome, "_", domain_source, ".bed"), sep = "\t")
  names(d) <- c("chrom", "start", "end", "name", "score", "strand")
  d$midpoint <- (d$start + d$end) / 2

  d <- d[, list(midpoint = c(0, midpoint, chrom_map$length[match(chrom, chrom_map$chrom)])), by = chrom]
  # domains <- d[, list(start = c(0, midpoint), end = c(midpoint, chrom_map$length[match(chrom, chrom_map$chrom)])), by = chrom]

  if (this_chrom == "genome")
  {
    boundaries <- (chrom_map$map_start[match(d$chrom, chrom_map$chrom)] - 1 +
      (match(br$chrom, chrom_map$chrom) - 1) * chrom_sep_bins) * bin_size + d$midpoint
  }
  else
  {
    sel <- d$chrom == this_chrom & ifelse(is.na(start), -Inf, start) <= d$midpoint & d$midpoint <= ifelse(is.na(end), Inf, end)
    boundaries <- d$midpoint[c(F, head(sel, -1)) | sel | c(tail(sel, -1), F)]
  }

  print(boundaries)
  return(boundaries)
}

square_plot_mark_domains <- function(dataset1, dataset2, domain_source, genome, chrom1, start1, end1, bin_size1, chrom2 = chrom1, start2 = start1, end2 = end1, bin_size2 = bin_size1, ...)
{
  d1 <- extract_domains(dataset1, domain_source, genome, chrom1, start1, end1, bin_size1)
  d1x <- head(rep(d1, each = 2), -1)
  d1y <- tail(rep(d1, each = 2), -1)

  d2 <- extract_domains(dataset2, domain_source, genome, chrom2, start2, end2, bin_size2)
  d2x <- tail(rep(d2, each = 2), -1)
  d2y <- head(rep(d2, each = 2), -1)

  trellis.focus("panel", 1, 1, clip.off = FALSE, highlight = FALSE)
  if (length(d1) > 0)
    panel.lines(x = d1x, y = current.panel.limits()$ylim[2] - d1y + current.panel.limits()$ylim[1], col = "royalblue")
  if (length(d2) > 0)
    panel.lines(x = d2x, y = current.panel.limits()$ylim[2] - d2y + current.panel.limits()$ylim[1], col = "royalblue")
  trellis.focus("panel", 1, 1, clip.off = TRUE, highlight = FALSE)
}

triangle_plot_mark_domains <- function(dataset1, domain_source, genome, this_chrom, start, end, bin_size1, ...)
{
  stopifnot(current.panel.limits()$xlim == current.panel.limits()$ylim)
  lim <- current.panel.limits()$xlim

  d1 <- extract_domains(dataset1, domain_source, genome, this_chrom, start, end, bin_size1)
  d1x <- pmin(pmax(tail(rep(d1, each = 2), -1), lim[1]), lim[2])
  d1y <- pmin(pmax(head(rep(d1, each = 2), -1), lim[1]), lim[2])

  if (length(d1) > 1)
  {
    if (head(d1, 1) < lim[1])
    {
      d1x <- tail(d1x, -1)
      d1y <- tail(d1y, -1)
    }

    if (tail(d1, 1) > lim[2])
    {
      d1x <- head(d1x, -1)
      d1y <- head(d1y, -1)
    }

    panel.lines(x = d1x, y = current.panel.limits()$ylim[2] - d1y + current.panel.limits()$ylim[1], col = "royalblue")
  }
}

supplement_limits <- function(dataset1, dataset2, bin_size1, bin_size2 = bin_size1, distance_fit = NA, mult1 = 1, mult2 = 1, what, filter = "none", zmin = NA, zmax = NA, ...)
{
  if (what == "p-value")
  {
    if (is.na(zmax))
    {
      zmax <- 4
      if (max(bin_size1, bin_size2) < 5000) zmax <- 3
      if (max(bin_size1, bin_size2) > 5000) zmax <- 8
    }
    at <- c(0, 5:(5 * zmax) / 5)
    suffix <- ", interaction -log10 p-value"
  }
  else if (what == "observed")
  {
    zmin <- 0
    if (is.na(zmax))
      zmax <- 2e-6 * max(bin_size1, bin_size2)
    zmax <- zmax * max(mult1, mult2)
    # if (zmax >= 1000)
    #   at <- 0:ceiling(zmax / 10) * 10
    # else
    at <- 0:50 / 50 * zmax
    suffix <- paste0(", observed Hi-C reads", ifelse(is.na(distance_fit), "", paste0(", distance fit from ", distance_fit)))
  }
  else if (what == "logobserved")
  {
    stopifnot(is.finite(zmin))
    # if (max(map2) <= 1)
    #   map2 <- map2 / median(map2, na.rm = T) * median(map1, na.rm = T)
    if (is.na(zmin))
      zmin <- log10(2e-6 * max(bin_size1, bin_size2))
    zmin <- zmin + log10(min(mult1, mult2))
    # if (is.na(zmax))
    #   zmax <- log10(9e-7 * max(bin_size1, bin_size2)^2)
    if (is.na(zmax))
      zmax <- -2.5
    zmax <- zmax + log10(max(mult1, mult2))
    if (zmax >= 2)
      at <- 0:ceiling(zmax * 10) / 5
    else
      at <- 0:50 / 50 * (zmax - zmin) + zmin
    suffix <- ", log10(observed)"
  }
  else if (what == "expected")
  {
    if (is.na(zmax))
      zmax <- 2e-6 * max(bin_size1, bin_size2)
    zmax <- zmax * max(mult1, mult2)
    if (zmax >= 1000)
      at <- 0:ceiling(zmax / 10) * 10
    else
      at <- 0:50 / 50 * zmax
    suffix <- ", expected Hi-C reads"
  }
  else if (what == "ratio")
  {
    if (is.na(zmax))
    {
      zmax <- 5
      if (max(bin_size1, bin_size2) < 5000) zmax <- 10
    }
    at <- c(0, 5:(5 * zmax) / 5)
    suffix <- ", observed/expected ratio"
  }
  else if (what == "minusLog10Padj")
  {
    if (is.na(zmax))
    {
      zmax <- 4
      if (max(bin_size1, bin_size2) < 5000) zmax <- 3
      if (max(bin_size1, bin_size2) > 5000) zmax <- 8
    }
    if (zmax >= 100)
      at <- 0:50 / 50 * zmax
    else
      at <- c(0, 5:(5 * zmax) / 5)
    suffix <- ", DESeq2 -log10 adjusted p-value"
  }
  else if (what == "log2FoldChange")
  {
    if (is.na(zmax))
      zmax <- 3
    zmax <- zmax * max(mult1, mult2)
    zmin <- -zmax
    at <- -25:25 / 25 * zmax
    suffix <- ", log2FoldChange"
  }
  else
    stop(paste0('unknown value for "what": ', what))

  if (filter == "onesignificant")
    suffix <- paste0(suffix, " filtered by p-value")
  else if (filter == "bothsignificant")
    suffix <- paste0(suffix, " filtered by p-value")
  else if (filter == "significant")
    suffix <- paste0(suffix, " [p<0.1]")

  if (identical(dataset2, "compartment"))
    suffix <- paste0("alization", suffix) # compartment...

  return(list(zmin = zmin, zmax = zmax, at = at, suffix = suffix))
}

colscale_to_colvec <- function(colscale)
{
  if (colscale == "orangered")
    return(c("white", "orange", "darkred", "black"))
  else if (colscale == "red")
    return(c("white", "red"))
  else if (colscale == "OrPu")
    return(c("#e66101", "white", "#5e3c99"))
  else if (colscale == "bluered")
    return(c("#336699", "white", "#cc3333"))
  stop(paste0('unknown value for "colscale": ', colscale))
}

colscale_to_colfun <- function(colscale)
{
  return(colorRampPalette(colscale_to_colvec(colscale)))
}

square_plot <- function(dataset1, dataset2, genome, chrom1, start1 = NA, end1 = NA, bin_size1, chrom2 = chrom1, start2 = start1, end2 = end1, bin_size2 = bin_size1, distance_fit = NA, what = "observed", filter = "none", raster = NA, zmin = NA, zmax = NA, secondary.zmax = NA, colscale = "orangered", mark_loops = F, mark_breakpoints = F, mark_domains = F, scale = "(do not propagate 'scale' to levelplot)", ...)
{
  mult <- extract_mult(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2)
  data <- supplement_limits(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2,
    distance_fit = distance_fit, mult1 = mult[1], mult2 = mult[2], what = what, filter = filter, zmin = zmin, zmax = zmax)
  zmin <- data$zmin
  zmax <- data$zmax
  at <- data$at
  suffix <- data$suffix

  data <- extract_plot_data(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2,
    distance_fit = distance_fit, mult1 = mult[1], mult2 = mult[2], chrom1 = chrom1, start1 = start1, end1 = end1,
    chrom2 = chrom2, start2 = start2, end2 = end2, what = what, filter = filter, zmin = zmin, zmax = zmax)
  map1 <- data$map1
  map2 <- data$map2
  midpointsx <- data$midpointsx
  midpointsy <- data$midpointsy
  xlim <- data$xlim
  ylim <- data$ylim
  xlab <- data$xlab

  if (is.na(raster))
    raster <- max(dim(map1), dim(map2)) > 202

  if (bin_size1 == bin_size2)
    panel.fun <- mypanel.levelplot.same_bin_sizes
  else
    panel.fun <- mypanel.levelplot.different_bin_sizes
  if (identical(dataset2, "compartment"))
    panel.fun <- mypanel.levelplot.compartment

  colfun <- colscale_to_colfun(colscale)

  fixed_levelplot(map1, midpointsx, midpointsy, main = ifelse(is.na(dataset2), paste0(dataset1, suffix), paste0(dataset1, " \\ ", dataset2, suffix)), xlab = xlab, at = at, colfun = colfun,
    secondary = map2, secondary.midpoints = midpointsy, secondary.zmax = secondary.zmax, xlim = xlim, ylim = ylim,
    panel = panel.fun, panel.internal = ifelse(raster, panel.levelplot.raster, panel.levelplot), ...)

  if (mark_loops)
    square_plot_mark_loops(dataset1, dataset2, genome, chrom1, start1, end1, bin_size1, bin_size2, ...)
  if (mark_breakpoints)
    square_plot_mark_breakpoints(genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2, chrom1 = chrom1, start1 = start1, end1 = end1, chrom2 = chrom2, start2 = start2, end2 = end2, ...)
  if (identical(mark_domains, T) || is.character(mark_domains))
    square_plot_mark_domains(dataset1 = dataset1, dataset2 = dataset2, domain_source = mark_domains, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2, chrom1 = chrom1, start1 = start1, end1 = end1, chrom2 = chrom2, start2 = start2, end2 = end2, ...)
}

square_balancer_plot <- function(dataset1, dataset2, genome, this_chrom, start = NA, end = NA, bin_size1, bin_size2 = bin_size1, distance_fit = NA, what = "observed", filter = "none", raster = NA, zmin = NA, zmax = NA, secondary.zmax = NA, colscale = "orangered", mark_loops = F, mark_breakpoints = F, scale = "(do not propagate 'scale' to levelplot)", ...)
{
  stopifnot(dataset1 == "HiC_DB_6-8h_VRG")
  stopifnot(dataset2 == "HiC_DB_6-8h_BAL")
  stopifnot(genome == "dm6")

  read_breakpoints <- function(fname)
  {
    bp <- fread(fname, sep = "\t")
    names(bp) <- c("chrom", "start", "end", "name", "score", "strand")
    bp$pos <- ifelse(xor(grepl("_A$", bp$name), bp$strand == "+"), bp$start, bp$end)
    return(bp)
  }

  bp1 <- read_breakpoints(paste0(balancer.dir, "/analysis/breakpoints_dm6.bed"))
  bp2 <- read_breakpoints(paste0(balancer.dir, "/analysis/breakpoints_dm6bal3.bed"))
  stopifnot(bp1$name == bp2$name)

  for (i in seq_len(nrow(bp1)))
  {
    chrom1 <- bp1$chrom[i]
    pos1 <- bp1$pos[i]
    chrom2 <- bp2$chrom[i]
    pos2 <- bp2$pos[i]
    stopifnot(bp1$strand[i] == "+")

    # if the requested interval does not overlap the current row in bp1, skip this value of i
    if (chrom1 != this_chrom || max(start - 1, bp1$start[i]) > min(end, bp1$end[i]))
      next
    print(cbind(bp1, bp2)[i])

    bal_to_vrg <- function(x)
    {
      if (bp2$strand[i] == "-")
        y <- rev(pos2 - (x - pos2))
      else
        y <- x
      return(y - pos2 + pos1)
    }
    vrg_to_bal <- function(x)
    {
      if (bp2$strand[i] == "-")
        y <- rev(pos1 - (x - pos1))
      else
        y <- x
      return(y - pos1 + pos2)
    }
    xylim1 <- c(start - 1, end)
    xylim2 <- sort(vrg_to_bal(xylim1))
    start2 <- xylim2[1] + 1
    end2 <- xylim2[2]
    # print(c(chrom2, start2, end2))

    mult <- extract_mult(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2)
    data <- supplement_limits(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2,
      distance_fit = distance_fit, mult1 = mult[1], mult2 = mult[2], what = what, filter = filter, zmin = zmin, zmax = zmax)
    # mult <- c(1, 1)
    # data <- supplement_limits(dataset1 = dataset1, dataset2 = NA, genome = "dm6", bin_size1 = bin_size1,
    #   distance_fit = distance_fit, mult1 = mult[1], mult2 = mult[2], what = "observed", zmax = 0.01)
    zmin <- data$zmin
    zmax <- data$zmax
    at <- data$at
    suffix <- data$suffix

    data <- extract_plot_data(dataset1 = dataset1, dataset2 = NA, genome = "dm6", bin_size1 = bin_size1, distance_fit = distance_fit,
      chrom1 = chrom1, start1 = start, end1 = end, zmin = zmin, zmax = zmax)
    map1 <- data$map1
    midpointsx <- data$midpointsx
    xlab1 <- data$xlab

    data2 <- extract_plot_data(dataset1 = dataset2, dataset2 = NA, genome = "dm6bal3", bin_size1 = bin_size1, distance_fit = distance_fit,
      chrom1 = chrom2, start1 = start2, end1 = end2, zmin = zmin, zmax = zmax)
    map2 <- data2$map1
    midpointsy <- data2$midpointsx
    xlab2 <- data2$xlab
    # take the reverse complement if necessary
    if (bp2$strand[i] == "-")
      map2 <- t(map2)[nrow(map2):1, ][, ncol(map2):1]

    if (is.na(raster))
      raster <- max(dim(map1), dim(map2)) > 202
    panel.fun <- mypanel.levelplot.different_bin_sizes
    colfun <- colscale_to_colfun(colscale)

    xlab <- paste0(strsplit(xlab1, ";")[[1]][1], " \\ ", strsplit(xlab2, ";")[[1]][1],
      ifelse(bp2$strand[i] == "-", " revcomp", ""), ";", strsplit(xlab2, ";")[[1]][2])
    fixed_levelplot(map1, midpointsx, midpointsy, main = paste0(dataset1, " \\ ", dataset2, suffix), xlab = xlab, at = at, colfun = colfun,
      secondary = map2, secondary.midpoints = bal_to_vrg(midpointsy), xlim = xylim1, ylim = xylim1,
      panel = panel.fun, panel.internal = ifelse(raster, panel.levelplot.raster, panel.levelplot), ...)

    # if (mark_loops)
    #   square_plot_mark_loops(dataset1, dataset2, genome, this_chrom, start, end, bin_size1, bin_size2, ...)

    if (mark_breakpoints)
    {
      # square_plot_mark_breakpoints(genome, this_chrom, start, end, bin_size1, bin_size2, ...)

      br1 <- extract_breakpoints("dm6", chrom1, start, end = end, bin_size1)
      panel.segments(unit(br1, "native"), unit(xylim1[2] - br1 + xylim1[1], "native"),
        unit(0, "npc"), unit(xylim1[2] - br1 + xylim1[1], "native"), lwd = 2, lty = "dashed", lineend = "butt", col = "magenta")
      panel.segments(br1, unit(xylim1[2] - br1 + xylim1[1], "native"), br1, unit(0, "npc"), lwd = 2, lty = "dashed", lineend = "butt", col = "magenta")

      br2 <- extract_breakpoints("dm6bal3", chrom2, start2, end = end2, bin_size1)
      panel.segments(unit(bal_to_vrg(br2), "native"), unit(xylim1[2] - bal_to_vrg(br2) + xylim1[1], "native"),
        unit(1, "npc"), unit(xylim1[2] - bal_to_vrg(br2) + xylim1[1], "native"), lwd = 2, lty = "dashed", lineend = "butt", col = "magenta")
      panel.segments(bal_to_vrg(br2), unit(xylim1[2] - bal_to_vrg(br2) + xylim1[1], "native"),
        bal_to_vrg(br2), unit(1, "npc"), lwd = 2, lty = "dashed", lineend = "butt", col = "magenta")

      arr_coords <- function(name, strand, pos)
      {
        if (grepl("_A$", name))
          v <- c(pos - (end - start + 1) / 5, pos)
        else
          v <- c(pos, pos + (end - start + 1) / 5)

        if (strand == "+")
          return(v)
        else
          return(rev(v))
      }

      ac <- arr_coords(bp1$name[i], bp1$strand[i], pos1)
      # print(ac)
      panel.arrows(ac[1], xylim1[2] - ac[1] + xylim1[1], ac[2], xylim1[2] - ac[2] + xylim1[1], lwd = 2, col = "magenta", length = 0.1)
      # ac <- bal_to_vrg(arr_coords(bp2$name[i], bp2$strand[i], pos2))
      # print(ac)
      # panel.arrows(ac[1], xylim1[2] - ac[1] + xylim1[1], ac[2], xylim1[2] - ac[2] + xylim1[1], lwd = 2, col = "darkgreen", length = 0.1)
    }
  }
}

square_plot_REFALT <- function(dataset1, genome, this_chrom, start = NA, end = NA, bin_size1, bin_size2 = bin_size1, distance_fit = NA, what = "observed", filter = "none", raster = NA, zmin = NA, zmax = NA, secondary.zmax = NA, colscale = "orangered", mark_loops = F, mark_breakpoints = F, scale = "(do not propagate 'scale' to levelplot)", normalize_by_coverage = F, ...)
{
  dataset2 <- dataset1
  map <- fread(paste0("zcat data/hic_raw_maps/", dataset1, "_", genome, "_", bin_size1, ".txt.gz"))
  map <- as.matrix(map)
  dimnames(map) <- NULL

  if (normalize_by_coverage)
  {
    n1 <- sqrt(apply(map, 1, sum))
    n2 <- sqrt(apply(map, 2, sum))
    map <- t(t(map / n1 * mean(n1)) / n2 * mean(n2))
  }

  mult <- extract_mult(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2)
  data <- supplement_limits(dataset1 = dataset1, dataset2 = dataset2, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2,
    distance_fit = distance_fit, mult1 = mult[1], mult2 = mult[2], what = what, filter = filter, zmin = zmin, zmax = zmax)
  zmin <- data$zmin
  zmax <- data$zmax
  at <- data$at
  suffix <- data$suffix

  data <- extract_plot_data(dataset1 = "HiC_DB_6-8h_VRG", dataset2 = "HiC_DB_6-8h_VRG", genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size2,
    distance_fit = distance_fit, mult1 = mult[1], mult2 = mult[2], chrom1 = this_chrom, start1 = start, end1 = end, what = what, filter = filter, zmin = zmin, zmax = zmax)
  midpointsx <- data$midpointsx
  midpointsy <- data$midpointsy
  xlim <- data$xlim
  ylim <- data$ylim
  xlab <- data$xlab

  # here the difference
  # data are saved by hiccup in a way that REF coordinate corresponds to a row, ALT coordinate to a column
  # we would like to plot VRG on x axis and BAL on y axis, hence t(...)
  map1 <- t(map[data$bins1x, data$bins1y])
  if (!is.na(zmax))
    map1 <- pmin(map1, zmax - 1e-6)
  if (!is.na(zmin))
    map1 <- pmax(map1, zmin + 1e-6) 
  # map1 <- map[seq_len(nrow(data$map1)), seq_len(ncol(data$map1))][data$bins1x, data$bins1y] * mult[1]

  if (is.na(raster))
    raster <- max(dim(map1)) > 202

  colfun <- colscale_to_colfun(colscale)

  fixed_levelplot(map1, midpointsx, midpointsy, main = paste0(dataset1, suffix), xlab = xlab, at = at, colfun = colfun,
    secondary = NULL, secondary.midpoints = midpointsy, xlim = xlim, ylim = ylim,
    panel = mypanel.levelplot.same_bin_sizes, panel.internal = ifelse(raster, panel.levelplot.raster, panel.levelplot), ...)

  if (mark_loops)
    square_plot_mark_loops(dataset1, dataset2, genome, this_chrom, start, end, bin_size1, bin_size2, ...)
  if (mark_breakpoints)
    square_plot_mark_breakpoints(genome, chrom1 = this_chrom, start1 = start, end1 = end, bin_size1 = bin_size1, bin_size2 = bin_size2, ...)
}

triangle_plot <- function(dataset1, genome, this_chrom, start = NA, end = NA, bin_size1, distance_fit = NA, what = "observed", filter = "none", raster = NA, zmin = NA, zmax = NA, colscale = "orangered", scale = "(do not propagate 'scale' to levelplot)", main = NULL, ...)
{
  mult <- extract_mult(dataset1 = dataset1, dataset2 = NA, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size1)
  data <- supplement_limits(dataset1 = dataset1, dataset2 = NA, bin_size1 = bin_size1, bin_size2 = bin_size1,
    distance_fit = distance_fit, mult1 = mult[1], what = what, filter = filter, zmin = zmin, zmax = zmax)
  zmin <- data$zmin
  zmax <- data$zmax
  at <- data$at

  if (is.null(main))
    main <- paste0(dataset1, data$suffix)

  data <- extract_plot_data(dataset1 = dataset1, dataset2 = NA, genome = genome, bin_size1 = bin_size1, bin_size2 = bin_size1, distance_fit = distance_fit, 
    mult1 = mult[1], chrom1 = this_chrom, start1 = start, end1 = end, what = what, filter = filter, zmin = zmin, zmax = zmax)
  map1 <- data$map1
  midpointsx <- data$midpointsx
  xlim <- data$xlim
  ylim <- data$ylim
  xlab <- data$xlab

  if (is.na(raster))
    raster <- max(dim(map1)) > 202
  colfun <- colscale_to_colfun(colscale)

  fixed_triangle_levelplot(map1, midpointsx, main = main, xlab = data$xlab, chrom = this_chrom,
    at = at, colfun = colfun, xlim = xlim, ylim = ylim,
    panel = mypanel.levelplot.triangle, panel.internal = ifelse(raster, panel.levelplot.raster, panel.levelplot), ...,
    dataset1 = dataset1, genome = genome, this_chrom = this_chrom, start = start, end = end, bin_size1 = bin_size1)
}

make_plot <- function(dataset1, dataset2, genome, this_chrom, start, end, bin_size1, bin_size2 = bin_size1, suffix = "", pdf.dim = 8, ...)
{
  pdf.fname <- paste0(getwd(), "/", dataset1, "_", dataset2, "_", this_chrom, "_", ifelse(this_chrom == "genome", "", ifelse(is.na(start), "", paste0(as.integer(start), "-", as.integer(end), "_"))), ifelse(bin_size1 == bin_size2, bin_size1 / 1e3, paste0(bin_size1 / 1e3, "kb_", bin_size2 / 1e3)), "kb", suffix, ".pdf")
  cat(paste0("Plotting to ", pdf.fname, "\n"))
  cairo_pdf(pdf.fname, width = pdf.dim, height = pdf.dim - 1, onefile = T)
  square_plot(dataset1, dataset2, genome, chrom1 = this_chrom, start1 = start, end1 = end, bin_size1 = bin_size1, bin_size2 = bin_size2, ...)
  dev.off()
}

make_balancer_plot <- function(dataset1, dataset2, genome, this_chrom, start, end, bin_size1, bin_size2 = bin_size1, suffix = "", pdf.dim = 8, ...)
{
  pdf.fname <- paste0(getwd(), "/balancer_", dataset1, "_", dataset2, "_", this_chrom, "_", ifelse(this_chrom == "genome", "", ifelse(is.na(start), "", paste0(as.integer(start), "-", as.integer(end), "_"))), ifelse(bin_size1 == bin_size2, bin_size1 / 1e3, paste0(bin_size1 / 1e3, "kb_", bin_size2 / 1e3)), "kb", suffix, ".pdf")
  cat(paste0("Plotting to ", pdf.fname, "\n"))
  cairo_pdf(pdf.fname, width = pdf.dim, height = pdf.dim - 1, onefile = T)
  square_balancer_plot(dataset1, dataset2, genome, this_chrom, start, end, bin_size1, bin_size2, ...)
  dev.off()
}

make_triangle_plot <- function(dataset1, genome, this_chrom, start, end, bin_size1, suffix = "", pdf.dim = 8, ...)
{
  pdf.fname <- paste0(getwd(), "/", dataset1, "_", this_chrom, "_", ifelse(this_chrom == "genome", "", ifelse(is.na(start), "", paste0(as.integer(start), "-", as.integer(end), "_"))), bin_size1 / 1e3, "kb", suffix, ".pdf")
  cat(paste0("Plotting to ", pdf.fname, "\n"))
  cairo_pdf(pdf.fname, width = pdf.dim, height = pdf.dim / 2 - 0.1, onefile = T)
  triangle_plot(dataset1, genome, this_chrom, start, end, bin_size1, ...)
  dev.off()
}

embeed_levelplot_into_triangle_ggplot <- function(plot, plot_limits, ggplot_limits)
{
  w <- convertUnit(unit(1, "npc"), "pt", "x", "dimension")
  viewport_x <- unit((mean(plot_limits) - ggplot_limits[1]) / diff(ggplot_limits), "npc")
  width_factor <- diff(plot_limits) / diff(ggplot_limits) * sqrt(0.5)
  pushViewport(viewport(x = viewport_x, width = w * width_factor, height = w * width_factor, angle = 45))

  gt <- grid.grabExpr(print(plot), warn = 0)
  names(gt$children) <- make.unique(names(gt$children))
  gt$childrenOrder <- make.unique(gt$childrenOrder)

  gt <- removeGrob(gt, ".background", grep = T, global = T)
  gt <- removeGrob(gt, ".border", grep = T, global = T)
  # gt <- removeGrob(gt, ".tick", grep = T, global = T)
  vp <- gt$childrenvp[[1]]$children[["plot_01.panel.1.1.vp"]]
  gt$childrenvp[[1]]$children[["plot_01.panel.1.1.vp"]] <- viewport(name = "plot_01.panel.1.1.vp", xscale = vp$xscale, yscale = vp$yscale)
  grid.draw(gt, recording = FALSE)
  # grid.rect(gp = gpar(col = "red", fill = NA, lwd = 2))

  popViewport()
  return(NULL)
}

drawDetails.triangle_ggplot_grob <- function(x, recording = TRUE)
{
  # check the viewport dimensions
  w <- convertUnit(unit(1, "npc"), "pt", "x", "dimension")
  h <- convertUnit(unit(0.5, "npc"), "pt", "y", "dimension")

  # calculate the expanded limits for the x axis, cf. scale_x_continuous(...)
  expand <- x$params$expand
  limits <- x$params$limits
  expansion <- expand[1] * diff(limits) + expand[2]
  expanded_limits <- limits + c(-1, 1) * expansion
  # print(expanded_limits)

  # then, extend the limits further according to the height of the plot, so that no empty triangles appear in upper corners
  extension <- diff(limits) * convertUnit(h, "pt", valueOnly = T) / convertUnit(w, "pt", valueOnly = T)
  # print(extension)
  # no idea why you need to add one more bin width here! is it a bug in ggbio?
  extended_limits <- expanded_limits + c(-1, 1) * extension + c(0, 1) * x$params$bin_size1
  # print(extended_limits)

  # round up to the bin range etc.
  x$params <- modifyList(x$params, list(chrom1 = x$params$this_chrom, start1 = extended_limits[1], end1 = extended_limits[2]))
  mult <- do.call(extract_mult, x$params)
  x$params <- modifyList(x$params, list(mult1 = mult[1], mult2 = mult[2]))
  data <- do.call(supplement_limits, x$params)
  at <- data$at
  x$params$zmin <- data$zmin
  x$params$zmax <- data$zmax
  # x$params$at <- data$at

  data <- do.call(extract_plot_data, x$params)
  map1 <- data$map1
  midpoints <- data$midpointsx
  xlim <- data$xlim
  ylim <- data$ylim

  # truncate the matrix at the calculated distance off the diagonal
  margin_top <- 2 * extension / x$params$bin_size1
  map1[col(map1) - row(map1) - 1L > margin_top] <- NA
  # truncate also on the left and right
  margin_left <- 2 * (expanded_limits[1] - data$xlim[1]) / x$params$bin_size1
  map1[row(map1) + col(map1) < margin_left] <- NA
  margin_right <- 2 * (data$xlim[2] - expanded_limits[2]) / x$params$bin_size1
  # as above, no idea why you need to add +1L here! is it a bug in ggbio?
  map1[nrow(map1) + ncol(map1) - row(map1) - col(map1) + 2L + 1L < margin_right] <- NA

  if (is.na(x$params$raster))
    x$params$raster <- diff(limits) / x$params$bin_size1 > 202
  colfun <- colscale_to_colfun(x$params$colscale)

  pl <- levelplot(t(map1), xlab = NULL, ylab = NULL, col.regions = colfun(length(at) - 1), at = at, xlim = xlim, ylim = ylim,
    row.values = midpoints, column.values = rev(midpoints), midpoints = midpoints,
    panel = mypanel.levelplot.triangle, panel.internal = ifelse(x$params$raster, panel.levelplot.raster, panel.levelplot),
    colorkey = F, scales = list(x = list(), draw = F), par.settings = list(axis.line = list(col = NA)),
    background = NA, mark_loops = F, mark_breakpoints = x$params$mark_breakpoints, lwd.breakpoints = x$params$lwd.breakpoints,
    lty.breakpoints = x$params$lty.breakpoints, lineend.breakpoints = x$params$lineend.breakpoints, col.breakpoints = x$params$col.breakpoints,
    mark_domains = F, dataset1 = NA, genome = x$params$genome,
    this_chrom = x$params$this_chrom, start = x$params$start1, end = x$params$end1, bin_size = x$params$bin_size1,
    xscale.components = function(...) { tmp.xscale <<- xscale.components.default(...); list(num.limit = tmp.xscale$num.limit) },
    yscale.components = function(...) { tmp.yscale <<- yscale.components.default(...); list(num.limit = tmp.yscale$num.limit) },
  )

  embeed_levelplot_into_triangle_ggplot(pl, data$xlim, limits)
}

triangle_ggplot <- function(dataset, genome, bin_size, which, expand = c(0.05, 0), use_mult = T, ...)
{
  stopifnot(class(which) == "GRanges")
  stopifnot(length(which) == 1)
  this_chrom <- as.character(seqnames(which))
  limits <- c(start(which), end(which))
  stopifnot(length(expand) == 2)

  grob <- grob(cl = "triangle_ggplot_grob")
  grob$params <- modifyList(list(distance_fit = NA, what = "observed", filter = "none", raster = NA, zmin = NA, zmax = NA, colscale = "orangered", scale = "(do not propagate 'scale' to levelplot)", main = NULL,
    mark_breakpoints = F, lwd.breakpoints = 1, lty.breakpoints = "dashed", lineend.breakpoints = "butt", col.breakpoints = "magenta"),
    list(dataset1 = dataset, genome = genome, dataset2 = NA, bin_size1 = bin_size, bin_size2 = bin_size, this_chrom = this_chrom,
    limits = limits, expand = expand, ...))

  if (use_mult)
    mult <- extract_mult(dataset1 = dataset, dataset2 = NA, genome = genome, bin_size1 = bin_size, bin_size2 = NA)
  else
    mult <- c(1.0, 1.0)
  data <- supplement_limits(dataset1 = dataset, dataset2 = NA, bin_size1 = bin_size, distance_fit = grob$params$distance_fit, mult1 = mult[1],
    what = grob$params$what, filter = grob$params$filter, zmin = grob$params$zmin, zmax = grob$params$zmax)
  zmin <- data$zmin
  zmax <- data$zmax
  # message("[", zmin, ", ", zmax, "]")
  colvec <- colscale_to_colvec(grob$params$colscale)

  p <- ggplot(data.frame(fill = numeric(0))) +
    geom_tile(aes(fill = fill)) +
    scale_fill_gradientn(colors = colvec, limits = c(zmin, zmax), name = NULL) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.background = element_blank()
    ) +
    xlab(NULL) +
    scale_x_continuous(limits = limits, expand = expand) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    annotation_custom(grob = grob, xmin = limits[1], xmax = limits[2], ymin = -1, ymax = 1)

  return(p)
}

export_track <- function(dataset, genome, this_chrom, position, bin_size, prefix = NULL, suffix = NULL, backscale = T, trans_contacts = F)
{
  require(GenomicRanges, quietly = T)
  require(rtracklayer, quietly = T)

  dataset <- gsub("_iterative", "", dataset)
  ensure_loaded(dataset, genome, bin_size)
  this_chrom_bin <- as.integer((position - 1L) %/% bin_size + 1)
  if (trans_contacts)
  {
    chrom_map <- chrom_maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]]
    m <- match(this_chrom, chrom_map$chrom)
    bin <- this_chrom_bin + chrom_map$map_start[m] - 1L
  }
  else
    bin <- this_chrom_bin
  midpoint <- (this_chrom_bin - 0.5) * bin_size + 0.5

  mult <- extract_mult(dataset1 = dataset, genome = genome, bin_size1 = bin_size)
  mult <- mult[1]
  if (!is.finite(mult) | !backscale)
    mult <- 1.0
  print(mult)

  # wig.fname <- paste0(getwd(), "/", ifelse(is.null(prefix), "", paste0(prefix, "_")), dataset, "_", this_chrom, "_", as.integer(midpoint), "_", ifelse(is.null(suffix), "", paste0(suffix, "_")), bin_size / 1e3, "kb_ratio.wig")
  # cat(paste0("Exporting track to ", wig.fname, "\n"))
  # map <- maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[this_chrom]]$observed[bin, ] / maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[this_chrom]]$expected[bin, ]
  # map[1:bin] <- maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[this_chrom]]$observed[1:bin, bin] / maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[this_chrom]]$expected[1:bin, bin]
  # len <- length(map)

  # map[is.na(map)] <- -1
  # # map[bin] <- 1e6
  # gr <- GRanges(this_chrom, IRanges((1:len - 1) * bin_size + 1, 1:len * bin_size), score = map)
  # export.wig(gr, wig.fname)

  wig.fname <- paste0(getwd(), "/", ifelse(is.null(prefix), "", paste0(prefix, "_")), dataset, "_", this_chrom, "_", as.integer(midpoint), "_", ifelse(is.null(suffix), "", paste0(suffix, "_")), bin_size / 1e3, "kb_observed.wig")
  cat(paste0("Exporting track to ", wig.fname, "\n"))

  bins <- map_bins[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][, c("chrom", "start")]
  bins[, end := start + bin_size - 1L]

  if (trans_contacts)
  {
    where <- "genome"
    gr <- GRanges(bins)
  }
  else
  {
    where <- this_chrom
    gr <- GRanges(subset(bins, chrom == this_chrom))
  }

  map <- maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[where]]$observed[bin, ] * mult
  map[1:bin] <- maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[where]]$observed[1:bin, bin] * mult
  map[is.na(map)] <- -1

  stopifnot(length(gr) == length(map))
  gr$score <- map
  export.wig(gr, wig.fname)

  # wig.fname <- paste0(getwd(), "/", ifelse(is.null(prefix), "", paste0(prefix, "_")), dataset, "_", this_chrom, "_", as.integer(midpoint), "_", ifelse(is.null(suffix), "", paste0(suffix, "_")), bin_size / 1e3, "kb_expected.wig")
  # cat(paste0("Exporting track to ", wig.fname, "\n"))
  # map <- maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[this_chrom]]$expected[bin, ] * mult
  # map[1:bin] <- maps[[paste0(dataset, "_", genome)]][[as.character(bin_size)]][[this_chrom]]$expected[1:bin, bin] * mult
  # len <- length(map)

  # map[is.na(map)] <- -1
  # # map[bin] <- 1e6
  # gr <- GRanges(this_chrom, IRanges((1:len - 1) * bin_size + 1, 1:len * bin_size), score = map)
  # export.wig(gr, wig.fname)
}

make_all_plots <- function(dataset1, dataset2, genome, bin_size1, bin_size2 = bin_size1, by = 5e6, pdf.dim = 8, suffix = "", ...)
{
  pdf.fname <- paste0(getwd(), "/", dataset1, "_", dataset2, "_", ifelse(bin_size1 == bin_size2, bin_size1 / 1e3, paste0(bin_size1 / 1e3, "kb_", bin_size2 / 1e3)), "kb", suffix, ".pdf")
  cat(paste0("Plotting to ", pdf.fname, "\n"))
  cairo_pdf(pdf.fname, width = pdf.dim, height = pdf.dim, onefile = T)

  chrom_map <- read.table(paste0(master.dir, "/data/hic_maps/chrom_map_", genome, "_", bin_size1, ".txt"), header = T, sep = "\t")

  for (this_chrom in chrom_lst)
  {
    len <- chrom_map$length[match(this_chrom, chrom_map$chrom)]
    s <- seq(1, max(1, len - 0.2 * by), 0.8 * by)
    s[length(s)] <- max(1, len - by + 1)
    for (ss in s)
    {
      print(c(this_chrom, ss, min(ss + by - 1, len)))
      print(square_plot(dataset1, dataset2, genome, this_chrom, ss, min(ss + by - 1, len), bin_size1, bin_size2, ...))
    }
  }

  dev.off()
}

export_tracks_from_bedfile <- function(dataset, genome, bin_size, bedfile)
{
  bed <- import.bed(bedfile)

  for (i in seq_along(bed))
  {
    print(bed$name[i])
    export_track("HiC_3-4h", genome, as.vector(seqnames(bed))[i], (start(bed)[i] + end(bed)[i]) / 2, 10e3, bed$name[i])
  }
}
