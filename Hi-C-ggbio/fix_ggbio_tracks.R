traverse_ggbio_gtable_map <- function(g, fun)
{
  if ("gtable" %in% class(g))
    c(fun(g), do.call(c, lapply(g$grobs, traverse_ggbio_gtable_map, fun = fun)))
  else if ("gTree" %in% class(g))
    do.call(c, lapply(g$children, traverse_ggbio_gtable_map, fun = fun))
  else
    list()
}

extract_margin_left <- function(g)
{
  m <- which(g$layout$name == "panel")
  if (length(m) > 0)
  {
    start <- 2L
    end <- max(g$layout$r[m]) - 1L
    stopifnot(start <= end + 1L)
    if (start == end + 1L)
      list(unit(0, "cm"))
    else
      list(sum(g$widths[start:end]))
  }
  else
    list()
}

extract_margin_right <- function(g)
{
  m <- which(g$layout$name == "panel")
  if (length(m) > 0)
  {
    start <- max(g$layout$r[m]) + 1L
    end <- length(g$widths) - 1L
    stopifnot(start <= end + 1L)
    if (start == end + 1L)
      list(unit(0, "cm"))
    else
      list(sum(g$widths[start:end]))
  }
  else
    list()
}

traverse_ggbio_gtable_modify <- function(g, fun, ...)
{
  if ("gtable" %in% class(g))
  {
    g$grobs <- lapply(g$grobs, traverse_ggbio_gtable_modify, fun = fun, ...)
    g <- fun(g, ...)
  }
  else if ("gTree" %in% class(g))
    g$children <- lapply(g$children, traverse_ggbio_gtable_modify, fun = fun, ...)

  return(g)
}

adjust_margin_left <- function(g, margin)
{
  m <- which(g$layout$name == "panel")
  if (length(m) > 0)
  {
    start <- 2L
    end <- max(g$layout$r[m]) - 1L
    stopifnot(start <= end + 1L)
    if (start == end + 1L)
      g$widths[1] <- margin
    else
      g$widths[1] <- margin - sum(g$widths[start:end])
  }

  g
}

adjust_margin_right <- function(g, margin)
{
  m <- which(g$layout$name == "panel")
  if (length(m) > 0)
  {
    start <- max(g$layout$r[m]) + 1L
    end <- length(g$widths) - 1L
    stopifnot(start <= end + 1L)
    if (start == end + 1L)
      g$widths[length(g$widths)] <- margin
    else
      g$widths[length(g$widths)] <- margin - sum(g$widths[start:end])
  }

  g
}

adjust_axis_bottom <- function(g)
{
  n <- length(g$grobs)
  gg_orig <- g$grobs[[n - 1]]$grobs[[2]]$grobs[[1]]$children[[1]]$children$layout

  gg <- gg_orig
  i <- which(gg$layout$name == "axis-b")
  gg$grobs[i] <- list(zeroGrob())
  gg$heights[i] = unit(0, "cm")
  g$grobs[[n - 1]]$grobs[[2]]$grobs[[1]]$children[[1]]$children$layout <- gg

  gg <- g$grobs[[n]]$grobs[[2]]
  j <- which(gg$layout$name == "axis-b")
  gg$grobs[j] <- gg_orig$grobs[i]
  gg$widths <- gg_orig$widths
  g$grobs[[n]]$grobs[[2]] <- gg

  g
}

fix_ggbio_tracks <- function(tr)
{
  stopifnot("Tracks" %in% class(tr))

  # # force the presence of x axis
  # tr@hasAxis[length(tr@hasAxis)] <- T

  g <- as(tr, "grob")

  max_margin_left <- do.call(max, traverse_ggbio_gtable_map(g, extract_margin_left))
  g <- traverse_ggbio_gtable_modify(g, adjust_margin_left, max_margin_left)

  max_margin_right <- do.call(max, traverse_ggbio_gtable_map(g, extract_margin_right))
  g <- traverse_ggbio_gtable_modify(g, adjust_margin_right, max_margin_right)

  # g <- adjust_axis_bottom(g)
  g
}
