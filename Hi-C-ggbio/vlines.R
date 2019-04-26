#' Map an arbitrary position within [start, end]
#' to the gene space (given by distribute_coords() )
#' df: data frame that contains genes-to-columns mapping (start -> start_)
#' bp: position to be mapped.
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
    bp.ir <- IRanges(bp, bp+1)
    if (overlapsAny(bp.ir, df.ir))
    {
        idx <- which(overlapsAny(df.ir, bp.ir))[1]
        coordinate = (df$start_[idx] + df$end_[idx])/2
    } else {
        idx <- precede(bp.ir, df.ir, select = "first")[1]
        coordinate = df$start_[idx]
    }
    coordinate
}




#'  internal function to draw vertical lines as grid system annotation
#'  note: this should be a global function (<<-), because its binding occurs only when the ggplot object is printed
drawDetails.vline_grob <<- function(x, recording = TRUE)
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
#'     <- if df is given, the plot produces a zig-zag shape (only meaningful
#'        in the genes-to-columns mapping panel).
#'  ymax: how far up to draw (in multiplicities of current plot height)
#'  gp: gpar() graphical parameters to use while drawing
ggplot_arbitrary_vertical_lines <- function(which, xintercepts, df = NULL, ymax = 1, gp = NULL)
{
    stopifnot(class(which) == "GRanges")
    stopifnot(length(which) == 1)
    
    x <- xintercepts
    
    if (is.null(x) || length(x)<1) return(NULL)
    else 
    {
        # assertions
        assert_that(is.numeric(x))
        if (!is.null(df))
        {
            assert_that(is.data.frame(df))
            assert_that(nrow(df)>0)
            assert_that('start' %in% colnames(df))
            assert_that('end' %in% colnames(df))
            assert_that('start_' %in% colnames(df))
            assert_that('end_' %in% colnames(df))
        }
        
        coords <- NULL
        for (i in seq_along(x))
        {
            if (!is.null(df))
            {
                if (x[i] <= min(df$start)) xbase <- start(wh)
                else if (x[i] >= max(df$end)) xbase <- end(wh)
                else xbase <- map_pos_to_gene_space(df,  x[i])
                coords <- rbind(coords, data.table(x = c(xbase, xbase, x[i], x[i]), 
                                                   y = c(0, 0.46 - 0.06*3, 0.96 - 0.06*3, ymax), id = i))
            } else {
                coords <- rbind(coords, data.table(x = c(x[i],x[i]), y=c(0,ymax), id = i))
            }
        }
        
        limits <- c(start(which), end(which))
        grob <- grob(cl = "vline_grob")
        grob$params <- list(limits = limits, ymax = ymax, gp = gp, coords = coords)
        
        return(annotation_custom(grob = grob, xmin = limits[1], xmax = limits[2]))
    }
}

