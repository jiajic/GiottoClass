# tileIterator Class ####

setClass(
    "tileIterator",
    slots = list(
        extent = "numeric",
        n = "numeric",
        tiles = "array"
    )
)

setMethod("initialize", signature("tileIterator"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    
    # initialize n tiles
    if (length(.Object@n) == 0L) {
        .Object@n <- 0
    }
    
    # return early if extent not provided
    if (length(.Object@extent) == 0L) {
        return(.Object)
    }
    
    # check extent validity
    if (length(.Object@extent) != 4L) {
        stop("tileIterator: invalid extent information", call. = FALSE)
    }
    
    # return early if n tiles = 0
    if (.Object@n == 0) {
        return(.Object)
    }
    
    # generate tile extent array
    n_desired <- .Object@n
    e <- terra::ext(.Object@extent)
    .Object@tiles <- .chunk_plan(e, min_chunks = n_desired)
    
    return(.Object)
})

setMethod("show", signature("tileIterator"), function(object) {
    cat("Object of class", class(object), "\n")
    
    # no extent, return early
    e <- object@extent
    if (length(e) == 0L) {
        cat("<empty>")
        return(invisible())
    }
    
    d <- dim(object)
    plist <- list(
        extent = sprintf(
            "%s (xmin, xmax, ymin, ymax)", 
            paste(.ext_to_num_vec(e), collapse = ", ")
        ),
        dim = paste(dim(x), collapse = " ")
    )
    print_list(plist)
})

setMethod("nrow", signature("tileIterator"), function(x) {
    nrow(x@tiles)
})

setMethod("ncol", signature("tileIterator"), function(x) {
    res <- ncol(x@tiles)
    if (is.na(res)) {
        res <- 0 # catch for when `x@tiles` is not an array
    }
    return(res)
})

setMethod("length", signature("tileIterator"), function(x) {
    nrow(x) * ncol(x)
})

setMethod("length<-", signature("tileIterator"), function(x, value) {
    x@n <- value
    return(initialize(x))
})

setMethod("dim", signature("tileIterator"), function(x) {
    c(nrow(x), ncol(x))
})

setMethod("ext", signature("tileIterator"), function(x, ...) {
    if (length(x@extent) == 0L) {
        stop("tileIterator: No extent set", call. = FALSE)
    }
    ext(x@extent, ...)
})

setMethod("ext<-", signature("tileIterator"), function(x, value) {
    x@extent <- .ext_to_num_vec(ext(value))
    return(initialize(x))
})

setMethod("[", signature(x = "tileIterator", i = "numeric", j = "missing", drop = "missing"), function(x, i) {
    i <- as.integer(i)
    if (any(i > length(x) | i <= 0)) stop("tileIterator: subscript out of bounds", call. = FALSE)
    
    i_idx <- floor(i / ncol(x)) + 1L
    no_resid <- i %% ncol(x) == 0L
    i_idx[no_resid] <- i_idx[no_resid] - 1L
    j_idx <- i %% ncol(x)
    j_idx[j_idx == 0L] <- ncol(x)
    
    mapply(function(i, j) {
        ext(x@tiles[i, j,])
    }, i_idx, j_idx)
})

setMethod("[", signature(x = "tileIterator", i = "numeric", j = "numeric", drop = "missing"), function(x, i, j) {
    mapply(function(i, j) {
        ext(x@tiles[i, j,])
    }, i, j)
})

setMethod("[", signature(x = "tileIterator", i = "missing", j = "missing", drop = "missing"), function(x) {
    x[seq_len(length(x))]
})

setMethod("plot", signature(x = "tileIterator", y = "missing"), function(x, ...) {
    if (length(x@tiles) == 0L) {
        stop("No tiles to plot.\nTry requesting tiles with `length()`")
    }
    
    .preview_chunk_plan(x[], mode = "poly", ...)
})

# helper functions ####

#' @name .get_dim_n_chunks
#' @title Get rows and cols needed to create at least n chunks from given extent
#' @description Algorithm to determine how to divide up a provided extent into
#' at least \code{n} different chunks. The chunks are arranged so as to prefer
#' being as square as posssible with the provided dimensions and minimum n chunks.
#' @param n minimum n chunks
#' @param e selection extent
#' @examples
#' e <- terra::ext(0, 100, 0, 100)
#' .get_dim_n_chunks(n = 5, e = e)
#' @seealso \code{\link{.chunk_plan}}
#' @return numeric vector of x and y stops needed
.get_dim_n_chunks = function(n, e) {
    # find x to y ratio as 'r'
    e = e[]
    r = (e[['xmax']] - e[['xmin']]) / (e[['ymax']] - e[['ymin']])
    
    # x * y = n = ... ry^2 = n
    y = ceiling(sqrt(n / r))
    x = ceiling(n / y)
    
    return(c(y, x))
}


#' @name .chunk_plan
#' @title Plan spatial chunking extents
#' @description
#' Generate the individual extents that will be used to spatially chunk a set of
#' data for piecewise and potentially parallelized processing. Chunks will be
#' generated first by row, then by column. The chunks try to be as square as
#' possible since downstream functions may require slight expansions of the
#' extents to capture all parts of selected polygons. Minimizing the perimeter
#' relative to area decreases waste.
#' @param extent terra SpatExtent that covers the region to spatially chunk
#' @param nrows,ncols numeric. nrow/ncol must be provided as a pair. Determines how many
#' rows and cols respectively will be used in spatial chunking. If NULL, min_chunks
#' will be used as an automated method of planning the spatial chunking
#' @param min_chunks numeric. minimum number of chunks to use.
#' @seealso \code{\link{.get_dim_n_chunks}}
#' @examples
#'  e <- ext(0, 100, 0, 100)
#'
#' a <- .chunk_plan(e, min_chunks = 9)
#' plot(e)
#' plot(ext(a[1,1,]), add = T)
#' plot(ext(a[1,2,]), add = T)
#' plot(ext(a[2,1,]), add = T)
#' @keywords internal
#' @return 3D array of extents to use, organized by row, col, and ext bounds
.chunk_plan = function(extent, min_chunks = NULL, nrows = NULL, ncols = NULL) {
    checkmate::assert_class(extent, 'SpatExtent')
    if(!is.null(nrows)) checkmate::assert_true(length(c(nrows, ncols)) == 2L)
    else {
        checkmate::assert_numeric(min_chunks)
        res = .get_dim_n_chunks(n = min_chunks, e = extent)
        nrows = res[1L]
        ncols = res[2L]
    }
    
    x_stops = seq(
        from = terra::xmin(extent), 
        to = terra::xmax(extent), 
        length.out = ncols + 1L
    )
    y_stops = seq(
        from = terra::ymin(extent), 
        to = terra::ymax(extent), 
        length.out = nrows + 1L
    )
    
    # vector of extent values
    e_vec <- c()
    for (i in seq_len(nrows)) {
        for (j in seq_len(ncols)) {
            e_vec <- c(
                e_vec, 
                x_stops[j], 
                x_stops[j + 1L], 
                y_stops[i], 
                y_stops[i + 1L]
            )
        }
    }
    
    a <- array(e_vec, dim = c(4, ncols, nrows))
    a <- aperm(a, perm = c(3, 2, 1 ))
    # reverse order of rows so tiles count from top to bottom
    a <- a[seq(from = nrow(a), to = 1),,]
    
    return(a)
}




#' @name .preview_chunk_plan
#' @title Plot a preview of the chunk plan
#' @description
#' Plots the output from \code{\link{.chunk_plan}} as a set of polygons to preview.
#' Can be useful for debugging. Invisibly returns the planned chunks as a SpatVector
#' of polygons
#' @param extent_list list of extents from \code{.chunk_plan}
#' @keywords internal
.preview_chunk_plan = function(extent_list, mode = c('poly', 'bound')) {
    checkmate::assert_list(extent_list, types = 'SpatExtent')
    mode = match.arg(mode, choices = c('poly', 'bound'))
    
    switch(mode,
           'poly' = {
               poly_list = sapply(extent_list, terra::as.polygons)
               poly_bind = do.call(rbind, poly_list)
               terra::plot(poly_bind, values = as.factor(seq_along(poly_list)))
               return(invisible(poly_bind))
           },
           'bound' = {
               xlim = c(extent_list[[1]]$xmin, extent_list[[length(extent_list)]]$xmax)
               ylim = c(extent_list[[1]]$ymin, extent_list[[length(extent_list)]]$ymax)
               # initiate plot
               plot(x = NULL, y = NULL, asp = 1L, xlim = xlim, ylim = ylim)
               # plot extent bounds
               for(e in extent_list) {
                   rect(e$xmin, e$ymin, e$xmax, e$ymax)
               }
               return(invisible())
           })
}
