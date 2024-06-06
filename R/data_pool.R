# giottoDataPool ####

giottoDataPool <- setClass(
    "giottoDataPool",
    slots = c(
        cell_ID = "character",
        feat_ID = "character",
        extent = "ANY", # stored as numeric
        data = "list"
    )
)

setMethod("show", signature("giottoDataPool"), function(object) {
    cat("<giottoDataPool>\n")
    elen <- length(object@extent)
    eprint <- if (elen < 1L) {
        "no spatial info"
    } else {
        sprintf("%d items; plot() to visualize", elen)
    }
    
    plist <- list(
        extent = eprint,
        spat_ids = length(object@cell_ID),
        feat_ids = length(object@feat_ID)
    )
    print_list(plist)

    data <- object@data
    if (length(data) == 0) {
        cat("data     : none")
        return(invisible())
    }
    data_pre <- format(c("data     :", rep("", length(data) - 1L)))
    obj_id <- format(names(data))
    obj_classes <- format(
        sprintf("`%s`", vapply(data, class, FUN.VALUE = character(1L)))
    )
    obj_names <- lapply(data, function(d) {
        res <- try(objName(d), silent = TRUE)
        if (inherits(res, "try-error")) res <- ""
        return(res)
    })
    data_lines <- sprintf(
        "%s * %s -- %s (\"%s\")", data_pre, obj_id, obj_classes, obj_names
    )
    cat(data_lines, sep = "\n")
})

setMethod("initialize", signature("giottoDataPool"), function(.Object, ...) {
    
    # get data
    data <- list(...)
    if (length(data) > 0L) {
        .Object@data <- c(.Object@data, data)
    }
    
    # return early if no data
    if (length(.Object@data) == 0L) return(.Object)
    
    # apply identifiers
    .Object@data <- .gdp_id(.Object@data)
    
    .scan_data <- function(fun) {
        info <- lapply(.Object@data, function(.data) {
            res <- try({
                fun(.data)
            }, silent = TRUE)
            if (inherits(res, "try-error")) res <- NULL
            return(res)
        })
        # remove NULLs
        info <- info[lengths(info) != 0L]
        return(info)
    }
    
    # check spat and feat IDs
    sids <- .scan_data(spatIDs)
    fids <- .scan_data(featIDs)
    exts <- .scan_data(terra::ext)
    
    # get common values
    sids <- Reduce(intersect, sids)
    fids <- Reduce(intersect, fids)
    exts <- lapply(exts, terra::as.polygons) %>%
        Reduce(f = "rbind")
    
    # if any were found, add to slot
    if (length(sids) > 0L) .Object@cell_ID <- sids
    if (length(fids) > 0L) .Object@feat_ID <- fids
    if (length(exts) > 0L) .Object@extent <- exts
    return(.Object)
})

.DollarNames.giottoDataPool <- function(x, pattern) {
    names(x@data)
}

setMethod("$", signature("giottoDataPool"), function(x, name) {
    x@data[[name]]
})

setMethod("+", signature("giottoDataPool", "ANY"), function(e1, e2) {
    e1@data <- c(e1@data, e2)
    e1 <- initialize(e1)
    return(e1)
})

setMethod("plot", signature("giottoDataPool"), function(x) {
    elen <- length(x@extent)
    if (elen == 0L) {
        cat("No spatial info\n")
        return(invisible())
    }
    terra::plot(x@extent, border = rainbow(elen))
})

setMethod("spatIDs", signature("giottoDataPool"), function(x) {
    x@cell_ID
})

setMethod("featIDs", signature("giottoDataPool"), function(x) {
    x@feat_ID
})

# provide list of objects
# return named list of objects
.gdp_id <- function(data) {
    ndata <- length(data)
    dnames <- names(data)

    # no data
    if (ndata == 0) return(NULL)
    # no names
    if (is.null(dnames)) {
        names(data) <- sprintf("data_%d", ndata)
        return(data)
    }
    # add missing names
    if (any(dnames == "")) {
        nmissing <- sum(dnames == "")
        max_idx <- gsub("data_", "", dnames) %>%
            as.numeric() %>%
            max(na.rm = TRUE)
        new_idx <- seq(from = max_idx + 1, length.out = nmissing)
        new_ids <- sprintf("data_%d", new_idx)
        dnames[dnames == ""] <- new_ids
        names(data) <- dnames
    }
    return(data)
}




# queryGDP ####

# no spatial subsets happen unless extent is given

#' @title Query a set of requested data from a giottoDataPool
#' @name queryGDP
#' @param x `giottoDataPool` object
#' @param key character. Output should be keyed on "spat" (spatial) or "feat"
#' (feature) information
#' @param cell_ID (optional) character vector. If provided, subset to only
#' the cell_IDs given
#' @param feat_ID (optional) character vector. If provided, subset to only
#' the feat_IDs given
#' @param atts character. Attributes/columns/info that must be present in the
#' output
#' @param minimal logical. Output data.table will only any intrinsic data 
#' (IDs, geometries, edge interactions, etc.) and those attributes specified
#' through `atts`
#' @param cast character. data type to output as
#' @export
queryGDP <- function(
        x,
        key = c("spat", "feat"),
        cell_ID = NULL,
        feat_ID = NULL,
        extent = NULL,
        atts = NULL,
        minimal = TRUE,
        cast = "data.table"
) {
    checkmate::assert_class(x, "giottoDataPool")
    x <- initialize(x) # update object
    key <- match.arg(key, choices = c("spat", "feat"))
    data <- x@data
    
    # cropping is requested
    if (!is.null(extent)) {
        .gdp_validate_crop(gdp_ext_polys = x@extent, extent)
        data <- .gdp_crop(data, extent)
    }
    
    # attribute selection to make minimal report
    if (isTRUE(minimal)) {
        data <- .do_data_select(data = data, key = key, atts = atts)
    }
    
    # TODO collapse special types
    
    # join to cast type
    
    
}

# internals ####

## cropping ####

# ensure that there is information in the crop region
# check this by cropping the extents that exist in the object and seeing
# how many survive the crop. Pass the check if at least 1 extent is still
# left.
.gdp_validate_crop <- function(gdp_ext_polys, extent) {
    if (nrow(crop(gdp_ext_polys, extent)) == 0L) {
        stop("crop region is empty", call. = FALSE)
    }
}

# perform crop on the data and return. There are also non-spatial objects
# in here, so we use try() to affect only the spatial objects.
.gdp_crop <- function(data, extent) {
    lapply(data, function(d) {
        res <- try(crop(d, extent), silent = TRUE)
        # if does not work with crop, return original object
        # we check this in currently a crude manner by simply seeing if
        # crop() fails
        if (inherits(res, "try-error")) res <- d
        return(res)
    })
}



# selecting specific attributes from data to extract for the final merged report
# is only possible when you have pre-knowledge of which object contains that
# attribute. Giotto does not have this information and instead expects the
# attribute to just be present in the final report (usually along with many
# other attributes that are not used.)
# To get around that, we search for attributes among the data we have.
#
# returns a list of character vectors of attributes, named by the object ID
# they can be found in.

.att_search <- function(data, key = "spat", atts = NULL) {
    atts_list <- lapply(data, .att, key)
    
    atts_location <- list()
    for (a in atts) {
        for (d in names(atts_list)) {
            if (a %in% atts_list[[d]]) {
                atts_location[[d]] <- c(atts_location[[d]], a)
            }
        }
    }

    # error if any atts not found
    found_atts <- unlist(atts_location, recursive = TRUE)
    missing_atts <- atts[!atts %in% found_atts]
    if (length(missing_atts) > 0L) {
        stop(sprintf(
            "Not all requested `atts` were found:\n%s",
            paste0(missing_atts, collapse = ", ")
        ), 
        call. = FALSE)
    }
    found_dups <- duplicated(found_atts)
    if (any(found_dups)) {
        warning(sprintf(
            "Some `atts` were found in more than one object:\n%s",
            paste0(found_atts[found_dups], collapse = ", ")
        ))
    }
    
    return(atts_location)
}

.do_data_select <- function(data, key = "spat", atts = NULL) {
    sel_list <- .att_search( # if no atts requested, `list()` returned
        data = data, key = key, atts = atts
    )
    
    atts_data_names <- names(sel_list)
    other_data_names <- names(data)[!names(data) %in% atts_data_names]
    
    if (length(atts_data_names) > 0) {
        data[atts_data_names] <- lapply(atts_data_names, function(dn) {
            d <- data[[dn]]
            .sel(x = d, key = key, atts = sel_list[[dn]])
        })
    }
    
    data[other_data_names] <- lapply(other_data_names, function(dn) {
        d <- data[[dn]]
        .sel(x = d, key = key, atts = NULL)
    })
    
    return(data)
}

.gtrueclass <- function(x) {
    if (inherits(x, "giottoSubobject")) return(class(x[]))
    else return(class(x))
}

.gtrueinherits <- function(x, what) {
    if (inherits(x, "giottoSubobject")) return(inherits(x[], what))
    else return(inherits(x, what))
}


## internal generics ####
setGeneric(".att", function(x, ...) standardGeneric(".att"))
setGeneric(".sel", function(x, ...) standardGeneric(".sel"))

### .att ####
setMethod(".att", signature("ANY"), function(x, ...) names(x))
setMethod(".att", signature("exprObj"), function(x, key = "spat", ...) {
    key <- match.arg(key, choices = c("spat", "feat"))
    switch(key,
           "spat" = colnames(x),
           "feat" = rownames(x)
    )
})

### .sel ####
# in cases where the minimal report is generated, only return the
# core information and any geom data along with any atts specifically
# requested. If there are no atts requested for an object with no geom
# info then a NULL will be returned instead
setMethod(".sel", signature("cellMetaObj"), function(x, atts, ...) {
    if (length(atts) == 0L) return(NULL)
    x[] <- x[][, unique(c("cell_ID", atts))]
    return(x)
})
setMethod(".sel", signature("featMetaObj"), function(x, atts, ...) {
    if (length(atts) == 0L) return(NULL)
    x[] <- x[][, unique(c("feat_ID", atts))]
    return(x)
})
setMethod(".sel", signature("dimObj"), function(x, atts, ...) {
    if (length(atts) == 0L) return(NULL)
    x[] <- x[][, unique(atts), drop = FALSE]
    return(x)
})
setMethod(".sel", signature("exprObj"), function(x, atts, key = "spat", ...) {
    if (length(atts) == 0L) return(NULL)
    switch(key,
           "spat" = {
               x[] <- x[][atts, drop = FALSE, ]
           },
           "feat" = {
               x[] <- x[][, atts, drop = FALSE]
           }
    )
    x[] <- x[][, unique(c("feat_ID", atts))]
    return(x)
})
setMethod(".sel", signature("giottoPoints"), function(x, atts, ...) {
    x[, c("feat_ID", "feat_ID_uniq", atts)]
})
setMethod(".sel", signature("giottoPolygon"), function(x, atts, ...) {
    x <- x[, c("poly_ID", atts)]
    x@spatVectorCentroids <- x@spatVectorCentroids[, c("poly_ID", atts)]
    # TODO overlaps
    return(x)
})
setMethod(".sel", signature("spatLocsObj"), function(x, atts, ...) {
    x # no subsetting needed
})
setMethod(".sel", signature("spatEnrObj"), function(x, atts, ...) {
    if (length(atts) == 0L) return(NULL)
    x[, atts]
})



### .gcombine ####

# accepts a list of objects
# collapse some object types first to utilize their faster joins
.gcombine <- function(x, key = "spat", cell_ID, feat_ID) {
    
    # remove any NULLs
    x <- x[!is.null(x)]
    
    # dbData collapsing
    is_dbdata <- lapply(x, .gtrueinherits, "dbData")
    if (length(is_dbdata) > 0L) {
        dbd <- x[is_dbdata]
        x <- x[!is_dbdata]
        
        # TODO collapse of dbdata objects
    }
    
    # data.frame collapsing
    is_df <- lapply(x, .gtrueinherits, "data.frame")
    if (length(is_df) > 0L) {
        df <- x[is_df]
        x <- x[!is_df]
    }
    
    # any other classes
    if (length(x) > 0L) stop("unrecognized data types")
    
}

# NOTE: using any signatures means that classes have to be imported in
# this package

# TODO figure out collapsing of specific types. 
# May be done with superclasses with .gtrueinherits()
# Then how to merge them all into a final object


