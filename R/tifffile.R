
#' @title OME Specification compatible images
#' @description
#' Class for working with images via RBioFormats and EBImage. Responsible for
#' delayed loading of the image at appropriate resolutions from a OME image.
#' @slot file_path full path to ome image
#' @slot extent named numeric vector in order of xmin, xmax, ymin, ymax.
#' This is the extent/coordinate mapping of the image
#' @slot dimorder character vector. Order of ome dimensions to use when 
#' slicing. Only dims present in the image should be provided.
#' Default is c("X", "Y", "C") but this may not be correct.
#' @slot resdims data.frame of pyramid levels and associated image dims info
#' colnames are `res`, `sizeX`, `sizeY`
#' @slot mapping list. Mapping of cropped region to represent in-memory vs the
#' true pixel dimensions of the file on-disk. Facilitates reading values only
#' from the cropped region.
setClass(
    "omeImage",
    slots = c(
        file_path = "character",
        extent = "ANY",
        dimorder = "character",
        resdims = "data.frame",
        mapping = "list"
    )
)



## omeImage init ####
setMethod(
    "initialize", signature("omeImage"), 
    function(.Object, ...
    ) {
        obj <- callNextMethod(.Object, ...)
        
        package_check(
            pkg_name = c("tifffile", "zarr"), 
            repository = c("pip:tifffile", "pip:zarr")
        )
        
        if (is.null(obj@extent)) {
            obj@extent <- .ext_to_num_vec(ext())
            names(obj@extent) <- c("xmin", "xmax", "ymin", "ymax")
        }
        
        obj@file_path <- p <- normalizePath(obj@file_path)
        # return early if no linked ome image
        if (length(p) == 0L) return(obj)
        
        # image and ome checks
        checkmate::assert_file_exists(p)
        if (!"ome" %in% file_extension(p)) {
            stop("Only .ome type images are allowed for `omeImage`")
        }
        
        # metadata
        # slicing order for .ome_zarr_get().
        # A default is set, but may be wrong
        obj@dimorder <- obj@dimorder %null% c("X", "Y", "C")
        # pyramids
        if (nrow(obj@resdims) == 0L) obj@resdims <- .ome_dim_info(p)
        # init on-disk vs in-mem spatial mapping (crop) 
        if (length(obj@mapping) == 0L) {
            pxdim <- .ome_px_dim(p)
            obj@mapping <- list(
                pxdim = pxdim,
                active = c(
                    xmin = 1, xmax = pxdim[["x"]],
                    ymin = 1, ymax = pxdim[["y"]]
                )
            )
        }
        # active (cropped) vs ondisk true dims
        m <- obj@mapping
        obj@mapping$ratios <- 
            m$active / c(rep(m$pxdim[["x"]], 2), rep(m$pxdim[["y"]], 2))
        
        
        return(obj)
    }
)




# python library ####

.pylib_tiffile <- function() {
    pylib <- getOption("giotto.tifffile", NULL)
    if (is.null(pylib)) {
        pylib <- reticulate::import("tifffile", delay_load = TRUE)
        options("giotto.tifffile" = pylib)
    }
    return(pylib)
}

.pylib_zarr <- function() {
    pylib <- getOption("giotto.zarr", NULL)
    if (is.null(pylib)) {
        pylib <- reticulate::import("zarr", delay_load = TRUE)
        options("giotto.zarr" = pylib)
    }
    return(pylib)
}


# image access ####

# Perform image access through tifffile and zarr.
# tifffile provides access to zarr to the on-disk ome.tif




.ome_tiff_zarr_store <- function(path) {
    package_check(
        pkg_name = c("tifffile", "zarr"), 
        repository = c("pip:tifffile", "pip:zarr")
    )
    TIF <- .pylib_tiffile()
    TIF$imread(path, aszarr = TRUE)
}

.ome_px_meta <- function(path) {
    m <- ometif_metadata(path, node = "Pixels")
    res <- as.list(m)[[1]]
    res <- as.list(res)
    names(res) <- rownames(m)
    return(res)
}

.ome_px_dim <- function(path) {
    pxmeta <- .ome_px_meta(path)
    pxdim <- as.numeric(c(pxmeta$SizeX, pxmeta$SizeY))
    names(pxdim) <- c("x", "y")
    return(pxdim)
}


#' @name ome_zarr_get
#' @title Get values from ome image zarr store
#' @description
#' Based on ome image dimensions (XYCZT) and pyramidal res level, retrieve
#' values from a zarr store. Works through reticulate python zarr. The reader
#' providing access to zarr is currently tifffile, but this may become
#' selecteable in the future depending on ome.zarr development.
#' @returns integer array with 2 or more dims
.ome_zarr_get <- function(
        path, res, shape = c("X", "Y", "C"), xrange, yrange, zplane, channels, time
) {
    ZARR <- .pylib_zarr()
    shape <- match.arg(
        unique(toupper(shape)), 
        choices = c("Y", "X", "C", "Z", "T"), 
        several.ok = TRUE
    )
    # convert for py vs R convention. No actual difference in arrays
    shape[match(c("X", "Y"), shape)] <- c("Y", "X")
    
    # slicing of zarr requires directly calling the builtins
    # otherwise it is not interpreted properly
    .py_slice <- reticulate::import_builtins()$slice
    .range_check <- function(x, dim_id) {
        if (length(x) != 2L) stop(
            wrap_txtf("%s must be a numeric of length 2 defining start/end",
                      dim_id)
        )
        as.integer(x) + c(-1L, 0L)
    }
    
    # get image pixel dimension characteristics
    px <- .ome_px_meta(path)
    
    # slicing info:
    # z[res][passed by `shape`]
    # res = pyramid res,
    # T = time Z = z plane, Y = y value, X = x value, C = channel
    # This is tifffile indexing convention (different from px$DimensionOrder)
    # 
    # These args are NOT positional. You may have more or fewer dimensions
    # to deal with depending on the structure of the image. The above simply
    # shows the expected dimension identities when all dimensions are provided
    # in the dataset. You may have fewer, such as only res, Y, X, C.
    # Omitting a slice for a dimension present in the image means taking the
    # whole range.
    
    # setup slicing args
    slice_args <- list(x = NULL) # earmark 1st arg for x
    res <- as.integer(res) - 1L # resolution
    
    # T
    if (!missing(time)) {
        time <- as.integer(time) - 1L
    } else {
        time <- seq_len(px$SizeT) - 1L
    }
    
    # Z
    if (!missing(zplane)) {
        zplane <- as.integer(zplane) - 1L
    } else {
        zplane <- seq_len(px$SizeZ) - 1L
    }
    
    # Y
    if (!missing(yrange)) {
        yrange <- .range_check(yrange, "yrange")
    } else {
        yrange <- c(0L, as.integer(px$SizeY))
    }
    
    # X
    if (!missing(xrange)) {
        xrange <- .range_check(xrange, "xrange")
    } else {
        xrange <- c(0L, as.integer(px$SizeX))
    }
     
    # C
    if (!missing(channels)) {
        channels <- as.integer(channels) - 1L
    } else {
        channels <- seq_len(px$SizeC) - 1L
    }
    
    # open image as zarr, close when function exits
    store <- .ome_tiff_zarr_store(path)
    z <- ZARR$open(store, mode = "r")
    on.exit(store$close, add = TRUE)
    slice_args$x <- z[res]
    
    ndim <- length(z[res]$shape) # ndims available in dataset
    if (ndim != length(shape)) {
        warning(wrap_txt(
            ".ome_zarr_get(): ndim of dataset != length of `shape`"
        ))
    }
    
    # assemble slice argslist
    for (d in shape) {
        slice_args[[length(slice_args) + 1L]] <- switch(d,
            "T" = time,
            "Z" = zplane,
            "Y" = .py_slice(yrange[1], yrange[2]),
            "X" = .py_slice(xrange[1], xrange[2]),
            "C" = channels
        )
    }
    
    a <- do.call(`[`, args = slice_args) # returns array
    vmsg(.is_debug = TRUE, ".ome_zarr_get:
         dim order", shape, "
         arraydims", dim(a)
    )
    
    return(a) # this will need an aperm()
}

# Array transposition for arrays pulled in with .ome_zarr_get()
#' @param a array to process
#' @param shape ome channels read order
#' @param target desired order of channels (more than available is fine)
#' @key
.ome_aperm <- function(a, 
    shape = c("X", "Y", "C"),
    target = c("X", "Y", "Z", "C", "T")
) {
    ndim <- length(dim(a))
    if (ndim > 3) {
        stop("readin array has more 3 dims. Not supported")
    }
    shape <- match.arg(
        unique(toupper(shape)), 
        choices = c("Y", "X", "C", "Z", "T"), 
        several.ok = TRUE
    )
    
    reorder <- match(target, shape)
    reorder <- reorder[!is.na(reorder)]
    return(aperm(a, reorder))
}








# converters ####


#' @title Convert ome.tif to tif
#' @name ometif_to_tif
#' @description
#' Simple converter from .ome.tif to .tif format. Utilizes the python
#' \pkg{tifffile} package. Performs image conversions one page at a time.
#' Wrap this in a for loop or lapply for more than one image or page.
#' @param input_file character. Filepath to ome.tif to convert
#' @param output_dir character. Directory to write .tif to. Defaults to a new
#' folder in the directory called "tif_exports"
#' @param page numeric. Which page of the tif to open (if needed). If provided,
#' a "_%04d" formatted suffix will be added to the output filename.
#' @param overwrite logical. Default = FALSE. Whether to overwrite if the
#' filename already exists.
#' @returns returns the written filepath invisibly
#' @family ometif utility functions
#' @export
ometif_to_tif <- function(
        input_file,
        output_dir = file.path(dirname(input_file), "tif_exports"),
        page,
        overwrite = FALSE
) {
    a <- list(input_file = input_file)
    
    # get tifffile py
    package_check("tifffile", repository = "pip:tifffile")
    ometif2tif_path <- system.file(
        "python", "ometif_convert.py",
        package = "GiottoClass"
    )
    reticulate::source_python(ometif2tif_path)
    # ensure output directory exists
    if (!checkmate::test_directory_exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # tif page
    if (!missing(page)) {
        a$page <- as.integer(page)
        fname_page <- sprintf("_%04d", a$page)
    } else {
        a$page <- 1L # default to page 1
        fname_page <- ""
    }
    a$page <- a$page - 1L # zero indexed
    
    # decide output filename
    fname <- sub(
        pattern = ".ome.tif$", replacement = "",
        x = basename(input_file)
    )
    fpath <- file.path(
        output_dir, paste0(fname, fname_page, ".tif")
    )
    a$output_file <- fpath
    
    # handle overwrites
    if (file.exists(fpath)) {
        if (isTRUE(overwrite)) {
            unlink(fpath, force = TRUE) # if overwrite, delete original
        } else {
            stop(fpath, "already exists. Set overwrite = TRUE to replace.\n",
                 call. = FALSE
            )
        }
    }
    do.call(ometif_2_tif, args = a)
    return(invisible(fpath))
}



# metadata ####


#' @name ometif_metadata
#' @title Read metadata of an ometif
#' @description Use the python package tifffile to get the the XML metadata
#' of a .ome.tif file. The R package xml2 is then used to work with it to
#' retrieve specific nodes in the xml data and extract data.
#' @param path character. filepath to .ome.tif image
#' @param node character vector. Specific xml node to get. More terms can be
#' added to get a node from a specific hierarchy.
#' @param output character. One of "data.frame" to return a data.frame of the
#' attributes information of the xml node, "xmL" for an xml2 representation
#' of the node, "list" for an R native list (note that many items in the
#' list may have overlapping names that make indexing difficult), or 
#' "structure" to invisibly return NULL, but print the structure of the XML
#' document or node.
#' @returns list of image metadata information
#' @family ometif utility functions
#' @export
ometif_metadata <- function(
        path, node = NULL, output = c("data.frame", "xml", "list", "structure")
) {
    checkmate::assert_file_exists(path)
    package_check(
        pkg_name = c("tifffile", "xml2"),
        repository = c("pip:tifffile", "CRAN:xml2")
    )
    
    TIF <- .pylib_tiffile()
    img <- TIF$TiffFile(path)
    output <- match.arg(
        output, choices = c("data.frame", "xml", "list", "structure")
    )
    x <- xml2::read_xml(img$ome_metadata)
    
    if (!is.null(node)) {
        node <- paste(node, collapse = "/")
        x <- xml2::xml_find_all(
            x, sprintf("//d1:%s", node), 
            ns = xml2::xml_ns(x)
        )
    }
    
    switch(output,
        "data.frame" = {
           x = Reduce("rbind", xml2::xml_attrs(x))
           rownames(x) <- NULL
           x <- as.data.frame(x)
           return(x)
        }, 
        "xml" = return(x),
        "list" = return(xml2::as_list(x)),
        "structure" = {
           xml2::xml_structure(x)
           return(invisible())
        }
    )
}


# get table of pyramid resolutions and dimensions information.
.ome_dim_info <- function(p) {
    pyr_meta <- GiottoClass::ometif_metadata(p, output = "list", "M")
    pyr_res <- unlist(pyr_meta) |> 
        strsplit(" ") |> 
        vapply(as.numeric, FUN.VALUE = numeric(2L)) |>
        t() |>
        data.frame()
    pyr_res <- cbind(
        GiottoClass::ometif_metadata(p, output = "data.frame", "M"),
        pyr_res
    ) |>
        data.table::setDT() |> 
        data.table::setnames(c("res", "sizeX", "sizeY"))
    px <- .ome_px_dim(p)
    basedim <- c(0, px)
    basedim <- t(basedim)
    colnames(basedim) <- c("res", "sizeX", "sizeY")
    pyr_res <- rbind(basedim, pyr_res)
    pyr_res$res <- as.integer(pyr_res$res) + 1L
    return(pyr_res)
}


# crop ####

# return numeric vector of xmin, xmax, ymin, ymax in res-specific px dims
.ome_pyramid_crop <- function(x, res) {
    checkmate::assert_class(x, "omeImage")
    res <- as.integer(res)
    ratios <- x@mapping$ratios
    pdims <- data.table::setDT(x@resdims)
    getres <- res # change name for easier DT indexing
    select_dims <- as.integer(pdims[res == getres, c("sizeX", "sizeY")])
    crop_vals <- c(rep(select_dims[[1]], 2), rep(select_dims[[2]], 2))
    crop_vals <- round(crop_vals * ratios)
    crop_vals[crop_vals == 0] <- 1
    return(crop_vals)
}


.ome_decide_res <- function(
        x, maxsize = getOption("giotto.plot_img_max_sample", 5e5)
) {
    checkmate::assert_class(x, "omeImage")
    checkmate::assert_numeric(maxsize)
    ratios <- x@mapping$ratios
    pdims <- data.table::as.data.table(x@resdims) # performs a copy
    ratiodiff <- c(
        ratios["xmax"] - ratios["xmin"], 
        ratios["ymax"] - ratios["ymin"]
    )
    sizeX <- sizeY <- active_area <- NULL
    pdims[, active_area := (sizeX * ratiodiff[1]) * (sizeY * ratiodiff[2])]
    pdims[order(active_area, decreasing = TRUE)]
    select_res <- pdims[active_area < maxsize, res][[1]]
    
    # no res (with crop) falls below maxsize
    if (length(select_res) == 0L) {
        warning(
            "No pyramid res low enough. Getting more than", maxsize, "values"
        )
        select_res <- pdims[, tail(res)]
    }
    return(select_res)
}
