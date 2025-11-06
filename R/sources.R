# classes ####

setClass(
    "giottoSource",
    contains = "VIRTUAL",
    slots = list(
        path = "character"
    ),
    prototype = list(
        path = NA_character_
    )
)

setClass("gDirSource",
    contains = "giottoSource",
    slots = list(
        content = "list",
        read = "function",
        write = "function"
    )
)


# generics ####

setGeneric(".gsource", function(x, ...) standardGeneric(".gsource"))
setGeneric(".gsource<-", function(x, ..., value) standardGeneric(".gsource<-"))
setGeneric(".gschema_name", function(x, ...) standardGeneric(".gschema_name"))

setMethod(".gsource", "giottoSource", function(x, ...) x@path)
setMethod(".gsource", "giotto", function(x, ...) x@source)
setMethod(".gsource<-", c("giotto", "ANY"), function(x, ..., value) {
    x@source <- value
    x
})

# update the content list from a gDirSource object with the default subdirectory
# names for specific formats.
gdirsource_defaults <- function(clist) {
    clist$stores$h5 <- clist$stores$h5 %null% "h5"
    clist$stores$spatvector <- clist$stores$spatvector %null% "spatvector"
    clist$stores$parquet <- clist$stores$parquet %null% "parquet"
    clist
}

setMethod("initialize", signature("gDirSource"), function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    p <- .Object@path <- normalizePath(.Object@path, mustWork = FALSE)
    json_path <- file.path(p, "giottodir.json")
    if (is.na(p)) {
        stop("[gDirSource] 'path' should be a directory path for the giotto project.\n")
    }

    # default content subdirectories
    content <- .Object@content
    content <- gdirsource_defaults(content)
    .Object@content <- content

    write_fun <- function() { # write a giotto dir json
        if (!dir.exists(p)) {
            message("Setting up Giotto project directory at:", p)
            dir.create(path = p)
        }
        jsonlite::write_json(.Object@content,
            path = json_path
        )
        invisible(TRUE)
    }
    .Object@write <- write_fun

    read_fun <- function() { # read a giotto dir json as gDirSource
        if (!dir.exists(p)) {
            stop("Path to Giotto project directory does not exist\n")
        }
        if (!file.exists(json_path)) {
            stop("giottodir.json not found in: ", p)
        }
        .Object@content <- jsonlite::fromJSON(json_path)
        .Object
    }
    .Object@read <- read_fun

    if (!dir.exists(p) || !file.exists(json_path)) {
        .Object@write()
    }

    .Object
})

setMethod("show", "gDirSource", function(object) {
    cat(sprintf("<%s>\n", class(object)))

    if (file.exists(object@path)) {
        object <- object@read()
        cat("stores:\n")
        print_list(object@content$stores)
        cat("\n")
        cat("artifacts:", length(object@content$artifacts), "\n")
        cat("versions:", length(object@content$versions), "\n")
    } else {
        cat("giottodir.json not written yet. Use `@write()`\n")
    }

    cat("\n")
    cat("* `@read(path)` to read from a Giotto directory json\n")
    cat("* `@write(path)` to write a Giotto directory json\n")
})

setMethod("$", "gDirSource", function(x, name) {
    file.path(x@path, x@content$stores[[name]])
})

#' @keywords internal
#' @export
.DollarNames.gDirSource <- function(x, pattern) names(x@content$stores)

setMethod("[", c(x = "gDirSource", i = "character", j = "missing", drop = "missing"), function(x, i, j, ..., drop) {
    file.path(x@path, x@content$stores[[i]])
})

# JSON records ####

# unique id generator
.make_uid <- function(n = 10) {
    paste(sample(c(letters, LETTERS, 1:9), n), collapse = "")
}
.timestamp <- function() {
    as.character(format(Sys.time()))
}

# these function return the giotto object

.giotto_json_edit_content <- function(gobject, fun) {
    gsrc <- gobject@source
    gsrc <- gsrc@read()
    gsrc@content <- fun(gsrc@content)
    gsrc <- initialize(gsrc)
    gsrc@write()
    gobject@source <- gsrc
    gobject
}

.giotto_json_add_artifact <- function(gobject, store, uid = .make_uid()) {
    checkmate::assert_character(store)
    checkmate::assert_character(uid)

    entry <- list(
        time = .timestamp(),
        store = store,
        version = NA_character_
    )

    .giotto_json_edit_content(gobject, function(x) {
        x$artifacts[[uid]] <- entry
        x
    })
}

.giotto_json_add_store <- function(gobject, name, path) {
    checkmate::assert_character(name)
    checkmate::assert_character(path)

    basepath <- gobject@source@path
    fullpath <- file.path(basepath, path)
    if (!dir.exists(fullpath)) {
        vmsg("creating", fullpath, "...")
        dir.create(fullpath)
    }

    entry <- list(
        path = path
    )

    .giotto_json_edit_content(gobject, function(x) {
        x$stores[[name]] <- entry
        x
    })
}

.giotto_json_add_project_version <- function(gobject,
    uid = paste0("giottosave_", .make_uid())) {
    checkmate::assert_character(uid)

    entry <- list(
        time = .timestamp()
    )

    .giotto_json_edit_content(gobject, function(x) {
        x$versions[[uid]] <- entry
        x
    })
}

.giotto_json_artifact_tag_version <- function(gobject, artifacts, version) {
    checkmate::assert_character(artifacts)
    checkmate::assert_character(version)

    .giotto_json_edit_content(gobject, function(x) {
        for (id in artifacts) {
            x$artifacts[[id]]$version <- version
        }
        x
    })
}

.eval_to_gdirsource <- function(x) {
    if (inherits(x, "gDirSource")) {
        src <- x
    } else if (inherits(x, "giotto")) {
        src <- x@source
    } else if (is.character(x)) {
        if (x == "") {
            stop("This is an in-memory giotto object.")
        }
        src <- new("gDirSource", path = x)
    } else {
        stop("`gDirSource` info can only be retrieved from a giotto object or the backing giotto directory",
             call. = FALSE)
    }
    src@read()
}

.giotto_json_artifacts <- function(x) {
    version <- store <- NULL # NSE var
    src <- .eval_to_gdirsource(x)
    artifacts <- src@content$artifacts
    artifacts <- data.table::rbindlist(artifacts, idcol = "id")
    stores <- src@content$stores
    stores <- data.frame(
        store = names(stores),
        path = unlist(stores),
        stringsAsFactors = FALSE
    )
    artifacts <- merge(artifacts, stores, by = "store")
    artifacts[, "fullpaths" := file.path(src@path, path, id)]
    artifacts
}

# pruning ####

.giotto_dir_prune <- function(x) {
    artifacts <- .giotto_json_artifacts(x)
    unversioned <- artifacts[is.na(version)]
    for (f in unversioned$fullpaths) {
        if (file.exists(f)) {
            file.remove(f)
        } else {
            vmsg("File not found (already deleted?):", f)
        }
    }
    vmsg("Removed ", nrow(unversioned), " unversioned artifacts")

    # update json
    keep_ids <- artifacts[!is.na(version), id]
    src <- .eval_to_gdirsource(x)
    src@content$artifacts <- src@content$artifacts[keep_ids]
    src@write()

    TRUE
}

# artifact writing ####

# gobject is giotto object with @source
# x is DelayedMatrix to save
# returns the updated matrix object and uid of the save
.delayedmatrix_save <- function(gobject, x) {
    # skip if no gDirSource
    gdir <- gobject@source
    if (is.character(gdir)) {
        if (gdir == "") return(list(x = x, uid = NA_character_))
    }
    # skip if not gDirSource -- may change for other backends
    if (!inherits(gdir, "gDirSource")) return(list(x = x, uid = NA_character_))

    uid <- .make_uid()
    backend <- getOption("giotto.disk_matrix_format", "h5")
    dirpath <- gdir[backend]
    savepath <- file.path(dirpath, uid)
    if (!dir.exists(dirpath)) dir.create(dirpath, recursive = TRUE)
    mem_formats <- c("matrix", "Matrix")

    # convert to DelayedArray if known mem format
    if (inherits(x, mem_formats)) {
        writer <- .delayedmatrix_seedwriter(backend)
        x <- writer(x, savepath)
    } else if (inherits(x, "DelayedArray")) {
        chihaya::saveDelayed(x, file = savepath)
        x <- chihaya::loadDelayed(savepath)
    } else {
        stop("[WriteMatrix] Unrecognized format to save:", class(x))
    }
    list(x = x, uid = uid)
}

# NOTE: writer fun MUST be in `x`, `filepath` arg ordering
.delayedmatrix_seedwriter <- function(backend = getOption("giotto.disk_matrix_format", "h5")) {
    match.arg(tolower(backend), c("h5", "tiledb"))
    switch(backend,
           "h5" = HDF5Array::writeHDF5Array,
           "tiledb" = TileDBArray::writeTileDBArray
    )
}

# schema name building ####

.perc_sep <- function(..., prefix) {
    out <- paste(..., sep = "%")
    if (!missing(prefix)) out <- paste0(prefix, ".", out)
    out
}

setMethod(".gschema_name", signature("exprObj"), function(x) {
    .perc_sep(
        prefix = class(x),
        x@spat_unit,
        x@feat_type,
        x@name
    )
})
