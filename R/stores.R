# store - a file location and/or save spec. Not exactly the data format itself
# data - values to read or write to the store. May be another store

# generics ####

# storeRead should return values or a representation to use
setGeneric("storeRead", function(store, ...) standardGeneric("storeRead"))
# storeWrite should not return anything
setGeneric("storeWrite", function(store, data, ...) standardGeneric("storeWrite"))

# classes ####
setOldClass("data.table")
setOldClass("data.frame")
setOldClass("matrix")
setClassUnion("memoryMatrixStore", members = c("matrix", "Matrix"))
setClassUnion("memoryStore", members = c("data.table", "data.frame", "memoryMatrixStore"))

# storeRead methods ####
setMethod("storeRead", signature("ANY"), function(store, ...) {
    if (nargs() == 1L) return(store)
    store[...]
})


# storeWrite methods ####

setMethod("storeWrite", signature("giotto", "ANY"), function(
        store, data, ...
    ) {
    gdir <- store@source
    # skip if no gDirSource
    if (identical(gdir, "")) return()
    # skip if not gDirSource -- may change for other backends
    storeWrite(gdir, data, ...)
})

setMethod("storeWrite", signature("gDirSource", "memoryMatrixStore"), function(
        store, data,
        remote_name = .make_uid(),
        backend = getOption("giotto.disk_matrix_format", "h5"),
        ...
    ) {
    dirpath <- store[backend]
    savepath <- file.path(dirpath, remote_name)
    if (!dir.exists(dirpath)) dir.create(dirpath, recursive = TRUE)
    storeWrite(savepath, data, backend = backend, ...)
})

setMethod("storeWrite", signature("gDirSource", "DelayedArray"), function(
        store, data, realize = FALSE, ...) {
    # no need to write again. Just pass through
    if (!realize) return(data)
    mem_method <- getMethod("storeWrite", c("gDirSource", "memoryMatrixStore"))
    mem_method(store, data, ...)
})

setMethod("storeWrite", signature("gDirSource", "ANY"), function(
        store, data, ...) {
    stop("[storeWrite] Unrecognized format to save:", class(data))
})

setMethod("storeWrite", signature("character", "memoryMatrixStore"), function(
        store, data,
        backend = getOption("giotto.disk_matrix_format", "h5"),
        ...
    ) {
    savepath <- store # for character stores, it is treated as the savepath
    vmsg(.is_debug = TRUE, c(.timestamp(), "writing to", savepath))
    writer <- .delayedmatrix_seedwriter(backend)
    writer(data, savepath)
})
