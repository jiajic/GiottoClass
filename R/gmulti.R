

setClass(
    "giottoMulti",
    slots = c(
        objects = "environment",  # environment of individual Giotto objects
        id_map = "list",   # Mapping of original IDs to new unified IDs
        active = "character", # Index of currently active Giotto object
        shared_data = "list" # Any data shared across all objects
    ),
    contains = "giotto"
)

setMethod("initialize", signature("giottoMulti"), function(.Object, objects, ...) {
    
    if (!missing(objects)) {
        checkmate::assert_list(objects)
        obj_names <- names(objects)
        if (is.null(obj_names)) {
            stop("objects must be a named list of giotto objects")
        }
        # append object(s) to environment
        for (name in obj_names) {
            .Object@objects[[name]] <- objects[[name]]
        }
    }
    
    # set all active if none specified
    if (is.na(.Object@active)) .Object@active <- names(.Object@objects)
    
    
})
