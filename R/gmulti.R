

setClass(
    "giottoMulti",
    slots = c(
        objects = "environment",  # environment of individual Giotto objects
        id_map = "list",   # Mapping of original IDs to new unified IDs
        active = "character", # Index of currently active Giotto object
        access = "data.frame", # table of spatunit and feat type to access
        meta = "data.frame" # Any data shared across all objects
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
    
    # if empty multi, return early
    if (length(.Object@objects) == 0) return(.Object)
    
    obj_ids <- names(.Object@objects)
    
    # set all active if none specified
    if (is.na(.Object@active)) .Object@active <- obj_ids
    
    # access spat/feat
    # - append in non-registered objects
    missing_obj <- obj_ids[!obj_ids %in% .Object@access$object]
    
    
    
    
    for (g in ls(.Object@objects)) {
        spat <- set_default_spat_unit(g)
        feat <- set_default_feat_type(g, spat_unit = spat)
    }
    
    
})
