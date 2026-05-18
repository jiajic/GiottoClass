# =============================================================================
# giottoMulti — multi-dataset container for shared expression-space analysis
# =============================================================================
#
# DESIGN NOTES (sketch, 2026-05-18)
# ---------------------------------
# `giottoMulti` represents N giotto objects that each keep their own SPATIAL
# information (different physical/embedding spaces) but participate in a SHARED
# expression-space analysis (joint normalization, dim reduction, NN graphs,
# clustering, integration).
#
# Inheritance:
#   gAny (virtual)
#   ├── giotto
#   └── giottoMulti
#
# Why a virtual base rather than `giottoMulti contains "giotto"`:
#   * `giottoMulti`'s spatial slots are intentionally empty — they live in
#     `@objects`. Inheriting from `giotto` would let spatial-domain methods
#     silently fall through to those empty slots. With `gAny` as a virtual
#     base, an undefined method fails loudly via no-method dispatch.
#   * Shared-domain methods are written once on `gAny` and apply to both.
#   * `giotto`'s validity rules don't have to accommodate empty spatial slots.
#
# Slot strategy: ALIGN where semantics are shared, ADD where new.
#   * Aligned with `giotto` (same name, same type): @expression,
#     @cell_metadata, @feat_metadata, @dimension_reduction, @nn_network,
#     @spatial_enrichment, @multiomics, @instructions, @parameters, etc.
#     Existing accessor generics (`getExpression`, `pDataDT`, ...) can be
#     promoted from `"giotto"` to `"gAny"` one at a time.
#   * NEW: @objects, @id_map, @active, @access — multi-specific.
#   * OMITTED from `giotto`: @spatial_locs, @spatial_info, @spatial_network,
#     @spatial_grid, @feat_info, @images, @join_info, @offset_file, @h5_file.
#     These are per-dataset and live in the children.
#
# Children in @objects:
#   * Named list of `giotto` objects (or disk-backed equivalents).
#   * `list` chosen over `environment` to keep value semantics — disk-backed
#     children get reference semantics from the storage layer regardless;
#     avoiding R-level reference semantics on the container removes a class
#     of surprise.
#
# id_map structure:
#   * `@id_map$cells`: data.table(object, local_id, global_id) — namespaces
#     cell IDs across children. Default global_id = paste(object, local_id, sep = "::").
#   * `@id_map$feats`: data.table(object, local_id, global_id) — features
#     usually overlap legitimately across datasets, so global_id often equals
#     local_id. Used to track which features are present in which datasets.
#
# active / access:
#   * `@active`: character vector of object names to route to by default for
#     per-object operations. `NA_character_` ≙ all objects.
#   * `@access`: data.frame(object, spat_unit, feat_type) giving the default
#     spat_unit/feat_type per child. Lets methods on `giottoMulti` resolve
#     "which spat_unit do I read from for object X" without forcing the
#     caller to specify per-child.
#
# Dispatch pattern (recommended, applied incrementally):
#   * Shared-domain method:
#       setMethod("getExpression", "gAny", function(gobject, ...) { ... })
#   * Spatial-domain / per-object method:
#       setMethod("getSpatialLocations", "giottoMulti",
#           function(gobject, object = NULL, ...) {
#               objs <- .gm_resolve_active(gobject, object)
#               lapply(objs, function(o) getSpatialLocations(gobject[[o]], ...))
#           })
#
# Open design questions (NOT resolved in this sketch):
#   * How does joinGiottoObjects interact with giottoMulti? Most natural
#     read: `giottoMulti(list(g1, g2, ...))` is the "preserve spaces"
#     constructor; `joinGiottoObjects()` remains the "merge into one
#     giotto" path. They are siblings, not nested.
#   * How does `subset()` work on giottoMulti? Subset children, subset cells
#     globally (with id_map lookup), or both?
#   * Do per-child cell-metadata columns get reflected upward into the joint
#     @cell_metadata, or kept separate and merged on demand?
#   * Where do integration-method parameters (Harmony, Seurat-anchor) live?
#     Probably @parameters, possibly a dedicated @integration slot if it
#     grows. Defer until we have a concrete use case.
# =============================================================================


# CLASS ####

#' @title S4 giottoMulti
#' @name giottoMulti-class
#' @description
#' Container for multiple `giotto` objects whose spatial information is kept
#' separate (one space per child) but whose expression-space analysis is
#' shared across all cells.
#'
#' @slot objects named `list` of `giotto` objects (children)
#' @slot id_map `list` with elements `cells` and `feats`, each a `data.table`
#'   mapping `(object, local_id) → global_id`
#' @slot active `character` vector of currently-active object names. Defaults
#'   to all object names.
#' @slot access `data.frame` of per-object default `spat_unit` and `feat_type`
#'
#' @slot expression shared expression matrices (rows = union of features,
#'   cols = global cell IDs)
#' @slot cell_metadata shared cell metadata (one row per global cell ID)
#' @slot feat_metadata shared feature metadata (one row per global feature)
#' @slot cell_ID shared cell ID lists (global IDs)
#' @slot feat_ID shared feature ID lists (global IDs)
#' @slot dimension_reduction shared joint dim-reductions (PCA, UMAP, harmony)
#' @slot nn_network shared joint NN graphs
#' @slot spatial_enrichment shared spatial enrichment results
#' @slot multiomics shared multi-omics info
#'
#' @slot instructions giotto-style instructions
#' @slot parameters analysis parameters (mirrors `giotto@parameters`)
#' @slot versions package versions
#' @slot misc miscellaneous
#'
#' @returns giottoMulti object
#' @exportClass giottoMulti
giottoMulti <- setClass(
    "giottoMulti",
    contains = "gAny",
    slots = c(
        # multi-specific
        objects             = "list",
        id_map              = "list",
        active              = "character",
        access              = "data.frame",

        # shared-domain (names aligned with giotto)
        expression          = "nullOrList",
        expression_feat     = "nullOrChar",
        cell_metadata       = "nullOrList",
        feat_metadata       = "nullOrList",
        cell_ID             = "nullOrList",
        feat_ID             = "nullOrList",
        spatial_enrichment  = "nullOrList",
        dimension_reduction = "nullOrList",
        nn_network          = "nullOrList",
        multiomics          = "ANY",

        # infrastructure
        instructions        = "nullOrInstructions",
        parameters          = "ANY",
        versions            = "list",
        misc                = "list"
    ),
    prototype = list(
        objects             = list(),
        id_map              = list(cells = NULL, feats = NULL),
        active              = NA_character_,
        access              = data.frame(
            object = character(),
            spat_unit = character(),
            feat_type = character(),
            stringsAsFactors = FALSE
        ),

        expression          = NULL,
        expression_feat     = NULL,
        cell_metadata       = NULL,
        feat_metadata       = NULL,
        cell_ID             = NULL,
        feat_ID             = NULL,
        spatial_enrichment  = NULL,
        dimension_reduction = NULL,
        nn_network          = NULL,
        multiomics          = NULL,

        instructions        = NULL,
        parameters          = list(),
        versions            = .versions_info(),
        misc                = list()
    )
)


# INITIALIZE ####

#' @noRd
setMethod("initialize", signature("giottoMulti"), function(.Object, objects = NULL, ...) {
    .Object <- callNextMethod(.Object, ...)

    if (is.null(objects) || length(objects) == 0L) return(.Object)

    checkmate::assert_list(objects, types = "giotto", names = "unique",
        .var.name = "objects")

    .Object@objects <- objects

    # active defaults to all
    if (length(.Object@active) == 1L && is.na(.Object@active)) {
        .Object@active <- names(objects)
    }

    # access table: one row per child with its default spat_unit / feat_type
    if (nrow(.Object@access) == 0L) {
        .Object@access <- .gm_default_access(objects)
    }

    # id_map: namespace cell IDs as "{object}::{local_id}"; feats default to
    # passthrough (overlap across datasets is real overlap)
    if (is.null(.Object@id_map$cells)) {
        .Object@id_map$cells <- .gm_build_cell_idmap(objects)
    }
    if (is.null(.Object@id_map$feats)) {
        .Object@id_map$feats <- .gm_build_feat_idmap(objects)
    }

    .Object
})


# CONSTRUCTOR ####

#' @title Create a giottoMulti object
#' @name createGiottoMulti
#' @description Container for multiple `giotto` objects analyzed in a shared
#' expression space. Each child keeps its own spatial information; shared
#' analyses (joint dim reduction, NN graphs, clustering) live on the parent.
#'
#' @param objects named `list` of `giotto` objects
#' @param active `character` vector of object names to mark active. Defaults
#'   to all.
#' @param instructions a `giottoInstructions` object (optional)
#'
#' @returns `giottoMulti`
#' @examples
#' \dontrun{
#' g1 <- GiottoData::loadGiottoMini("visium")
#' g2 <- GiottoData::loadGiottoMini("viz")
#' mg <- createGiottoMulti(list(visium = g1, viz = g2))
#' }
#' @export
createGiottoMulti <- function(objects, active = NULL, instructions = NULL) {
    checkmate::assert_list(objects, types = "giotto", names = "unique")
    args <- list(objects = objects)
    if (!is.null(active)) args$active <- active
    if (!is.null(instructions)) args$instructions <- instructions
    do.call(new, c("giottoMulti", args))
}


# INTROSPECTION ####

#' @noRd
setMethod("names", "giottoMulti", function(x) names(x@objects))

#' @noRd
setMethod("length", "giottoMulti", function(x) length(x@objects))

#' @noRd
setMethod("[[", signature(x = "giottoMulti", i = "ANY", j = "missing"),
    function(x, i, j, ...) x@objects[[i]])

#' @noRd
setReplaceMethod("[[", signature(x = "giottoMulti", i = "ANY", j = "missing", value = "giotto"),
    function(x, i, j, ..., value) {
        x@objects[[i]] <- value
        # invalidate id_map for this object — caller is expected to refresh
        # via .gm_rebuild_idmap(x). Doing it eagerly would be surprising for
        # large objects; leave it explicit.
        x
    }
)


# SHOW ####

#' @noRd
setMethod("show", "giottoMulti", function(object) {
    cat(sprintf("An object of class %s\n", class(object)))
    cat(sprintf("  %d object(s): %s\n",
        length(object),
        paste(names(object), collapse = ", ")))
    if (length(object@active) > 0L && !any(is.na(object@active))) {
        cat(sprintf("  active: %s\n",
            paste(object@active, collapse = ", ")))
    }
    if (!is.null(object@id_map$cells)) {
        cat(sprintf("  %d total cells (global)\n",
            nrow(object@id_map$cells)))
    }
    if (!is.null(object@id_map$feats)) {
        cat(sprintf("  %d total features (global)\n",
            nrow(object@id_map$feats)))
    }
    invisible(NULL)
})


# INTERNAL HELPERS ####

#' @noRd
.gm_default_access <- function(objects) {
    rows <- lapply(names(objects), function(nm) {
        g <- objects[[nm]]
        data.frame(
            object = nm,
            spat_unit = tryCatch(set_default_spat_unit(g), error = function(e) NA_character_),
            feat_type = tryCatch(set_default_feat_type(g), error = function(e) NA_character_),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

#' @noRd
.gm_build_cell_idmap <- function(objects, sep = "::") {
    parts <- lapply(names(objects), function(nm) {
        ids <- tryCatch(spatIDs(objects[[nm]]), error = function(e) character())
        if (length(ids) == 0L) return(NULL)
        data.table::data.table(
            object = nm,
            local_id = ids,
            global_id = paste(nm, ids, sep = sep)
        )
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) return(NULL)
    data.table::rbindlist(parts)
}

#' @noRd
.gm_build_feat_idmap <- function(objects) {
    parts <- lapply(names(objects), function(nm) {
        ids <- tryCatch(featIDs(objects[[nm]]), error = function(e) character())
        if (length(ids) == 0L) return(NULL)
        # default passthrough: feature names are the same global vocabulary
        data.table::data.table(
            object = nm,
            local_id = ids,
            global_id = ids
        )
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) return(NULL)
    data.table::rbindlist(parts)
}

#' Resolve which children a per-object method should operate on.
#'
#' @param x giottoMulti
#' @param object NULL (use @active), or character vector of object names
#' @returns character vector of object names
#' @noRd
.gm_resolve_active <- function(x, object = NULL) {
    if (is.null(object)) {
        if (length(x@active) == 1L && is.na(x@active)) return(names(x))
        return(x@active)
    }
    checkmate::assert_character(object)
    bad <- setdiff(object, names(x))
    if (length(bad) > 0L) {
        stop("unknown object(s): ", paste(bad, collapse = ", "), call. = FALSE)
    }
    object
}


# DISPATCH PATTERN — EXAMPLES (not yet exhaustive) ####
#
# Two representative methods. The point is the shape, not coverage —
# we add more incrementally as needs arise.

# --- shared-domain example -------------------------------------------------
# `getExpression` reads from @expression on either class. Promote to gAny so
# one definition serves both. Existing `setMethod("getExpression", "giotto", ...)`
# in slot_accessors.R should be migrated when we're confident the change
# doesn't break anything spatial-aware inside that method.
#
# (NOT activating here — placeholder showing the pattern.)
#
# setMethod("getExpression", "gAny", function(gobject, ...) {
#     # implementation identical to the current "giotto" version
# })


# --- per-object example ----------------------------------------------------
# Spatial-domain accessors return per-child lists rather than try to combine.
# Caller picks which child to act on via `object =`, or operates on all
# active children.
#
# (Stub — wire up once we settle the per-child API ergonomics.)
#
# setMethod("getSpatialLocations", "giottoMulti",
#     function(gobject, object = NULL, ...) {
#         objs <- .gm_resolve_active(gobject, object)
#         out <- lapply(objs, function(nm) {
#             getSpatialLocations(gobject[[nm]], ...)
#         })
#         names(out) <- objs
#         out
#     })