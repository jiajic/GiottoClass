#' @include classes-view.R
#' @include classes.R
#' @include generics.R
NULL

# =============================================================================
# methods-view.R — public API for giottoView
#
# Views are subset/narrowing recipes only — no transforms. Transforms live on
# `giottoSpace` (see methods-space.R). A view may reference a slotted space
# via `@space`; crop extents are then meaningful in that frame.
#
# Methods added here (record-style on view receiver):
#   subset(view, predicate)        NSE
#   crop(view, extent)             record viewCrop
#   selectSamples(view, ...)       record viewSampleSelect (gmulti-only)
#
# `+` for `giottoView + giottoView` is reserved (errors for now).
# =============================================================================


# Constructor ####
# `giottoView()` and `giottoView(space = "x")` dispatch through the generic
# on `signature(gobject = "missing", name = "missing")` (defined below
# alongside the accessor methods). `...` carries the optional `space`
# argument through to the method.


# Internal helper ####
.view_record_step <- function(view, step) {
    view@steps <- c(view@steps, list(step))
    view
}

# Substitute env-resident scalar / vector values into `pred` so the
# predicate becomes self-contained. Functions and missing names are left
# alone — they're resolved at eval time via the standard scope chain.
.eager_substitute_env <- function(pred, env) {
    all_vars <- all.vars(pred)
    sub_list <- list()
    for (v in all_vars) {
        if (exists(v, envir = env, inherits = TRUE)) {
            val <- tryCatch(get(v, envir = env, inherits = TRUE),
                error = function(e) NULL)
            if (!is.null(val) && !is.function(val)) {
                sub_list[[v]] <- val
            }
        }
    }
    if (length(sub_list) == 0L) return(pred)
    do.call("substitute", list(pred, sub_list))
}


# subset() — NSE predicate capture on views ####

#' @title Subset (filter) a giottoView
#' @name subset-view
#' @description
#' Record a predicate-style filter step on a [giottoView-class]. The
#' predicate expression is captured unevaluated (NSE) at call site; at
#' resolution time it is evaluated against `spatValues(g, feats = <names
#' in predicate>)`, with any additional arguments (`spat_unit`, `feat_type`,
#' `negate`, ...) forwarded through.
#'
#' Cell-keyed propagation to other slots is automatic via cell_ID relations
#' — no flag needed.
#'
#' @param x a `giottoView`
#' @param subset an unevaluated predicate expression (captured via NSE)
#' @param ... additional arguments forwarded to the underlying subset /
#'   spatValues call at resolution time (`spat_unit`, `feat_type`, `negate`,
#'   `feat_ids`, `cell_ids`, ...)
#' @returns the view, with the filter step recorded
#' @examples
#' v <- subset(giottoView(), cluster == "A")
#' v
NULL

#' @rdname subset-view
#' @export
setMethod("subset", signature(x = "giottoView"),
    function(x, subset, ...) {
        pred <- substitute(subset)
        # Walk the call stack to find the user-level frame that contains
        # the predicate's free variables. S4 dispatch + pipe + testthat
        # wrappers each insert extra frames; `parent.frame()` alone gets
        # the immediate dispatch frame, which usually has no user locals.
        env <- parent.frame()
        for (i in seq_len(8L)) {
            f <- tryCatch(parent.frame(i), error = function(e) NULL)
            if (is.null(f)) break
            if (any(vapply(all.vars(pred),
                function(v) exists(v, envir = f, inherits = FALSE),
                logical(1L)))) { env <- f; break }
        }
        # Eagerly substitute env-resident scalar / vector values into the
        # predicate so the recipe is self-contained — independent of the
        # caller's local bindings, and immune to subsequent mutation of
        # the captured vars. Keep `env` as the user's frame (not baseenv)
        # so function lookups (`median`, `c`, `%in%`, etc.) resolve via
        # the proper package chain at evaluation time.
        pred <- .eager_substitute_env(pred, env)
        .view_record_step(x, new("viewFilter",
            predicate  = pred,
            env        = env,
            scope_args = list(...)
        ))
    }
)


# crop() — extent crop on views ####

#' @rdname crop
#' @export
setMethod("crop", signature(x = "giottoView", y = "ANY"),
    function(x, y, ...) {
        .view_record_step(x, new("viewCrop", extent = y))
    }
)


# selectSamples() — gmulti-only sample selector ####

#' @title Select samples within a gmulti-scoped view
#' @name selectSamples
#' @description
#' Record a sample-selection step on a [giottoView-class] that will be
#' consumed against a [giottoMulti-class]. Picks which children participate
#' in resolution. The step is resolved FIRST, before any other step. Errors
#' at resolution time if the parent is not a `giottoMulti`.
#'
#' @param x a `giottoView`
#' @param ... `character` child names (or a single `character` vector)
#' @returns the view, with the sample-select step recorded
#' @examples
#' v <- selectSamples(giottoView(), "sample1", "sample2")
#' v
#' @export
setGeneric("selectSamples",
    function(x, ...) standardGeneric("selectSamples"))

#' @rdname selectSamples
#' @export
setMethod("selectSamples", signature(x = "giottoView"),
    function(x, ...) {
        samples <- unlist(list(...), use.names = FALSE)
        checkmate::assert_character(samples, min.len = 1L, any.missing = FALSE)
        .view_record_step(x, new("viewSampleSelect", samples = samples))
    }
)


# Composition (+) ####
# Cross-view composition is reserved for a follow-up.

#' @noRd
setMethod("+", signature(e1 = "giottoView", e2 = "giottoView"),
    function(e1, e2) {
        stop("giottoView + giottoView: cross-view composition is not yet ",
            "implemented.", call. = FALSE)
    }
)


# Accessors ####

#' @title Slotted views on a giotto object
#' @name giottoView
#' @description
#' List, retrieve, attach, or remove [giottoView-class] objects slotted into
#' a [giotto-class] object's `@view` slot.
#'
#' * `giottoView()` — construct an empty standalone view
#' * `giottoView(g, "name")` — retrieve a slotted view by name
#' * `giottoView(g, "name") <- v` — slot in (or replace) a view
#' * `giottoView(g, "name") <- NULL` — remove a view
#' * `giottoViews(g)` — list slotted view names
#'
#' Views are subset/narrowing recipes; for coordinate-frame recipes see
#' [giottoSpace-class].
#'
#' @param gobject a `giotto` object (or omitted for the constructor)
#' @param name `character(1)`. The slot key.
#' @param value a `giottoView`, or `NULL` to remove.
#' @returns the view, an updated gobject, or a character vector of view names
#' @examples
#' g <- giotto()
#' giottoView(g, "demo") <- giottoView()
#' giottoViews(g)
NULL

#' @rdname giottoView
#' @export
setGeneric("giottoView",
    function(gobject, name, ...) standardGeneric("giottoView"))

#' @rdname giottoView
#' @export
setGeneric("giottoView<-",
    function(gobject, name, ..., value) standardGeneric("giottoView<-"))

#' @rdname giottoView
#' @export
setGeneric("giottoViews",
    function(gobject, ...) standardGeneric("giottoViews"))

#' @rdname giottoView-class
#' @export
setMethod("giottoView", signature(gobject = "missing", name = "missing"),
    function(gobject, name, space = NA_character_, ...) {
        if (!is.na(space)) {
            checkmate::assert_character(space, len = 1L, any.missing = FALSE)
        }
        new("giottoView", space = as.character(space))
    }
)

#' @rdname giottoView
#' @export
setMethod("giottoView", signature(gobject = "gAny", name = "character"),
    function(gobject, name, ...) {
        checkmate::assert_character(name, len = 1L)
        v <- gobject@view[[name]]
        if (is.null(v)) {
            stop("no slotted giottoView named '", name, "'. ",
                "Available: ", paste(giottoViews(gobject), collapse = ", "),
                call. = FALSE)
        }
        v
    }
)

#' @rdname giottoView
#' @export
setMethod("giottoView", signature(gobject = "gAny", name = "missing"),
    function(gobject, name, ...) {
        nm <- giottoViews(gobject)
        if (length(nm) == 0L) return(NULL)
        if (length(nm) == 1L) return(gobject@view[[nm]])
        stop("multiple views slotted; specify `name`. ",
            "Available: ", paste(nm, collapse = ", "), call. = FALSE)
    }
)

#' @rdname giottoView
#' @export
setMethod("giottoView<-",
    signature(gobject = "gAny", name = "character", value = "giottoView"),
    function(gobject, name, ..., value) {
        checkmate::assert_character(name, len = 1L)
        value@name <- name
        if (is.null(gobject@view)) gobject@view <- list()
        gobject@view[[name]] <- value
        gobject
    }
)

#' @rdname giottoView
#' @export
setMethod("giottoView<-",
    signature(gobject = "gAny", name = "character", value = "NULL"),
    function(gobject, name, ..., value) {
        if (is.null(gobject@view) || !name %in% names(gobject@view)) {
            return(gobject)
        }
        gobject@view[[name]] <- NULL
        gobject
    }
)

#' @rdname giottoView
#' @export
setMethod("giottoViews", signature(gobject = "gAny"),
    function(gobject, ...) {
        nm <- names(gobject@view)
        if (is.null(nm)) character() else nm
    }
)


# Mutation escape hatches (stubs) ####

#' @title materialize a giottoView into a new gobject
#' @name materialize
#' @description
#' Resolve a [giottoView-class] (optionally with a slotted [giottoSpace-class]
#' frame) against a gobject and return a new gobject containing the projected
#' subobjects. Use this when downstream work needs to produce structured
#' outputs (spatial networks, dim reductions) on top of the projected data —
#' those outputs live in the materialised gobject, never in the parent.
#'
#' Read-only contract: the input gobject is not mutated.
#'
#' \strong{Status:} stub. The actual resolution engine lands in a follow-up.
#'
#' @param gobject a `giotto` object
#' @param view either a `giottoView` or a `character(1)` slot key
#' @param space `character(1)` optional — name of a slotted `giottoSpace` to
#'   resolve in. If `NULL`, uses the view's own `@space` reference (which
#'   may itself be `NA`).
#' @param coordinator a [viewCoordinator-class]-inheriting object brokering
#'   IDs and joins between storage backings. Defaults to the coordinator
#'   selected from `gobject@source` (in-memory for non-disk gobjects).
#' @param ... reserved
#' @returns a new `giotto` object reflecting the resolved view
#' @export
setGeneric("materialize",
    function(gobject, view, ...) standardGeneric("materialize"))

#' @rdname materialize
#' @export
setMethod("materialize",
    signature(gobject = "giotto", view = "giottoView"),
    function(gobject, view, space = NULL, coordinator = NULL, ...) {
        if (is.null(coordinator)) {
            coordinator <- .default_view_coordinator(gobject)
        }
        # Normalise space to a giottoSpace (or NULL) once at the entry
        # point so per-subobject resolution doesn't re-look-up by name
        space_obj <- .resolve_view_space(gobject, view, space)
        # Per-call cache shared across all slot walks within this
        # materialize. surviving_cell_ids computed at most once per call.
        cache <- .new_resolver_cache()

        out <- gobject

        # Tabular slots — narrow by surviving cell_IDs only
        out <- .materialize_walk(out, "cell_metadata",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "expression",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "dimension_reduction",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "spatial_enrichment",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "feat_metadata",
            view, space_obj, coordinator, cache)

        # Spatial slots — narrow + transform + crop
        out <- .materialize_walk(out, "spatial_locs",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "spatial_info",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "feat_info",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "images",
            view, space_obj, coordinator, cache)

        # Networks (spatial_network, nn_network) intentionally not walked:
        # they're built from a particular cell state and don't carry
        # spatial coords; view/space resolution would be misleading.

        out
    }
)

#' @rdname materialize
#' @export
setMethod("materialize",
    signature(gobject = "giotto", view = "character"),
    function(gobject, view, space = NULL, coordinator = NULL, ...) {
        v <- giottoView(gobject, view)
        materialize(gobject, v, space = space, coordinator = coordinator, ...)
    }
)


# materialize on giottoMulti ####
# 1. Apply selectSamples FIRST — narrow children before any per-child
#    work touches storage (matters at 4B-points-per-multi scale).
# 2. Per-surviving-child materialize with the child-scoped giottoSpace.
# 3. Narrow joint shared slots (multi-level @cell_metadata, @expression,
#    @dimension_reduction, @spatial_enrichment, @feat_metadata) via the
#    existing resolveSubobject dispatch — spatValues works on multi now,
#    so the joint-level predicates resolve against joint slots and the
#    surviving global cell_IDs narrow each joint subobject.

#' @rdname materialize
#' @export
setMethod("materialize",
    signature(gobject = "giottoMulti", view = "giottoView"),
    function(gobject, view, space = NULL, coordinator = NULL, ...) {
        if (is.null(coordinator)) {
            coordinator <- .default_view_coordinator(gobject)
        }
        space_obj <- .resolve_view_space(gobject, view, space)
        cache <- .new_resolver_cache()

        # Resolve selectSamples FIRST — narrow children before any
        # per-child work touches storage.
        selected <- .resolve_sample_select(gobject, view)
        if (length(selected) == 1L && is.na(selected)) {
            selected <- names(gobject@objects)
        } else {
            selected <- intersect(selected, names(gobject@objects))
        }

        out <- gobject
        out@objects <- gobject@objects[selected]

        # Per-surviving-child materialize with the child-scoped space.
        out@objects <- setNames(lapply(selected, function(samp) {
            child <- out@objects[[samp]]
            child_space <- .scope_space_to_sample(space_obj, samp)
            materialize(child, view, space = child_space,
                coordinator = coordinator, ...)
        }), selected)

        # Narrow joint shared slots — uses the same resolveSubobject
        # dispatch as the giotto path; spatValues-on-multi resolves
        # predicates against joint slots and returns global cell_IDs
        # which then filter each joint subobject's metaDT / matrix.
        out <- .materialize_walk(out, "cell_metadata",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "expression",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "dimension_reduction",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "spatial_enrichment",
            view, space_obj, coordinator, cache)
        out <- .materialize_walk(out, "feat_metadata",
            view, space_obj, coordinator, cache)

        out
    }
)

#' @rdname materialize
#' @export
setMethod("materialize",
    signature(gobject = "giottoMulti", view = "character"),
    function(gobject, view, space = NULL, coordinator = NULL, ...) {
        v <- giottoView(gobject, view)
        materialize(gobject, v, space = space, coordinator = coordinator, ...)
    }
)

# Walk one slot list (potentially nested by spat_unit / feat_type) calling
# resolveSubobject on each subobject. The slot is a `nullOrList`; structure
# is recursive — list of lists of subobjects. Apply the resolver leaf-wise.
# `cache` (optional env from .new_resolver_cache) memoises surviving_cell_ids
# across all subobjects walked within one materialize call.
#' @keywords internal
#' @noRd
.materialize_walk <- function(gobject, slot_name, view, space, coordinator,
                              cache = NULL) {
    x <- methods::slot(gobject, slot_name)
    if (is.null(x) || length(x) == 0L) return(gobject)
    methods::slot(gobject, slot_name) <- .materialize_apply(
        x, gobject, view, space, coordinator, cache)
    gobject
}

.materialize_apply <- function(node, gobject, view, space, coordinator,
                               cache = NULL) {
    if (is.list(node) && !isS4(node)) {
        return(lapply(node, .materialize_apply, gobject = gobject,
            view = view, space = space, coordinator = coordinator,
            cache = cache))
    }
    if (isS4(node) && inherits(node, "giottoSubobject")) {
        return(resolveSubobject(node, gobject, view, space, coordinator,
            .cache = cache))
    }
    node
}


#' @title Attach a view-scoped derived value back to the parent gobject
#' @name attach_derived
#' @description
#' Reattach a column-style derivation (cluster assignments, module scores,
#' QC metrics) computed under a view back to the parent gobject. The value is
#' stored alongside the existing metadata, tagged with the view name so
#' provenance is explicit.
#'
#' \strong{Status:} stub.
#'
#' @param gobject a `giotto` object
#' @param value the derived value
#' @param name `character(1)`. Column / artefact name.
#' @param view `character(1)`. The slotted view name.
#' @param ... reserved
#' @returns the gobject with the derivation attached
#' @export
setGeneric("attach_derived",
    function(gobject, value, name, view, ...) standardGeneric("attach_derived"))

#' @rdname attach_derived
#' @export
setMethod("attach_derived", signature(gobject = "giotto"),
    function(gobject, value, name, view, ...) {
        stop("attach_derived(): not yet implemented.", call. = FALSE)
    }
)


# Show methods ####

#' @noRd
setMethod("show", signature("giottoView"), function(object) {
    cat("<giottoView>\n")
    nm <- if (is.na(object@name)) "<ephemeral>" else object@name
    cat(sprintf("  name : %s\n", nm))
    sp <- if (is.na(object@space)) "<native>" else object@space
    cat(sprintf("  space: %s\n", sp))
    if (length(object@steps) == 0L) {
        cat("  (empty — pipe through `subset()`, `crop()`, `selectSamples()`)\n")
        return(invisible(NULL))
    }
    cat("  steps:\n")
    for (s in object@steps) {
        cat(sprintf("    - %s\n", .view_step_label(s)))
    }
})

#' @noRd
setMethod("show", signature("viewStep"), function(object) {
    cat(sprintf("<%s> %s\n", class(object)[[1L]], .view_step_label(object)))
})

# Compact one-line label for a viewStep — used by show methods.
.view_step_label <- function(step) {
    if (inherits(step, "viewFilter")) {
        scope <- if (length(step@scope_args) > 0L) {
            scope_str <- paste(names(step@scope_args), "=",
                vapply(step@scope_args, function(a) {
                    tryCatch(deparse(a, nlines = 1L)[[1L]],
                        error = function(e) "?")
                }, character(1L)), collapse = ", ")
            sprintf(" [%s]", scope_str)
        } else ""
        return(sprintf("filter:    %s%s",
            tryCatch(deparse(step@predicate, nlines = 1L)[[1L]],
                error = function(e) "<predicate>"),
            scope))
    }
    if (inherits(step, "viewCrop")) {
        return(sprintf("crop:      %s",
            tryCatch(deparse(step@extent, nlines = 1L)[[1L]],
                error = function(e) "<extent>")))
    }
    if (inherits(step, "viewSampleSelect")) {
        return(sprintf("samples:   %s",
            paste(step@samples, collapse = ", ")))
    }
    sprintf("<%s>", class(step)[[1L]])
}
