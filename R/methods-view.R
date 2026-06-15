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


# Indirect-usage routing: lets generics like subset() / crop() on a
# `giotto` accept `view = <name|giottoView>` and record the step rather
# than executing eagerly. Returns:
#   - the gobject (with the named view slotted / appended) when `view`
#     is a character name
#   - the modified giottoView object when `view` is a recipe (caller
#     continues building before slotting later)
.record_view_on_gobject <- function(gobject, view, step) {
    # view contract: character(1) name of a slotted view, or NULL.
    # Inline giottoView objects were considered and rejected (see
    # vignettes/DESIGN_gmulti_federation.md). If programmatic composition
    # is needed, build the view, slot it under a name, then reference it:
    #   v <- giottoView() |> subset(...) |> crop(...)
    #   giottoView(g, "tmp") <- v
    #   subset(g, ..., view = "tmp")
    checkmate::assert_string(view, .var.name = "view")
    existing <- if (view %in% giottoViews(gobject)) {
        giottoView(gobject, view)
    } else {
        giottoView()
    }
    new_view <- .view_record_step(existing, step)
    giottoView(gobject, view) <- new_view
    gobject
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


# crop() — region narrowing on views ####
# Note: name kept as crop() for ergonomics (familiar verb), but the semantic
# is relate-based membership, not geometric clipping. The view records the
# region + relation; resolution narrows cells whose centroid satisfies the
# relation against the region. For polygon regions and non-default
# relations, terra::is.related provides the precise check; rectangular
# extents take an AABB short-circuit. Geometry of surviving subobjects is
# NOT modified — only the cell set narrows.

#' @rdname crop
#' @param relation `character(1)`. Spatial relation evaluated against cell
#'   centroids. One of `"intersects"` (default), `"within"`, `"contains"`,
#'   `"covers"`, `"covered_by"`, `"overlaps"`, `"touches"`, `"crosses"`,
#'   `"disjoint"`. Passed through to [terra::is.related].
#' @export
setMethod("crop", signature(x = "giottoView", y = "ANY"),
    function(x, y, relation = "intersects", ...) {
        checkmate::assert_character(relation, len = 1L, any.missing = FALSE)
        .view_record_step(x, new("viewCrop", region = y, relation = relation))
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
#' @param slots optional `character` vector of slot names to narrow.
#'   When `NULL` (default), all slot lists are walked (`cell_metadata`,
#'   `expression`, `dimension_reduction`, `spatial_enrichment`,
#'   `feat_metadata`, `spatial_locs`, `spatial_info`, `feat_info`,
#'   `images`). When supplied, only the listed slots are walked — the
#'   rest are left untouched on the returned object. Useful for
#'   internal helpers that only consume a subset of slots and want to
#'   share one resolver pass without paying for irrelevant slots.
#' @param ... reserved
#' @returns a new `giotto` object reflecting the resolved view
#' @export
setGeneric("materialize",
    function(gobject, view, ...) standardGeneric("materialize"))


# All slots `materialize()` knows how to walk, in the canonical order
# (tabular → spatial → images). Used as the default slot set when
# `slots = NULL` and to validate caller-supplied slot names.
.materialize_default_slots <- c(
    "cell_metadata", "expression", "dimension_reduction",
    "spatial_enrichment", "feat_metadata",
    "spatial_locs", "spatial_info", "feat_info", "images"
)


# Validate and order a caller-supplied slot vector against the
# canonical walk order. NULL → all default slots. Unknown slot names
# error.
#' @keywords internal
#' @noRd
.materialize_slot_filter <- function(slots) {
    if (is.null(slots)) return(.materialize_default_slots)
    bad <- setdiff(slots, .materialize_default_slots)
    if (length(bad) > 0L) {
        stop(sprintf(
            "[materialize] unknown slot(s): %s. Available: %s",
            paste(bad, collapse = ", "),
            paste(.materialize_default_slots, collapse = ", ")
        ), call. = FALSE)
    }
    intersect(.materialize_default_slots, slots)  # canonical order
}

# Internal implementation: materialize on a giotto with an already-resolved
# giottoView object. Called from the public character-signature method
# (after slot lookup) and from the giottoMulti per-child loop (where the
# view object is already in hand). Not user-facing; the public API is
# the character-signature method below.
#' @keywords internal
#' @noRd
.materialize_giotto_resolved <- function(gobject, view,
                                          space = NULL,
                                          coordinator = NULL,
                                          slots = NULL,
                                          ...) {
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

    # Walk the (possibly filtered) slot list in canonical order:
    # tabular → spatial → images. Slot names not in `slots` are
    # left untouched on the returned gobject.
    for (slot_name in .materialize_slot_filter(slots)) {
        out <- .materialize_walk(out, slot_name,
            view, space_obj, coordinator, cache)
    }

    # Networks (spatial_network, nn_network) intentionally not walked:
    # they're built from a particular cell state and don't carry
    # spatial coords; view/space resolution would be misleading.

    out
}

#' @rdname materialize
#' @export
setMethod("materialize",
    signature(gobject = "giotto", view = "character"),
    function(gobject, view, space = NULL, coordinator = NULL,
             slots = NULL, ...) {
        v <- giottoView(gobject, view)
        .materialize_giotto_resolved(gobject, v,
            space = space, coordinator = coordinator,
            slots = slots, ...)
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

# Internal implementation: materialize on a giottoMulti with an
# already-resolved giottoView object. Called from the public
# character-signature method (after slot lookup). Not user-facing.
#' @keywords internal
#' @noRd
.materialize_gmulti_resolved <- function(gobject, view,
                                          space = NULL,
                                          coordinator = NULL,
                                          slots = NULL,
                                          ...) {
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
    # `slots` filter is forwarded so per-child narrowing matches the
    # joint-level scope. Uses the internal resolved-view helper directly
    # — the public materialize() dispatch only accepts character views,
    # but per-child iteration already holds the giottoView object.
    out@objects <- setNames(lapply(selected, function(samp) {
        child <- out@objects[[samp]]
        child_space <- .scope_space_to_sample(space_obj, samp)
        .materialize_giotto_resolved(child, view, space = child_space,
            coordinator = coordinator, slots = slots, ...)
    }), selected)

    # Narrow joint shared slots. Joint-level walk respects the
    # `slots` filter: only multi-level cell_metadata / expression /
    # dim_reduction / spatial_enrichment / feat_metadata are
    # legitimately joint, so we intersect with that subset.
    joint_candidates <- c("cell_metadata", "expression",
        "dimension_reduction", "spatial_enrichment", "feat_metadata")
    joint_slots <- intersect(.materialize_slot_filter(slots),
        joint_candidates)
    for (slot_name in joint_slots) {
        out <- .materialize_walk(out, slot_name,
            view, space_obj, coordinator, cache)
    }

    out
}

#' @rdname materialize
#' @export
setMethod("materialize",
    signature(gobject = "giottoMulti", view = "character"),
    function(gobject, view, space = NULL, coordinator = NULL,
             slots = NULL, ...) {
        v <- giottoView(gobject, view)
        .materialize_gmulti_resolved(gobject, v,
            space = space, coordinator = coordinator,
            slots = slots, ...)
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
        rel <- if (identical(step@relation, "intersects")) ""
            else sprintf(" [%s]", step@relation)
        region_lbl <- if (inherits(step@region, "SpatVector")) {
            sprintf("<SpatVector: %d geoms>", length(step@region))
        } else {
            tryCatch(deparse(step@region, nlines = 1L)[[1L]],
                error = function(e) "<region>")
        }
        return(sprintf("crop:      %s%s", region_lbl, rel))
    }
    if (inherits(step, "viewSampleSelect")) {
        return(sprintf("samples:   %s",
            paste(step@samples, collapse = ", ")))
    }
    sprintf("<%s>", class(step)[[1L]])
}
