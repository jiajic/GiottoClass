#' @include methods-view.R
NULL

# federatedReadHandle: a lazy wrapper around per-sample substores
#
# Holds a list of per-sample data fragments (matrices, data.tables, exprObjs,
# cellMetaObjs, arrow queries, parquet substores, ...) alongside the
# information needed to fold them into a single output object on demand.
#
# Motivation: gmulti joint federation today eagerly cbinds / rbinds per-child
# fragments into a unified joint subobject. For in-memory fixtures that's
# fine, but at atlas scale the children may be file-backed parquet
# expression stores — concatenating them into a dense matrix is expensive,
# and downstream consumers (duckdb, sedonadb, arrow) can lower the
# federation into a single query plan more efficiently than we can.
#
# This class defers the materialization decision to the consumer. Holding
# the federated handle is cheap (just a list of references); calling
# `materialize()` on it (or `as.*` coercions) produces the requested
# concrete object via `@combine`.
#
# Phase 5 scope: structural class + constructor + show + length + names +
# materialize. Wiring into the federation helpers is opt-in; existing
# eager paths keep working unchanged. Future work (duckdb / sedonadb
# lowering of expression / cmeta queries) reaches for this class as the
# canonical handle shape.
#
# See `vignettes/DESIGN_gmulti_federation.md` §6 for the design.


#' @title S4 federatedReadHandle
#' @name federatedReadHandle-class
#' @description
#' Lazy wrapper around per-sample data fragments. Holds the list of
#' fragments plus a `combine` function that knows how to fold them into
#' a single output object of the requested type.
#'
#' @slot substores `list` — per-sample fragments (one entry per sample).
#'   Fragments can be any class the `@combine` function knows how to
#'   consume: matrices, data.tables, exprObjs, arrow queries, parquet
#'   substores, etc.
#' @slot keys `character` — sample names. Same length as `@substores`;
#'   names of `@substores` may be NULL for unnamed lists, in which case
#'   `@keys` is the authoritative ordering.
#' @slot output_class `character(1)` — hint for the default
#'   materialization shape (`"matrix"`, `"data.table"`, `"exprObj"`,
#'   `"cellMetaObj"`, ...). The consumer may override via the
#'   `materialize()` `as = ` argument.
#' @slot combine `function` — folds `@substores` into a single object.
#'   Signature: `function(substores, keys, ...)`. Receives the fragment
#'   list and the key vector; returns the merged object.
#' @slot meta `list` — opaque metadata for the combine function (e.g.
#'   feature intersection mode for expression matrices, column overlap
#'   rules for cmeta data.tables). Combine functions consume this via
#'   `attr(combine_result, "meta")` or by inspecting the slot directly
#'   when called through `materialize()`.
#'
#' @returns a `federatedReadHandle`
#' @exportClass federatedReadHandle
setClass("federatedReadHandle",
    slots = c(
        substores    = "list",
        keys         = "character",
        output_class = "character",
        combine      = "function",
        meta         = "list"
    ),
    prototype = list(
        substores    = list(),
        keys         = character(),
        output_class = NA_character_,
        combine      = function(substores, keys, ...) {
            stop("[federatedReadHandle] @combine not set", call. = FALSE)
        },
        meta         = list()
    )
)


#' @title Construct a federatedReadHandle
#' @name federatedReadHandle
#' @description Wrap a per-sample list of fragments + a folder function
#' into a federatedReadHandle. The handle is cheap to hold; consumers
#' call [materialize()] (or `[`/`[[`-extraction) when they need a
#' concrete object.
#'
#' @param substores per-sample list of data fragments
#' @param keys character vector of sample names (defaults to
#'   `names(substores)`)
#' @param output_class default output class hint (one of `"matrix"`,
#'   `"data.table"`, `"exprObj"`, `"cellMetaObj"`, ...). `NA` means the
#'   consumer must request via `materialize(handle, as = ...)`.
#' @param combine `function(substores, keys, ...)` that folds
#'   `@substores` into one object
#' @param meta optional `list` of metadata for the combine function
#' @returns a `federatedReadHandle`
#' @export
federatedReadHandle <- function(substores, keys = names(substores),
        output_class = NA_character_, combine, meta = list()) {
    if (missing(substores)) {
        stop("[federatedReadHandle] `substores` is required", call. = FALSE)
    }
    checkmate::assert_list(substores)
    if (is.null(keys)) keys <- character()
    checkmate::assert_character(keys, len = length(substores))
    checkmate::assert_string(output_class, na.ok = TRUE)
    if (missing(combine)) {
        stop("[federatedReadHandle] `combine` is required", call. = FALSE)
    }
    checkmate::assert_function(combine)
    checkmate::assert_list(meta)
    new("federatedReadHandle",
        substores = substores,
        keys = keys,
        output_class = output_class,
        combine = combine,
        meta = meta
    )
}


#' @noRd
setMethod("length", "federatedReadHandle",
    function(x) length(x@substores))

#' @noRd
setMethod("names", "federatedReadHandle",
    function(x) x@keys)

#' @noRd
setMethod("[[", signature(x = "federatedReadHandle", i = "ANY", j = "missing"),
    function(x, i, j, ...) {
        if (is.character(i)) {
            idx <- match(i, x@keys)
            if (is.na(idx)) {
                stop(sprintf("[federatedReadHandle] key '%s' not in handle",
                    i), call. = FALSE)
            }
            return(x@substores[[idx]])
        }
        x@substores[[i]]
    }
)


#' @noRd
setMethod("show", "federatedReadHandle", function(object) {
    cat(sprintf("federatedReadHandle (%d substore%s)\n",
        length(object),
        if (length(object) == 1L) "" else "s"))
    if (length(object) > 0L) {
        cls <- vapply(object@substores,
            function(s) class(s)[[1L]], character(1L))
        df <- data.frame(key = object@keys, class = cls,
            stringsAsFactors = FALSE)
        print(df, row.names = FALSE)
    }
    if (!is.na(object@output_class)) {
        cat(sprintf("  default output class: %s\n", object@output_class))
    }
    if (length(object@meta) > 0L) {
        cat(sprintf("  meta: %s\n",
            paste(names(object@meta), collapse = ", ")))
    }
    invisible(NULL)
})


# materialize() — fold the handle into a concrete object ####

#' @rdname materialize
#' @param as optional `character(1)` output class override. When `NULL`
#'   (default), uses the handle's `@output_class` hint. Passed through
#'   to the handle's combine function as `as`.
#' @export
setMethod("materialize",
    signature(gobject = "federatedReadHandle", view = "missing"),
    function(gobject, view, as = NULL, ...) {
        target <- as %||% gobject@output_class
        if (!is.na(target)) {
            checkmate::assert_string(target)
        }
        gobject@combine(
            substores = gobject@substores,
            keys = gobject@keys,
            as = if (is.na(target)) NULL else target,
            meta = gobject@meta,
            ...
        )
    }
)
