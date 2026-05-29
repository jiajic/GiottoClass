# =============================================================================
# giottoView — read-only subset / narrowing recipe
# =============================================================================
#
# DESIGN NOTES (sketch, 2026-05-28)
# ---------------------------------
# `giottoView` is a composable standalone S4 class describing a deferred,
# read-only NARROWING of a `giotto` (or `giottoMulti`) object. A view is a
# *recipe*, not a snapshot: each time it is consumed it is re-resolved
# against the current state of the underlying gobject.
#
# Spatial transforms (positioning) live in [giottoSpace-class], NOT here.
# Views handle subsets, crops, and sample-selection only. The two compose at
# the consumer-function API: `plot(g, space = "atlas", view = "tumor")`.
#
# A view may optionally reference a named space via `@space`. This is what
# defines the coordinate frame in which extent-based crops are meaningful:
# `crop(c(0, 100, 0, 100))` records as a viewCrop step against the view's
# `@space` reference, so the resolver knows to position the data first.
#
# Shape:
#   v <- giottoView(space = "atlas") |>
#       subset(cluster == "A") |>             # cell-keyed predicate
#       crop(c(0, 100, 0, 100))               # crop in atlas frame
#
#   v_multi <- giottoView() |>
#       selectSamples("a", "b") |>            # gmulti-only child filter
#       subset(cluster == "A")
#
#   giottoView(g, "tumor_focus") <- v
#
# Step taxonomy (subset-flavor only — transforms live on giottoSpace):
#   viewStep                  virtual base
#   ├── viewFilter            predicate-style row filter (NSE)
#   ├── viewCrop              extent-based crop (extent meaningful in @space)
#   └── viewSampleSelect      gmulti-only child filter
#
# Cell-keyed propagation:
#   A `subset()` predicate is evaluated against `spatValues(g)` for the
#   columns it references. Surviving cell_IDs propagate to cell-keyed slots
#   (expression, spatial_locs, dim_reduction, polys with a `cell_ID` column)
#   automatically via the existing relational structure — no flag needed.
#   Polygons without a `cell_ID` linkage are out of scope for views; link
#   them first or handle them in a separate step.
#
# Read-only contract (signature-as-contract):
#   * Views never mutate the underlying gobject.
#   * Functions that accept a `view =` parameter return their result.
#   * Functions that mutate the gobject do not accept a `view =` parameter.
#   * `attach_derived()` — explicit escape hatch from view-scoped column
#     output back to the gobject.
#   * `materialize()` — explicit escape hatch from a view to a new
#     standalone gobject.
#
# See `R/classes-space.R` for the spatial-transform recipe (`giottoSpace`).
# See `R/methods-view.R` and `R/methods-space.R` for the constructors,
# composition, accessors, show methods, and stubs.
# =============================================================================


# viewStep taxonomy (subset-flavor) ####

#' @title viewStep virtual class
#' @description Base class for individual steps in a `giottoView` recipe.
#' Concrete subclasses: `viewFilter`, `viewCrop`, `viewSampleSelect`. Spatial
#' transforms live on [giottoSpace-class] via `spaceTransform`, not here.
#' @keywords internal
#' @noRd
setClass(
    "viewStep",
    contains = "VIRTUAL",
    slots = list(misc = "list"),
    prototype = list(misc = list())
)

#' @title viewFilter
#' @description Captures a predicate expression (e.g. `cluster == "A"`) to be
#' evaluated lazily against `spatValues(g, feats = <names_in_predicate>)`
#' at resolution time. `scope_args` carries through any additional arguments
#' passed alongside the predicate at `subset()` call site (e.g. `spat_unit`,
#' `feat_type`, `negate`) so they reach the underlying subset / spatValues
#' call. Cell-keyed propagation to other slots is automatic via cell_ID
#' relations; no flag needed.
#' @keywords internal
#' @noRd
setClass(
    "viewFilter",
    contains = "viewStep",
    slots = list(
        predicate  = "ANY",
        env        = "ANY",
        scope_args = "list"
    ),
    prototype = list(predicate = NULL, env = NULL, scope_args = list())
)

#' @title viewCrop
#' @description Spatial-region narrowing step. `region` is the target
#' geometry: a numeric extent vector (`c(xmin, xmax, ymin, ymax)`), a
#' `SpatExtent`, or a `SpatVector` polygon. `relation` is the membership
#' relation evaluated against cell centroids — `"intersects"` (default),
#' `"within"`, `"contains"`, `"covers"`, etc. — any relation supported by
#' [terra::is.related].
#'
#' The region is meaningful in the coordinate frame of the view's `@space`
#' reference (or the gobject's native frame if `@space` is NA). At
#' resolution time: apply the referenced space's transforms first, then
#' filter cells by the relation.
#'
#' Implementation note: for rectangular `region`s (numeric / SpatExtent),
#' an AABB short-circuit is used in-memory. For polygon `region`s, an AABB
#' pre-filter narrows candidates before the precise relate check, then
#' [terra::is.related] gives the final survival set. On disk-backed
#' subobjects, the parquet store's existing crop machinery handles
#' inverse-affine back-projection + AABB pushdown + half-plane filter
#' automatically (GiottoDisk side).
#' @keywords internal
#' @noRd
setClass(
    "viewCrop",
    contains = "viewStep",
    slots = list(
        region = "ANY",
        relation = "character"
    ),
    prototype = list(region = NULL, relation = "intersects")
)

#' @title viewSampleSelect
#' @description gmulti-only step selecting which children participate.
#' Resolved FIRST, before any other steps.
#' @keywords internal
#' @noRd
setClass(
    "viewSampleSelect",
    contains = "viewStep",
    slots = list(samples = "character"),
    prototype = list(samples = character())
)


# giottoView ####

#' @title S4 giottoView class
#' @name giottoView-class
#' @description A `giottoView` is a composable, read-only, lazy subset recipe
#' over a `giotto` (or `giottoMulti`) object. It records predicate filters,
#' extent crops, and sample selectors. Spatial transforms (positioning) live
#' on [giottoSpace-class], not here.
#'
#' At resolution time, the recorded steps are applied against the current
#' gobject state without mutating it. Slotted views (added via
#' `giottoView(g, "name") <- v`) travel with the gobject through save/load.
#'
#' Compose by piping through `subset()`, `crop()`, `selectSamples()`:
#'
#' ```r
#' v <- giottoView(space = "atlas") |>
#'     subset(cluster == "A") |>
#'     crop(c(0, 100, 0, 100))
#'
#' giottoView(g, "tumor_focus") <- v
#' ```
#'
#' @slot steps `list` of `viewStep` objects.
#' @slot space `character(1)`. Optional reference to a slotted space name
#'   (see [giottoSpace-class]). The coordinate frame in which crop extents
#'   are interpreted. `NA_character_` means the gobject's native frame.
#' @slot name `character(1)`. `NA_character_` until the view is slotted into
#'   a gobject; thereafter holds the slot key.
#' @slot source `ANY`. Reserved pointer / fingerprint. `NULL` for standalone.
#' @slot misc `list`. Provenance, cache keys, version stamps.
#' @returns `giottoView`
#' @examples
#' giottoView()
#' @export
#' @exportClass giottoView
setClass(
    "giottoView",
    slots = list(
        steps  = "list",
        space  = "character",
        name   = "character",
        source = "ANY",
        misc   = "list"
    ),
    prototype = list(
        steps  = list(),
        space  = NA_character_,
        name   = NA_character_,
        source = NULL,
        misc   = list()
    )
)
