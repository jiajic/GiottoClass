# =============================================================================
# giottoSpace — coordinate-frame recipe (parallel space opt-in)
# =============================================================================
#
# DESIGN NOTES (sketch, 2026-05-28)
# ---------------------------------
# `giottoSpace` is a composable standalone S4 class describing a coordinate
# frame for a `giotto` (single-sample) or `giottoMulti` (multi-sample) object.
# Slotted spaces are named alternate coordinate frames that consumer functions
# opt into via `space = "name"`.
#
# Unlike views (which are read-only narrowings), spaces are NOT subject to
# the read-only contract — analyses run in a non-native space are fine; the
# coordinate frame just differs. Mutations still target the underlying data
# in its native frame.
#
# Sample scope: transforms in a `giottoSpace` are SAMPLE-UNIFORM. Within a
# single sample, all spatial elements (cells, polys, points, image,
# spatlocs) move together. Per-element overrides are deliberately not
# supported here — sample-level is the granularity that matches the typical
# spatialomics alignment workflow.
#
# Shape:
#   # single-sample (giotto)
#   s <- giottoSpace() |> affine(M)
#   giottoSpace(g, "tilted") <- s
#
#   # multi-sample (giottoMulti)
#   sa <- giottoSpace("sample_a") |> affine(M_a)
#   sb <- giottoSpace("sample_b") |> affine(M_b)
#   atlas <- sa + sb                       # combine sample recipes
#   giottoSpace(mg, "atlas") <- atlas
#
# `+` composition:
#   * same-sample (`sample_a` + `sample_a`) → steps concatenated in order
#   * different-sample → samples merged into one giottoSpace keyed by name
#
# Step taxonomy:
#   spaceTransform — records a deferred call to one of the GiottoClass
#   spatial transform generics (`affine`, `spin`, `spatShift`, `flip`,
#   `rescale`, `shear`, `zoom`). At resolution time the receiving gobject
#   (or child) is spliced as the first argument and `do.call()` dispatches
#   to the existing transform method.
#
# Storage on gobject:
#   `gobject@spaces` — named list of `giottoSpace`. Sentinel sample name
#   `:default:` is used for single-giotto entries (no explicit sample).
#
# See `R/classes-view.R` for the subset/narrowing recipe.
# See `R/methods-space.R` for the constructors, `+` composition, record
# methods on the existing transform generics, accessor, show.
# =============================================================================


# Sentinel for sample-anonymous (single-giotto) space construction.
.space_default_sample <- ":default:"


# spaceTransform step ####

#' @title spaceTransform
#' @description Captures a deferred call to one of the GiottoClass spatial
#' transform generics (`affine`, `spin`, `spatShift`, `flip`, `rescale`,
#' `shear`, `zoom`). At resolution time the receiving object is spliced as
#' the first argument and `do.call()` dispatches to the existing transform
#' method for that class.
#' @keywords internal
#' @noRd
setClass(
    "spaceTransform",
    slots = list(
        op   = "character",
        args = "list",
        misc = "list"
    ),
    prototype = list(
        op = NA_character_, args = list(), misc = list()
    )
)


# giottoSpace ####

#' @title S4 giottoSpace class
#' @name giottoSpace-class
#' @description A `giottoSpace` is a composable, opt-in coordinate-frame
#' recipe over a `giotto` (single-sample) or `giottoMulti` (multi-sample)
#' object. It records a sequence of spatial transform steps keyed by sample
#' name; consumer functions opt into a named space via `space = "name"` and
#' the recorded transforms are applied to position the data in that frame.
#'
#' Transforms in a space are SAMPLE-UNIFORM — within a single sample, all
#' spatial elements move together. Per-element overrides are not supported;
#' use [materialize()] and per-element transforms post-hoc for that.
#'
#' Compose with `+` to combine per-sample recipes into a multi-sample space:
#'
#' ```r
#' sa <- giottoSpace("sample_a") |> affine(M_a)
#' sb <- giottoSpace("sample_b") |> affine(M_b)
#' atlas <- sa + sb
#' giottoSpace(mg, "atlas") <- atlas
#' ```
#'
#' Same-sample composition concatenates steps in order:
#'
#' ```r
#' s <- (giottoSpace("sample_a") |> spin(30)) +
#'      (giottoSpace("sample_a") |> spatShift(dx = 10))
#' ```
#'
#' @slot samples named `list`. Keys are sample names (`:default:` for
#'   single-giotto context); values are lists of `spaceTransform` step
#'   objects to apply in order.
#' @slot name `character(1)`. `NA_character_` until slotted into a gobject.
#' @slot source `ANY`. Reserved pointer / fingerprint. `NULL` for standalone.
#' @slot misc `list`. Provenance, cache keys.
#' @returns `giottoSpace`
#' @examples
#' giottoSpace()
#' giottoSpace("sample_a")
#' @export
#' @exportClass giottoSpace
setClass(
    "giottoSpace",
    slots = list(
        samples = "list",
        name    = "character",
        source  = "ANY",
        misc    = "list"
    ),
    prototype = list(
        samples = list(),
        name    = NA_character_,
        source  = NULL,
        misc    = list()
    ),
    validity = function(object) {
        if (length(object@samples) > 0L) {
            nms <- names(object@samples)
            if (is.null(nms) || any(is.na(nms)) || any(nms == "")) {
                return("@samples must be a named list (sample name -> step list)")
            }
        }
        TRUE
    }
)
