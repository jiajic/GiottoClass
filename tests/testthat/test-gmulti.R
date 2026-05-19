# Tests for giottoMulti — sketch-level coverage only.
# Verifies the basic class machinery: construction, introspection, id_map,
# activeObjects accessor, and spatIDs/featIDs dispatch on global IDs.

.mk_minimal <- function(ncell, nfeat) {
    m <- matrix(0, nrow = nfeat, ncol = ncell)
    rownames(m) <- paste0("f", seq_len(nfeat))
    colnames(m) <- paste0("c", seq_len(ncell))
    createGiottoObject(expression = m, verbose = FALSE)
}

test_that("empty giottoMulti constructs and shows", {
    mg <- new("giottoMulti")
    expect_s4_class(mg, "giottoMulti")
    expect_true(is(mg, "gAny"))
    expect_length(mg, 0L)
    expect_null(idMap(mg, "cells"))
    expect_output(show(mg), "giottoMulti")
})

test_that("populated giottoMulti exposes children", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    expect_s4_class(mg, "giottoMulti")
    expect_true(is(mg, "gAny"))
    expect_identical(names(mg), c("a", "b"))
    expect_length(mg, 2L)
    expect_s4_class(mg[["a"]], "giotto")
    expect_identical(mg[[1]], mg[["a"]])
})

test_that("id_map namespaces cells globally and leaves feats passthrough", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    cells <- idMap(mg, "cells")
    expect_identical(nrow(cells), 8L)
    expect_identical(unique(cells$object), c("a", "b"))
    expect_true(all(grepl("^[ab]::c[0-9]+$", cells$global_id)))

    feats <- idMap(mg, "feats")
    # 4 features per object x 2 objects = 8 rows in long form
    expect_identical(nrow(feats), 8L)
    # but feature names are shared (passthrough), so global = local
    expect_identical(feats$global_id, feats$local_id)
})

test_that("spatIDs returns global by default, local on request, filterable", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    expect_identical(
        spatIDs(mg),
        c("a::c1", "a::c2", "a::c3", "a::c4", "a::c5",
          "b::c1", "b::c2", "b::c3")
    )
    expect_identical(spatIDs(mg, local = TRUE),
        c("c1", "c2", "c3", "c4", "c5", "c1", "c2", "c3"))
    expect_identical(spatIDs(mg, object = "a"),
        c("a::c1", "a::c2", "a::c3", "a::c4", "a::c5"))
    expect_identical(spatIDs(mg, object = "a", local = TRUE),
        c("c1", "c2", "c3", "c4", "c5"))
})

test_that("featIDs returns uniques by default", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    expect_identical(featIDs(mg), c("f1", "f2", "f3", "f4"))
    expect_identical(length(featIDs(mg, uniques = FALSE)), 8L)
})

test_that("activeObjects get/set works and validates", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    # default: all
    expect_identical(activeObjects(mg), c("a", "b"))

    activeObjects(mg) <- "a"
    expect_identical(activeObjects(mg), "a")

    activeObjects(mg) <- NULL
    expect_identical(activeObjects(mg), c("a", "b"))

    expect_error(activeObjects(mg) <- "nope", "unknown object")
})

test_that("[[<- replaces a child", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    g3 <- .mk_minimal(2, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    mg[["b"]] <- g3
    expect_identical(mg[["b"]], g3)
    # NOTE: id_map is not eagerly refreshed; documented behavior
})

test_that("gAny inheritance does not change giotto dispatch", {
    g <- .mk_minimal(5, 4)
    expect_true(is(g, "giotto"))
    expect_true(is(g, "gAny"))
    # spatIDs("giotto", ...) still wins
    expect_identical(spatIDs(g), c("c1", "c2", "c3", "c4", "c5"))
})


# Shared-domain accessors promoted to S4 (gAny) ####

test_that("getExpression works identically on giotto via gAny method", {
    g <- .mk_minimal(5, 4)
    e <- getExpression(g)
    expect_s4_class(e, "exprObj")
    expect_identical(dim(e[]), c(4L, 5L))

    em <- getExpression(g, output = "matrix")
    expect_true(inherits(em, c("matrix", "Matrix")))
})

test_that("getExpression on giottoMulti reads from parent's shared slot", {
    g1 <- .mk_minimal(5, 4)
    mg <- createGiottoMulti(list(a = g1))

    # parent slot is empty; getExpression should error usefully
    expect_error(getExpression(mg))

    # populate parent's shared slot, exercising the gAny method's slot read
    mg@expression <- g1@expression
    e <- getExpression(mg)
    expect_s4_class(e, "exprObj")
    expect_identical(dim(e[]), c(4L, 5L))
})

test_that("getCellMetadata works identically on giotto via gAny method", {
    g <- .mk_minimal(5, 4)
    cm <- getCellMetadata(g)
    expect_s4_class(cm, "cellMetaObj")
    expect_identical(nrow(cm[]), 5L)
})

test_that("getCellMetadata on giottoMulti reads from parent's shared slot", {
    g1 <- .mk_minimal(5, 4)
    mg <- createGiottoMulti(list(a = g1))
    # populate both shared slots — set_default_spat_unit reads @expression to
    # resolve defaults
    mg@expression <- g1@expression
    mg@cell_metadata <- g1@cell_metadata
    cm <- getCellMetadata(mg)
    expect_s4_class(cm, "cellMetaObj")
    expect_identical(nrow(cm[]), 5L)
})


# Spatial-domain accessors — per-child dispatch on giottoMulti ####

test_that("getSpatialLocations on giottoMulti returns named per-child list", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    out <- getSpatialLocations(mg)
    expect_type(out, "list")
    expect_identical(names(out), c("a", "b"))
    expect_s4_class(out$a, "spatLocsObj")
    expect_s4_class(out$b, "spatLocsObj")
    expect_identical(nrow(out$a[]), 5L)
    expect_identical(nrow(out$b[]), 3L)
})

test_that("getSpatialLocations honors object= to subset children", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))

    out <- getSpatialLocations(mg, object = "b")
    expect_identical(names(out), "b")
    expect_identical(nrow(out$b[]), 3L)
})

test_that("setSpatialLocations on giottoMulti requires object= and routes", {
    g1 <- .mk_minimal(5, 4)
    g2 <- .mk_minimal(3, 4)
    mg <- createGiottoMulti(list(a = g1, b = g2))
    sl_b <- getSpatialLocations(g2)

    # missing object should error
    expect_error(setSpatialLocations(mg, x = sl_b),
        "must name the child")

    # length > 1 should error
    expect_error(
        setSpatialLocations(mg, x = sl_b, object = c("a", "b")),
        "length 1"
    )

    # round-trip: writing back the child's own spatlocs returns a giottoMulti
    sl_a <- getSpatialLocations(g1)
    mg2 <- setSpatialLocations(mg, x = sl_a, object = "a", verbose = FALSE)
    expect_s4_class(mg2, "giottoMulti")
    out_a <- getSpatialLocations(mg2, object = "a")$a
    expect_identical(nrow(out_a[]), 5L)
})
