# Tests for giottoView + giottoSpace classes, the resolver engine, and
# the JIT view/space integration in getters.
#
# silence deprecated internal functions
rlang::local_options(lifecycle_verbosity = "quiet")
options("giotto.use_conda" = FALSE)

# fixture — visium mini with leiden clusters in metadata
.fixture_giotto <- function() {
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    updateGiottoObject(g)
}


# --- giottoView class ------------------------------------------------------

test_that("giottoView() constructs empty view", {
    v <- giottoView()
    expect_s4_class(v, "giottoView")
    expect_length(v@steps, 0L)
    expect_true(is.na(v@name))
    expect_true(is.na(v@space))
})

test_that("giottoView(space = ...) records space reference", {
    v <- giottoView(space = "atlas")
    expect_identical(v@space, "atlas")
})

test_that("subset() records a viewFilter step", {
    v <- giottoView() |> subset(cluster == "A")
    expect_length(v@steps, 1L)
    expect_s4_class(v@steps[[1L]], "viewFilter")
    expect_identical(deparse(v@steps[[1L]]@predicate), 'cluster == "A"')
})

test_that("subset() forwards scope_args", {
    v <- giottoView() |> subset(x > 0, spat_unit = "cell",
        feat_type = "rna", negate = TRUE)
    sa <- v@steps[[1L]]@scope_args
    expect_identical(sa$spat_unit, "cell")
    expect_identical(sa$feat_type, "rna")
    expect_true(sa$negate)
})

test_that("crop() records a viewCrop step", {
    v <- giottoView() |> crop(c(0, 100, 0, 100))
    expect_length(v@steps, 1L)
    expect_s4_class(v@steps[[1L]], "viewCrop")
    expect_equal(v@steps[[1L]]@extent, c(0, 100, 0, 100))
})

test_that("selectSamples() records a viewSampleSelect step", {
    v <- giottoView() |> selectSamples("a", "b")
    expect_length(v@steps, 1L)
    expect_s4_class(v@steps[[1L]], "viewSampleSelect")
    expect_identical(v@steps[[1L]]@samples, c("a", "b"))
})

test_that("steps compose in order under pipe", {
    v <- giottoView() |>
        subset(x > 0) |>
        crop(c(0, 100, 0, 100)) |>
        selectSamples("a")
    expect_length(v@steps, 3L)
    expect_s4_class(v@steps[[1L]], "viewFilter")
    expect_s4_class(v@steps[[2L]], "viewCrop")
    expect_s4_class(v@steps[[3L]], "viewSampleSelect")
})


# --- giottoSpace class -----------------------------------------------------

test_that("giottoSpace() constructs empty sample-anonymous space", {
    s <- giottoSpace()
    expect_s4_class(s, "giottoSpace")
    expect_named(s@samples, ":default:")
    expect_length(s@samples[[1L]], 0L)
})

test_that("giottoSpace(sample) constructs sample-bound space", {
    s <- giottoSpace("sample_a")
    expect_named(s@samples, "sample_a")
})

test_that("transform generics record on giottoSpace", {
    M <- diag(c(1, 1, 1))
    s <- giottoSpace() |> spin(30) |> affine(M) |> spatShift(dx = 10)
    steps <- s@samples[[1L]]
    expect_length(steps, 3L)
    expect_identical(vapply(steps, function(x) x@op, character(1L)),
        c("spin", "affine", "spatShift"))
})

test_that("spin/affine record (0,0) anchor by default", {
    s <- giottoSpace() |> spin(45)
    expect_equal(s@samples[[1L]][[1L]]@args$x0, 0)
    expect_equal(s@samples[[1L]][[1L]]@args$y0, 0)
})

test_that("user-supplied anchor overrides default", {
    s <- giottoSpace() |> spin(45, x0 = 100, y0 = 200)
    expect_equal(s@samples[[1L]][[1L]]@args$x0, 100)
    expect_equal(s@samples[[1L]][[1L]]@args$y0, 200)
})


# --- + composition on giottoSpace -----------------------------------------

test_that("+ on same-sample concatenates step lists", {
    s <- (giottoSpace("a") |> spin(30)) +
         (giottoSpace("a") |> spatShift(dx = 10))
    expect_length(s@samples, 1L)
    expect_length(s@samples[["a"]], 2L)
})

test_that("+ on different samples merges keyed", {
    s <- (giottoSpace("a") |> spin(30)) +
         (giottoSpace("b") |> spin(45))
    expect_named(s@samples, c("a", "b"))
    expect_length(s@samples[["a"]], 1L)
    expect_length(s@samples[["b"]], 1L)
})

test_that("+ giottoView+giottoView errors (not yet implemented)", {
    expect_error(giottoView() + giottoView(),
        "not yet implemented")
})


# --- Accessors -------------------------------------------------------------

test_that("giottoView<- slots in and giottoView() retrieves by name", {
    g <- giotto()
    v <- giottoView() |> subset(cluster == "A")
    giottoView(g, "tumor") <- v
    expect_identical(giottoViews(g), "tumor")
    out <- giottoView(g, "tumor")
    expect_s4_class(out, "giottoView")
    expect_identical(out@name, "tumor")
})

test_that("giottoView(g, name) <- NULL removes", {
    g <- giotto()
    giottoView(g, "a") <- giottoView()
    giottoView(g, "b") <- giottoView()
    giottoView(g, "a") <- NULL
    expect_identical(giottoViews(g), "b")
})

test_that("giottoSpace accessor and lookup", {
    g <- giotto()
    giottoSpace(g, "atlas") <- giottoSpace() |> spin(30)
    expect_identical(giottoSpaces(g), "atlas")
    out <- giottoSpace(g, "atlas")
    expect_s4_class(out, "giottoSpace")
    expect_identical(out@name, "atlas")
})

test_that("missing slotted name errors clearly", {
    g <- giotto()
    expect_error(giottoView(g, "missing"), "no slotted giottoView")
    expect_error(giottoSpace(g, "missing"), "no slotted giottoSpace")
})


# --- Migration on updateGiottoObject --------------------------------------

test_that("updateGiottoObject() adds @view and @spaces for pre-0.7.0", {
    g <- giotto()
    g@versions$gclass <- "0.6.0"
    g <- updateGiottoObject(g)
    expect_true(methods::.hasSlot(g, "view"))
    expect_true(methods::.hasSlot(g, "spaces"))
    expect_null(g@view)
    expect_null(g@spaces)
})

test_that("save/load round-trip preserves @view and @spaces", {
    g <- giotto()
    giottoView(g, "demo") <- giottoView() |> subset(x > 0)
    giottoSpace(g, "tilted") <- giottoSpace() |> spin(45)

    td <- tempfile("gv-")
    on.exit(unlink(file.path(dirname(td), basename(td)), recursive = TRUE),
        add = TRUE)

    saveGiotto(g, foldername = basename(td), dir = dirname(td),
        verbose = FALSE, overwrite = TRUE)
    g2 <- loadGiotto(file.path(dirname(td), basename(td)), verbose = FALSE)

    expect_identical(giottoViews(g2), "demo")
    expect_identical(giottoSpaces(g2), "tilted")
})


# --- Coordinator: dataTableCoordinator ---------------------------------------

test_that("dataTableCoordinator() constructs", {
    p <- dataTableCoordinator()
    expect_s4_class(p, "dataTableCoordinator")
    expect_s4_class(p, "viewCoordinator")
})

test_that(".default_view_coordinator returns dataTableCoordinator for in-memory", {
    g <- giotto()
    p <- GiottoClass:::.default_view_coordinator(g)
    expect_s4_class(p, "dataTableCoordinator")
})

test_that("prepareIds() for dataTableCoordinator is identity", {
    ids <- c("a", "b", "c")
    expect_identical(prepareIds(dataTableCoordinator(), ids), ids)
})


# --- materialize() end-to-end on visium mini ------------------------------

test_that("materialize() with empty view returns equivalent gobject", {
    g <- .fixture_giotto()
    g2 <- materialize(g, giottoView())
    expect_equal(nrow(pDataDT(g2)), nrow(pDataDT(g)))
    expect_equal(
        nrow(getSpatialLocations(g2, output = "data.table")),
        nrow(getSpatialLocations(g, output = "data.table"))
    )
})

test_that("materialize() narrows tabular slots by subset predicate", {
    g <- .fixture_giotto()
    n_total <- length(spatIDs(g))
    n_target <- sum(pDataDT(g)$leiden_clus == "1")

    v <- giottoView() |> subset(leiden_clus == "1")
    g2 <- materialize(g, v)

    expect_lt(nrow(pDataDT(g2)), n_total)
    expect_equal(nrow(pDataDT(g2)), n_target)
    expect_equal(ncol(getExpression(g2, output = "matrix")), n_target)
})

test_that("materialize() narrows spatial slots via cell_ID cascade", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus == "1")
    g2 <- materialize(g, v)

    n_filter <- nrow(pDataDT(g2))
    expect_equal(
        nrow(getSpatialLocations(g2, output = "data.table")), n_filter)
    expect_equal(
        length(spatIDs(getPolygonInfo(g2, return_giottoPolygon = TRUE))),
        n_filter)
})

test_that("materialize() with %in% and env-resident value works (NSE)", {
    g <- .fixture_giotto()
    targets <- c("1", "2")
    v <- giottoView() |> subset(leiden_clus %in% targets)
    g2 <- materialize(g, v)
    expected <- sum(pDataDT(g)$leiden_clus %in% targets)
    expect_equal(nrow(pDataDT(g2)), expected)
})

test_that("materialize() with expression-column predicate routes via spatValues", {
    g <- .fixture_giotto()
    # pick a gene known to be in the panel by literal name to avoid NSE
    gene <- "Gfap"
    skip_if_not(gene %in% rownames(getExpression(g, output = "matrix")),
        sprintf("gene %s not in panel", gene))

    expected <- sum(getExpression(g, output = "matrix")[gene, ] > 0)
    v <- giottoView() |> subset(Gfap > 0)
    g2 <- materialize(g, v)
    expect_equal(nrow(pDataDT(g2)), expected)
})

test_that("materialize() with crop narrows via spatLocs extent", {
    g <- .fixture_giotto()
    sl <- getSpatialLocations(g, output = "data.table")
    ext <- c(4000, 5500, -5000, -3500)
    expected <- sum(sl$sdimx >= ext[1L] & sl$sdimx <= ext[2L] &
                    sl$sdimy >= ext[3L] & sl$sdimy <= ext[4L])

    v <- giottoView() |> crop(ext)
    g2 <- materialize(g, v)
    expect_equal(nrow(pDataDT(g2)), expected)
    expect_equal(
        nrow(getSpatialLocations(g2, output = "data.table")), expected)
})

test_that("materialize() with space transforms spatial coords only", {
    g <- .fixture_giotto()
    giottoSpace(g, "tilted") <- giottoSpace() |> spin(30)

    g2 <- materialize(g, giottoView(), space = "tilted")
    sl_native <- getSpatialLocations(g, output = "data.table")
    sl_tilted <- getSpatialLocations(g2, output = "data.table")

    expect_equal(nrow(sl_tilted), nrow(sl_native))
    expect_false(isTRUE(all.equal(sl_tilted$sdimx, sl_native$sdimx)))
    # tabular slots unchanged in row count
    expect_equal(nrow(pDataDT(g2)), nrow(pDataDT(g)))
})

test_that("materialize() with filter + space combines both", {
    g <- .fixture_giotto()
    giottoSpace(g, "tilted") <- giottoSpace() |> spin(30)

    v <- giottoView() |> subset(leiden_clus %in% c("1", "2"))
    g2 <- materialize(g, v, space = "tilted")

    expected <- sum(pDataDT(g)$leiden_clus %in% c("1", "2"))
    sl_tilted <- getSpatialLocations(g2, output = "data.table")
    sl_native <- getSpatialLocations(g, output = "data.table")
    expect_equal(nrow(sl_tilted), expected)
    expect_false(isTRUE(all.equal(sl_tilted$sdimx, sl_native$sdimx[1:expected])))
})


# --- JIT getter integration -----------------------------------------------

test_that("getCellMetadata respects view = name", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    n_base <- nrow(getCellMetadata(g, output = "data.table"))
    n_view <- nrow(getCellMetadata(g, view = "x", output = "data.table"))
    expect_lt(n_view, n_base)
    expect_equal(n_view, sum(pDataDT(g)$leiden_clus == "1"))
})

test_that("getCellMetadata accepts ad-hoc view object", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus == "2")
    n <- nrow(getCellMetadata(g, view = v, output = "data.table"))
    expect_equal(n, sum(pDataDT(g)$leiden_clus == "2"))
})

test_that("getExpression view narrows columns", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    n_base <- ncol(getExpression(g, output = "matrix"))
    n_view <- ncol(getExpression(g, view = "x", output = "matrix"))
    expect_lt(n_view, n_base)
})

test_that("getSpatialLocations view + space combined", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    giottoSpace(g, "tilted") <- giottoSpace() |> spin(30)

    sl_base <- getSpatialLocations(g, output = "data.table")
    sl_v <- getSpatialLocations(g, view = "x", output = "data.table")
    sl_s <- getSpatialLocations(g, space = "tilted", output = "data.table")
    sl_vs <- getSpatialLocations(g, view = "x", space = "tilted",
        output = "data.table")

    expect_equal(nrow(sl_v), sum(pDataDT(g)$leiden_clus == "1"))
    expect_false(isTRUE(all.equal(sl_s$sdimx, sl_base$sdimx)))
    expect_equal(nrow(sl_vs), nrow(sl_v))
    # rotated x for the view+space combo differs from the unrotated view
    expect_false(isTRUE(all.equal(sl_vs$sdimx, sl_v$sdimx)))
})

test_that("getPolygonInfo view narrows; both output forms agree", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    gp_full <- getPolygonInfo(g, view = "x", return_giottoPolygon = TRUE)
    sv_full <- getPolygonInfo(g, view = "x")  # SpatVector default
    expect_equal(length(spatIDs(gp_full)), nrow(sv_full))
})

test_that("getFeatureMetadata view is no-op (feat-keyed)", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    n_base <- nrow(getFeatureMetadata(g, output = "data.table"))
    n_view <- nrow(getFeatureMetadata(g, view = "x", output = "data.table"))
    expect_equal(n_view, n_base)
})

test_that("getter without view/space returns unchanged baseline", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    expect_equal(
        nrow(getCellMetadata(g, output = "data.table")),
        length(spatIDs(g))
    )
})


# --- Resolver cache --------------------------------------------------------

test_that(".cached_surviving_cell_ids memoises within a cache env", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus == "1")
    cache <- GiottoClass:::.new_resolver_cache()

    a <- GiottoClass:::.cached_surviving_cell_ids(g, v, NULL,
        dataTableCoordinator(), cache)
    expect_true(exists("surviving_ids", envir = cache))
    b <- GiottoClass:::.cached_surviving_cell_ids(g, v, NULL,
        dataTableCoordinator(), cache)
    expect_identical(a, b)
})

test_that(".cached_surviving_cell_ids with NULL cache works", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus == "1")
    ids <- GiottoClass:::.cached_surviving_cell_ids(g, v, NULL,
        dataTableCoordinator(), NULL)
    expect_type(ids, "character")
    expect_equal(length(ids), sum(pDataDT(g)$leiden_clus == "1"))
})
