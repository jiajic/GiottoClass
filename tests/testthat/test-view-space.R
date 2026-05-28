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

# fixture — two-sample giottoMulti with each child's cells under distinct
# local IDs (s2_* in the second child) so global IDs are unambiguous.
.fixture_gmulti <- function() {
    g1 <- .fixture_giotto()
    g2 <- .fixture_giotto()
    cm2 <- pDataDT(g2)
    new_ids <- paste0("s2_", cm2$cell_ID)
    g2@cell_metadata$cell$rna@metaDT$cell_ID <- new_ids
    sl <- g2@spatial_locs$cell$raw
    sl@coordinates$cell_ID <- new_ids
    g2@spatial_locs$cell$raw <- sl
    sv <- g2@spatial_info$cell@spatVector
    sv$poly_ID <- new_ids
    g2@spatial_info$cell@spatVector <- sv
    g2@spatial_info$cell@unique_ID_cache <- new_ids
    e <- g2@expression$cell$rna$raw
    colnames(e@exprMat) <- new_ids
    g2@expression$cell$rna$raw <- e
    g2@cell_ID$cell <- new_ids
    createGiottoMulti(list(a = g1, b = g2))
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


# --- spatValues with view = ----------------------------------------------

test_that("spatValues view = name narrows returned rows", {
    g <- .fixture_giotto()
    giottoView(g, "tumor") <- giottoView() |> subset(leiden_clus == "1")
    sv <- spatValues(g, feats = "leiden_clus", view = "tumor")
    expect_equal(nrow(sv), sum(pDataDT(g)$leiden_clus == "1"))
    expect_true(all(sv$leiden_clus == "1"))
})

test_that("spatValues view = ad-hoc giottoView object works", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus %in% c("2", "3"))
    sv <- spatValues(g, feats = "leiden_clus", view = v)
    expect_equal(nrow(sv), sum(pDataDT(g)$leiden_clus %in% c("2", "3")))
    expect_true(all(sv$leiden_clus %in% c("2", "3")))
})

test_that("spatValues view = NULL is identity (matches raw)", {
    g <- .fixture_giotto()
    raw <- spatValues(g, feats = "leiden_clus")
    same <- spatValues(g, feats = "leiden_clus", view = NULL)
    expect_identical(raw, same)
})

test_that("spatValues empty view recipe returns same as raw", {
    g <- .fixture_giotto()
    raw <- spatValues(g, feats = "leiden_clus")
    via_empty <- spatValues(g, feats = "leiden_clus", view = giottoView())
    expect_equal(nrow(via_empty), nrow(raw))
})

test_that("spatValues view = ... matches getCellMetadata view = ... narrowing", {
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    sv <- spatValues(g, feats = "leiden_clus", view = "x")
    cm <- getCellMetadata(g, view = "x", output = "data.table")
    # both paths produce the same cell_ID set
    expect_setequal(sv$cell_ID, cm$cell_ID)
})

test_that("spatValues view = ... is consistent with materialize -> spatValues raw", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus == "1")
    direct <- spatValues(g, feats = "leiden_clus", view = v)
    g_m <- materialize(g, v)
    via_materialize <- spatValues(g_m, feats = "leiden_clus")
    # direct narrowing should match the materialized-then-raw path
    expect_setequal(direct$cell_ID, via_materialize$cell_ID)
})

test_that("spatValues view that filters to zero cells returns empty data.table", {
    g <- .fixture_giotto()
    v <- giottoView() |> subset(leiden_clus == "_nonexistent_cluster_")
    sv <- spatValues(g, feats = "leiden_clus", view = v)
    expect_equal(nrow(sv), 0L)
    expect_true("cell_ID" %in% colnames(sv))
})

test_that("spatValues view = composed predicate AND-narrows correctly", {
    g <- .fixture_giotto()
    # two subset steps chained — both should apply (intersection semantics)
    v <- giottoView() |>
        subset(leiden_clus %in% c("1", "2")) |>
        subset(total_expr > median(total_expr))
    sv <- spatValues(g, feats = "leiden_clus", view = v)
    cm <- pDataDT(g)
    n_expected <- sum(
        cm$leiden_clus %in% c("1", "2") &
            cm$total_expr > median(cm$total_expr))
    expect_equal(nrow(sv), n_expected)
})

test_that("spatValues view = on giottoMulti narrows joint output", {
    mg <- .fixture_gmulti()
    sv_raw <- spatValues(mg, feats = "leiden_clus")
    v <- giottoView() |> subset(leiden_clus == "1")
    sv_v <- spatValues(mg, feats = "leiden_clus", view = v)
    expect_lt(nrow(sv_v), nrow(sv_raw))
    expect_true(all(sv_v$leiden_clus == "1"))
    # global cell_ID format still present
    expect_true(any(grepl("^a::", sv_v$cell_ID)) ||
                any(grepl("^b::", sv_v$cell_ID)))
})

test_that("spatValues view re-entry guard prevents recursion", {
    # Set the option as if we're mid-resolution; an outer spatValues
    # call with view should drop the view arg rather than recurse.
    g <- .fixture_giotto()
    giottoView(g, "x") <- giottoView() |> subset(leiden_clus == "1")
    options(giotto.spatValues_view_active = TRUE)
    on.exit(options(giotto.spatValues_view_active = FALSE), add = TRUE)
    sv_v <- spatValues(g, feats = "leiden_clus", view = "x")
    sv_raw <- spatValues(g, feats = "leiden_clus")
    # with guard active, view= is dropped; result matches raw
    expect_equal(nrow(sv_v), nrow(sv_raw))
})

test_that("spatValues space = NULL is no-op (currently accepted but not value-transforming)", {
    g <- .fixture_giotto()
    giottoSpace(g, "tilted") <- giottoSpace() |> spin(30)
    sv_native <- spatValues(g, feats = "leiden_clus")
    sv_space <- spatValues(g, feats = "leiden_clus", space = "tilted")
    # space does NOT transform value columns (leiden_clus is a label);
    # rows and values are unchanged
    expect_equal(nrow(sv_native), nrow(sv_space))
    expect_setequal(sv_native$leiden_clus, sv_space$leiden_clus)
})


# --- gmulti dispatch ------------------------------------------------------

test_that("giottoView accessors work on giottoMulti via gAny", {
    mg <- .fixture_gmulti()
    v <- giottoView() |> subset(leiden_clus == "1")
    giottoView(mg, "tumor") <- v
    expect_identical(giottoViews(mg), "tumor")
    out <- giottoView(mg, "tumor")
    expect_s4_class(out, "giottoView")
})

test_that("giottoSpace accessors work on giottoMulti via gAny", {
    mg <- .fixture_gmulti()
    s <- (giottoSpace("a") |> spin(30)) + (giottoSpace("b") |> spin(45))
    giottoSpace(mg, "atlas") <- s
    expect_identical(giottoSpaces(mg), "atlas")
    out <- giottoSpace(mg, "atlas")
    expect_named(out@samples, c("a", "b"))
})

test_that(".scope_space_to_sample picks the right key for a child", {
    s <- (giottoSpace("a") |> spin(30)) + (giottoSpace("b") |> spin(45))
    sa <- GiottoClass:::.scope_space_to_sample(s, "a")
    expect_named(sa@samples, GiottoClass:::.space_default_sample)
    expect_length(sa@samples[[1L]], 1L)
    expect_equal(sa@samples[[1L]][[1L]]@op, "spin")
    expect_equal(sa@samples[[1L]][[1L]]@args$angle, 30)

    sb <- GiottoClass:::.scope_space_to_sample(s, "b")
    expect_equal(sb@samples[[1L]][[1L]]@args$angle, 45)
})

test_that(".scope_space_to_sample falls back to :default: key", {
    s <- giottoSpace() |> spin(15)
    out <- GiottoClass:::.scope_space_to_sample(s, "any_sample_name")
    expect_equal(out@samples[[1L]][[1L]]@args$angle, 15)
})

test_that(".scope_space_to_sample returns NULL when no matching key", {
    s <- giottoSpace("only_x") |> spin(15)
    out <- GiottoClass:::.scope_space_to_sample(s, "missing")
    expect_null(out)
})

test_that("materialize on giottoMulti narrows children via selectSamples", {
    mg <- .fixture_gmulti()
    v <- giottoView() |> selectSamples("a")
    out <- materialize(mg, v)
    expect_identical(names(out@objects), "a")
})

test_that("materialize on giottoMulti applies view per-child", {
    mg <- .fixture_gmulti()
    v <- giottoView() |> subset(leiden_clus == "1")
    out <- materialize(mg, v)
    expect_named(out@objects, c("a", "b"))
    # each child has been narrowed
    n_a <- nrow(pDataDT(out@objects$a))
    n_b <- nrow(pDataDT(out@objects$b))
    expect_lt(n_a, length(spatIDs(mg@objects$a)))
    expect_lt(n_b, length(spatIDs(mg@objects$b)))
})

test_that("materialize on giottoMulti scopes space per-child", {
    mg <- .fixture_gmulti()
    s <- (giottoSpace("a") |> spin(30)) + (giottoSpace("b") |> spin(45))
    giottoSpace(mg, "atlas") <- s
    out <- materialize(mg, giottoView(), space = "atlas")

    sl_a_native <- getSpatialLocations(mg@objects$a, output = "data.table")
    sl_b_native <- getSpatialLocations(mg@objects$b, output = "data.table")
    sl_a_post <- getSpatialLocations(out@objects$a, output = "data.table")
    sl_b_post <- getSpatialLocations(out@objects$b, output = "data.table")

    # both children transformed — but by different angles
    expect_false(isTRUE(all.equal(sl_a_post$sdimx, sl_a_native$sdimx)))
    expect_false(isTRUE(all.equal(sl_b_post$sdimx, sl_b_native$sdimx)))
    # different angles → ratios of transformed-to-native should differ
    # (we don't compute the exact rotation expectation here, just that
    # the two children's transforms are not identical)
    expect_false(isTRUE(all.equal(
        sl_a_post$sdimx - sl_a_native$sdimx,
        sl_b_post$sdimx - sl_b_native$sdimx)))
})

test_that("spatValues on giottoMulti finds features in joint cell_metadata", {
    mg <- .fixture_gmulti()
    sv <- spatValues(mg, feats = "leiden_clus")
    expect_equal(nrow(sv), sum(lengths(lapply(mg@objects, spatIDs))))
    expect_true(all(c("cell_ID", "leiden_clus") %in% colnames(sv)))
    # global cell_IDs use the sample::local_id format
    expect_true(any(grepl("^a::", sv$cell_ID)))
    expect_true(any(grepl("^b::", sv$cell_ID)))
})

test_that("spatValues on giottoMulti finds features in joint expression", {
    mg <- .fixture_gmulti()
    gene <- rownames(getExpression(mg, output = "matrix"))[1L]
    sv <- spatValues(mg, feats = gene)
    expect_equal(nrow(sv), sum(lengths(lapply(mg@objects, spatIDs))))
    expect_true(gene %in% colnames(sv))
})

test_that("materialize on giottoMulti narrows joint shared slots", {
    mg <- .fixture_gmulti()
    v <- giottoView() |> subset(leiden_clus == "1")
    n_total <- nrow(pDataDT(mg))
    n_target <- sum(pDataDT(mg)$leiden_clus == "1")

    out <- materialize(mg, v)

    # joint cell_metadata narrowed
    expect_equal(nrow(pDataDT(out)), n_target)
    expect_lt(nrow(pDataDT(out)), n_total)
    # joint expression narrowed in column count
    expect_equal(ncol(getExpression(out, output = "matrix")), n_target)
})

test_that("materialize via slotted view name dispatches on multi", {
    mg <- .fixture_gmulti()
    giottoView(mg, "x") <- giottoView() |> selectSamples("a")
    out <- materialize(mg, "x")
    expect_identical(names(out@objects), "a")
})
