# hnswKNN: approximate kNN search returning a dbscan-compatible shape

.nn_mat <- function(n = 500L, d = 10L, seed = 1L) {
    set.seed(seed)
    matrix(stats::rnorm(n * d), nrow = n, ncol = d)
}

test_that("hnswKNN returns the same shape as dbscan::kNN", {
    m <- .nn_mat()
    k <- 15L

    hn <- hnswKNN(m, k = k)
    ex <- dbscan::kNN(m, k = k, sort = TRUE)

    expect_identical(class(hn), class(ex))
    expect_true(all(c("id", "dist", "k", "sort", "metric") %in% names(hn)))
    expect_identical(dim(hn$id), dim(ex$id))
    expect_identical(dim(hn$dist), dim(ex$dist))
    expect_identical(hn$k, k)
    expect_type(hn$id, "integer")
})

test_that("hnswKNN excludes self and returns sorted distances", {
    m <- .nn_mat()
    hn <- hnswKNN(m, k = 10L)

    expect_false(any(hn$id == seq_len(nrow(m))))
    expect_true(all(apply(hn$dist, 1L, function(r) !is.unsorted(r))))
})

test_that("hnswKNN drops the correct entry when self is not column 1", {
    # Duplicated coordinates put a zero-distance twin alongside the self-hit,
    # so the self can land at any column. A column-major mask extract
    # misaligns rows here while looking correct when self is always first.
    m <- .nn_mat(n = 100L, d = 5L)
    md <- rbind(m, m) # every row has an exact duplicate

    hn <- hnswKNN(md, k = 5L)

    expect_false(any(hn$id == seq_len(nrow(md))))
    expect_identical(dim(hn$id), c(nrow(md), 5L))
    expect_false(anyNA(hn$id))
})

test_that("hnswKNN recall is high against exact search", {
    m <- .nn_mat(n = 1000L, d = 15L)
    k <- 20L

    hn <- hnswKNN(m, k = k)
    ex <- dbscan::kNN(m, k = k, sort = TRUE)
    recall <- mean(vapply(seq_len(nrow(m)), function(i) {
        length(intersect(hn$id[i, ], ex$id[i, ])) / k
    }, numeric(1L)))

    expect_gt(recall, 0.95)
})

test_that("hnswKNN is reproducible at the default n_threads", {
    # The default is serial precisely because hnswlib's index build is only
    # reproducible single-threaded. If this ever fails, the default changed.
    m <- .nn_mat(n = 800L, d = 12L)
    expect_identical(hnswKNN(m, k = 10L), hnswKNN(m, k = 10L))
})

test_that("dbscan::sNN consumes an hnswKNN result", {
    m <- .nn_mat()
    hn <- hnswKNN(m, k = 15L)

    snn <- dbscan::sNN(x = hn, k = 15L, kt = NULL)
    expect_true(all(c("shared", "id", "dist") %in% names(snn)))
    expect_identical(dim(snn$shared), dim(hn$id))
})

test_that("hnswKNN rejects k >= nrow(x)", {
    m <- .nn_mat(n = 20L)
    expect_error(hnswKNN(m, k = 20L), "must be less than")
})

test_that("engine resolution follows the declared space", {
    # kNN/sNN params default to "auto"; the giotto method fills it from `space`
    expect_identical(kNNNetworkParam()@engine, "auto")
    expect_identical(sNNNetworkParam()@engine, "auto")
    expect_identical(kNNNetworkParam(engine = "dbscan")@engine, "dbscan")

    # spatial wrapper declares dbscan for its 2-3 dimensional coordinates
    g <- GiottoData::loadGiottoMini("visium", verbose = FALSE)
    sn <- createSpatialKNNnetwork(g, k = 4L, return_gobject = FALSE,
        output = "data.table")
    expect_s3_class(sn, "data.table")
})
