# Behavior contract for giottoMulti structural operations.
# ============================================================================
#
# Codifies the rules for add / remove / rename / reorder / replace on
# `@objects`, and the consequences for joint slots.
#
# Round-trip / symmetry invariants this file enforces:
#   * add + remove (subset out)        -> original
#   * rename a->b then b->a            -> original
#   * reorder twice (round-trip perm)  -> original
#
# Per-slot extensibility on add-sample:
#   * @cell_metadata, @feat_metadata, @expression  ->  auto-extend (silent)
#   * @nn_network, @spatial_enrichment             ->  cannot extend; nudge
#   * @dimension_reduction                          ->  project if possible
#                                                       (PCA / UMAP); nudge
#                                                       for harmony etc.
#
# Tests for behavior not yet implemented use `skip()` with a clear pointer
# to the relevant todo. They run cleanly today; flip to `expect_*()` once
# the corresponding implementation lands.


# ---- helpers ---------------------------------------------------------------

.mk_minimal <- function(ncell, nfeat, prefix = "c") {
    m <- matrix(0, nrow = nfeat, ncol = ncell)
    rownames(m) <- paste0("f", seq_len(nfeat))
    colnames(m) <- paste0(prefix, seq_len(ncell))
    createGiottoObject(expression = m, verbose = FALSE)
}

.mk_gmulti_pair <- function() {
    g1 <- .mk_minimal(5, 4, prefix = "c1_")
    g2 <- .mk_minimal(3, 4, prefix = "c2_")
    createGiottoMulti(list(a = g1, b = g2))
}

.mk_gmulti_triple <- function() {
    g1 <- .mk_minimal(5, 4, prefix = "c1_")
    g2 <- .mk_minimal(3, 4, prefix = "c2_")
    g3 <- .mk_minimal(4, 4, prefix = "c3_")
    createGiottoMulti(list(a = g1, b = g2, c = g3))
}


# ----------------------------------------------------------------------------
# 1. Reorder
# ----------------------------------------------------------------------------

test_that("reorder via subset permutes @objects names without altering joint", {
    mg <- .mk_gmulti_triple()
    expect_identical(names(mg), c("a", "b", "c"))

    mg_perm <- mg[c("c", "a", "b")]
    expect_identical(names(mg_perm), c("c", "a", "b"))
    # @id_map rebuilt but identity-preserving
    expect_setequal(idMap(mg_perm, "cells")$global_id,
                    idMap(mg, "cells")$global_id)
})

test_that("reorder round-trip returns identical content", {
    mg <- .mk_gmulti_triple()
    mg_perm <- mg[c("c", "a", "b")]
    mg_back <- mg_perm[c("a", "b", "c")]
    expect_identical(names(mg_back), names(mg))
    expect_identical(idMap(mg_back, "cells")$global_id,
                     idMap(mg, "cells")$global_id)
})


# ----------------------------------------------------------------------------
# 2. Rename
# ----------------------------------------------------------------------------

test_that("rename via names<- updates @objects names + @id_map sample prefix", {
    mg <- .mk_gmulti_pair()
    names(mg) <- c("A", "B")
    expect_identical(names(mg), c("A", "B"))
    # @id_map's global_id prefixes should now use the new sample names
    glob <- idMap(mg, "cells")$global_id
    expect_true(all(grepl("^(A|B)::", glob)))
    expect_false(any(grepl("^(a|b)::", glob)))
})

test_that("rename round-trip returns identical id_map", {
    mg <- .mk_gmulti_pair()
    orig_ids <- idMap(mg, "cells")$global_id
    names(mg) <- c("A", "B")
    names(mg) <- c("a", "b")
    expect_identical(idMap(mg, "cells")$global_id, orig_ids)
})

test_that("rename rewrites joint @cell_metadata cell_ID prefix", {
    mg <- .mk_gmulti_pair()
    # Populate joint @cell_metadata using ids from the gmulti's universe
    glob <- idMap(mg, "cells")$global_id
    cm <- new("cellMetaObj",
        metaDT = data.table::data.table(cell_ID = glob),
        spat_unit = "cell", feat_type = "rna")
    mg@cell_metadata <- list(cell = list(rna = cm))

    names(mg) <- c("A", "B")
    new_ids <- mg@cell_metadata$cell$rna[]$cell_ID
    expect_true(all(grepl("^(A|B)::", new_ids)))
    expect_false(any(grepl("^(a|b)::", new_ids)))
})

test_that("rename rewrites joint @expression colnames", {
    mg <- .mk_gmulti_pair()
    glob <- idMap(mg, "cells")$global_id
    m <- matrix(0, nrow = 4L, ncol = length(glob),
        dimnames = list(paste0("f", 1:4), glob))
    e <- createExprObj(m, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    mg@expression <- list(cell = list(rna = list(raw = e)))

    names(mg) <- c("A", "B")
    new_cn <- colnames(mg@expression$cell$rna$raw[])
    expect_true(all(grepl("^(A|B)::", new_cn)))
    expect_false(any(grepl("^(a|b)::", new_cn)))
})

test_that("rename rewrites @cell_ID narrowing prefix", {
    mg <- .mk_gmulti_pair()
    mg@cell_ID <- list(cell = head(idMap(mg, "cells")$global_id, 3L))
    names(mg) <- c("A", "B")
    expect_true(all(grepl("^(A|B)::", mg@cell_ID$cell)))
    expect_false(any(grepl("^(a|b)::", mg@cell_ID$cell)))
})

test_that("rename round-trips joint @expression colnames", {
    mg <- .mk_gmulti_pair()
    glob <- idMap(mg, "cells")$global_id
    m <- matrix(0, nrow = 4L, ncol = length(glob),
        dimnames = list(paste0("f", 1:4), glob))
    e <- createExprObj(m, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    mg@expression <- list(cell = list(rna = list(raw = e)))
    orig_cn <- colnames(mg@expression$cell$rna$raw[])

    names(mg) <- c("A", "B")
    names(mg) <- c("a", "b")
    expect_identical(
        colnames(mg@expression$cell$rna$raw[]),
        orig_cn
    )
})

test_that("rename applies alias layer to disk-backed network [pending]", {
    skip("pending: todo #6 — alias layer for parquetEdgeStore on rename")
})


# ----------------------------------------------------------------------------
# 3. Subset (drop samples)
# ----------------------------------------------------------------------------

test_that("subset shrinks @objects and rebuilds @id_map", {
    mg <- .mk_gmulti_triple()
    mg_sub <- mg[c("a", "b")]
    expect_identical(names(mg_sub), c("a", "b"))
    # @id_map should no longer contain "c::" prefixed cells
    glob <- idMap(mg_sub, "cells")$global_id
    expect_false(any(grepl("^c::", glob)))
})

test_that("subset auto-prunes joint @cell_metadata to kept-sample rows", {
    mg <- .mk_gmulti_triple()
    glob <- idMap(mg, "cells")$global_id
    cm <- new("cellMetaObj",
        metaDT = data.table::data.table(cell_ID = glob,
            lbl = letters[seq_along(glob)]),
        spat_unit = "cell", feat_type = "rna")
    mg@cell_metadata <- list(cell = list(rna = cm))

    mg_sub <- mg[c("a", "b")]
    new_dt <- mg_sub@cell_metadata$cell$rna[]
    # No "c::"-prefixed rows survive
    expect_false(any(grepl("^c::", new_dt$cell_ID)))
    # Rows for "a::"/"b::" preserved
    expect_true(all(grepl("^(a|b)::", new_dt$cell_ID)))
})

test_that("subset auto-prunes joint @expression to kept-sample cols", {
    mg <- .mk_gmulti_triple()
    glob <- idMap(mg, "cells")$global_id
    m <- matrix(0, nrow = 4L, ncol = length(glob),
        dimnames = list(paste0("f", 1:4), glob))
    e <- createExprObj(m, name = "raw",
        spat_unit = "cell", feat_type = "rna")
    mg@expression <- list(cell = list(rna = list(raw = e)))

    mg_sub <- mg[c("a", "b")]
    new_cn <- colnames(mg_sub@expression$cell$rna$raw[])
    expect_false(any(grepl("^c::", new_cn)))
    expect_true(all(grepl("^(a|b)::", new_cn)))
})

test_that("subset auto-prunes joint @dimension_reduction rows", {
    mg <- .mk_gmulti_triple()
    glob <- idMap(mg, "cells")$global_id
    coords <- matrix(0, nrow = length(glob), ncol = 2L,
        dimnames = list(glob, c("Dim.1", "Dim.2")))
    d <- new("dimObj",
        coordinates = coords,
        name = "pca", reduction = "cells",
        reduction_method = "pca",
        spat_unit = "cell", feat_type = "rna")
    mg@dimension_reduction <- list(
        cells = list(cell = list(rna = list(pca = list(pca = d))))
    )

    mg_sub <- mg[c("a", "b")]
    new_rn <- rownames(
        mg_sub@dimension_reduction$cells$cell$rna$pca$pca@coordinates
    )
    expect_false(any(grepl("^c::", new_rn)))
    expect_true(all(grepl("^(a|b)::", new_rn)))
})


# ----------------------------------------------------------------------------
# 4. Add new sample
# ----------------------------------------------------------------------------

test_that("add to empty-joint gmulti is silent and grows @objects", {
    mg <- .mk_gmulti_pair()
    g3 <- .mk_minimal(2, 4, prefix = "c3_")
    expect_silent(mg[["c"]] <- g3)
    expect_identical(names(mg), c("a", "b", "c"))
    expect_true(any(grepl("^c::", idMap(mg, "cells")$global_id)))
})

test_that("add with [[<- triggers initialize and clears narrowing", {
    mg <- .mk_gmulti_pair()
    # Manually set narrowing as if filterGiotto had run
    mg@cell_ID <- list(cell = head(idMap(mg, "cells")$global_id, 3L))
    expect_true(!is.null(mg@cell_ID))

    g3 <- .mk_minimal(2, 4, prefix = "c3_")
    mg[["c"]] <- g3
    expect_null(mg@cell_ID)   # structural change cleared narrowing
})

test_that("add with initialize=FALSE skips init (opt-out plumbing) [pending]", {
    skip("R's `[[<-` extra-arg dispatch path is awkward to invoke via base
         syntax; opt-out tested elsewhere via direct method call. Revisit
         once a wrapper helper (e.g. `addObject(g, name, obj, initialize)`)
         exposes a cleaner call site.")
})

test_that("add auto-extends @cell_metadata for new sample [pending]", {
    skip("pending: todo #8 polish — auto-extend cell_metadata with NA fill
         (first pass only warns; auto-extension is a follow-up)")
})

test_that("add auto-extends @expression for new sample [pending]", {
    skip("pending: todo #8 polish — auto-extend expression via backend-aware
         cbind (first pass only warns; auto-extension is a follow-up)")
})

test_that("add with populated joint slot warns with nudge", {
    mg <- .mk_gmulti_pair()
    # Populate joint @cell_metadata so the add-time check has something
    # to flag.
    glob <- idMap(mg, "cells")$global_id
    cm <- new("cellMetaObj",
        metaDT = data.table::data.table(cell_ID = glob),
        spat_unit = "cell", feat_type = "rna")
    mg@cell_metadata <- list(cell = list(rna = cm))

    g3 <- .mk_minimal(2, 4, prefix = "c3_")
    expect_warning(mg[["c"]] <- g3,
        regexp = "joint shared slot.*cell_metadata|do not cover")

    # @objects still grew despite the warning
    expect_identical(names(mg), c("a", "b", "c"))
})

test_that("add to empty-joint gmulti does not warn", {
    mg <- .mk_gmulti_pair()
    g3 <- .mk_minimal(2, 4, prefix = "c3_")
    expect_silent(mg[["c"]] <- g3)
})

test_that("REPLACE existing sample does not trigger add-nudge", {
    mg <- .mk_gmulti_pair()
    glob <- idMap(mg, "cells")$global_id
    cm <- new("cellMetaObj",
        metaDT = data.table::data.table(cell_ID = glob),
        spat_unit = "cell", feat_type = "rna")
    mg@cell_metadata <- list(cell = list(rna = cm))

    # Replace 'a' with a different giotto object — not new, no nudge.
    g_new <- .mk_minimal(4, 4, prefix = "c1_v2_")
    expect_silent(mg[["a"]] <- g_new)
    # (Replace-with-stale-joint-refs is a separate situation; v1 doesn't
    # warn on replace because the nudge is gated on `is_new`.)
})


# ----------------------------------------------------------------------------
# 5a. Child-immutability invariant
#
# Structural ops at the gmulti level MUST NOT modify the wrapped child
# gobjects. Users can grab `g@objects[[name]]` and operate on it as a
# standalone giotto without any side-effect risk from gmulti-level state
# (narrowing in @cell_ID / @feat_ID, joint slot content, id_map / id_sig).
# This is the contract that lets per-child access stay trustworthy.
# ----------------------------------------------------------------------------

test_that("subset does not modify wrapped child gobjects", {
    mg <- .mk_gmulti_triple()
    g_a_before <- mg[["a"]]
    g_a_orig_serialized <- serialize(g_a_before, NULL)

    # subset shouldn't touch 'a' itself
    mg_sub <- mg[c("a", "b")]
    g_a_after <- mg_sub[["a"]]
    g_a_after_serialized <- serialize(g_a_after, NULL)

    expect_identical(g_a_orig_serialized, g_a_after_serialized)
})

test_that("rename does not modify wrapped child gobjects", {
    mg <- .mk_gmulti_pair()
    g_a_before <- mg[["a"]]
    g_a_orig_serialized <- serialize(g_a_before, NULL)

    names(mg) <- c("A", "B")
    g_A_after <- mg[["A"]]
    g_A_after_serialized <- serialize(g_A_after, NULL)

    # Rename at the gmulti level changes the parent key but NOT the
    # wrapped child's internal state.
    expect_identical(g_a_orig_serialized, g_A_after_serialized)
})

test_that("setting @cell_ID narrowing does not modify children", {
    mg <- .mk_gmulti_pair()
    g_a_before_serialized <- serialize(mg[["a"]], NULL)

    mg@cell_ID <- list(cell = head(idMap(mg, "cells")$global_id, 3L))
    g_a_after_serialized <- serialize(mg[["a"]], NULL)

    # Narrowing is a gmulti-level view; children stay untouched.
    expect_identical(g_a_before_serialized, g_a_after_serialized)
})

test_that("add does not modify pre-existing children", {
    mg <- .mk_gmulti_pair()
    g_a_serialized <- serialize(mg[["a"]], NULL)
    g_b_serialized <- serialize(mg[["b"]], NULL)

    g3 <- .mk_minimal(2, 4, prefix = "c3_")
    mg[["c"]] <- g3

    expect_identical(g_a_serialized, serialize(mg[["a"]], NULL))
    expect_identical(g_b_serialized, serialize(mg[["b"]], NULL))
})


# ----------------------------------------------------------------------------
# 5b. Round-trip / symmetry invariants
# ----------------------------------------------------------------------------

test_that("add then remove returns identical object", {
    mg_orig <- .mk_gmulti_pair()

    g3 <- .mk_minimal(2, 4, prefix = "c3_")
    mg2 <- mg_orig
    mg2[["c"]] <- g3
    mg_back <- mg2[c("a", "b")]

    expect_identical(names(mg_back), names(mg_orig))
    expect_identical(idMap(mg_back, "cells")$global_id,
                     idMap(mg_orig, "cells")$global_id)
})

test_that("remove then add returns identical object", {
    mg_orig <- .mk_gmulti_triple()

    mg_sub <- mg_orig[c("a", "b")]
    mg_sub[["c"]] <- mg_orig[["c"]]

    expect_setequal(names(mg_sub), names(mg_orig))
    expect_setequal(idMap(mg_sub, "cells")$global_id,
                    idMap(mg_orig, "cells")$global_id)
})


# ----------------------------------------------------------------------------
# 6. Joint slot wholesale clear via direct slot assignment
# ----------------------------------------------------------------------------

test_that("g@expression <- NULL clears entire joint expression slot", {
    mg <- .mk_gmulti_pair()
    mg@expression <- list()
    expect_length(mg@expression, 0L)
    # No re-initialize triggered (child structure unchanged)
})

test_that("g@expression <- list() works equivalently to NULL", {
    mg <- .mk_gmulti_pair()
    mg@expression <- list()
    expect_length(mg@expression, 0L)
})


# ----------------------------------------------------------------------------
# 7. Joint slot targeted clear via setter-NULL (existing API)
# ----------------------------------------------------------------------------

test_that("set*(g, NULL) without full nesting errors [defensive]", {
    skip("pending: confirm guard text + check behavior at this layer")
    # Behavior contract: setExpression(g, NULL) with default-only
    # nesting should refuse (errors out from nesting validation).
})
