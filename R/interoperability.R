## sp classes ####



## gef object ####
#' @title Convert gef to Giotto
#' @name gefToGiotto
#' @description Converts .gef file (output stereo-seq pipeline) into
#' giotto subcellular object
#'
#' @param gef_file path to .gef file
#' @param bin_size bin size to select from .gef file
#' @param h5_file name to create and on-disk HDF5 file
#' @param verbose be verbose
#'
#' @details Function in beta. Converts .gef object to Giotto object.
#'
#' There are six possible choices for bin_size:
#' bin1, bin10, bin20, bin50, bin100, bin200.
#'
#' See SAW pipeline for additional information about the gef file.
#' @returns giotto object
#' @export

gefToGiotto <- function(
    gef_file, 
    bin_size = "bin100", 
    verbose = FALSE,
    h5_file = NULL) {
    # data.table vars
    genes <- gene_idx <- x <- y <- sdimx <- sdimy <- cell_ID <- bin_ID <-
        count <- i.bin_ID <- NULL

    # package check
    package_check(pkg_name = "rhdf5", repository = "Bioc")
    if (!file.exists(gef_file)) stop("File path to .gef file does not exist")

    # check if proper bin_size is selected. These are determined in SAW pipeline
    wrap_msg("1. gefToGiotto() begin... \n")
    bin_size_options <- c("bin1", "bin10", "bin20", "bin50", "bin100", "bin200")
    if (!(bin_size %in% bin_size_options)) {
        stop("Please select valid bin size, see ?gefToGiotto for details.")
    }

    # 1. read .gef file at specific bin size
    geneExpData <- rhdf5::h5read(
        file = gef_file,
        name = paste0("geneExp/", bin_size)
    )
    geneDT <- data.table::as.data.table(geneExpData[["gene"]])

    exprDT <- data.table::as.data.table(geneExpData[["expression"]])
    exprDT[, count := lapply(.SD, as.integer), .SDcols = "count"]
    data.table::setorder(exprDT, x, y) # sort by x, y coords (ascending)
    geneDT <- data.table::as.data.table(geneExpData[["gene"]])

    if (isTRUE(verbose)) wrap_msg("finished reading in .gef", bin_size, "\n")

    # 2. create spatial locations
    if (isTRUE(verbose)) wrap_msg("2. create spatial_locations... \n")
    cell_locations <- unique(exprDT[, c("x", "y")], by = c("x", "y"))
    cell_locations[, bin_ID := as.factor(seq_len(nrow(cell_locations)))]
    cell_locations[, cell_ID := paste0("cell_", bin_ID)]
    data.table::setcolorder(cell_locations, c("x", "y", "cell_ID", "bin_ID"))
    # ensure first non-numerical col is cell_ID
    if (isTRUE(verbose)) wrap_msg(nrow(cell_locations), " bins in total \n")
    if (isTRUE(verbose)) wrap_msg("finished spatial_locations \n")

    # 3. create expression matrix
    if (isTRUE(verbose)) wrap_msg("3. create expression matrix... \n")
    exprDT[, genes := as.character(rep(x = geneDT$gene, geneDT$count))]
    exprDT[, gene_idx := as.integer(factor(exprDT$genes,
        levels = unique(exprDT$genes)
    ))]

    # merge on x,y and populate based on bin_ID values in cell_locations
    exprDT[cell_locations, cell_ID := i.bin_ID, on = .(x, y)]
    exprDT$cell_ID <- as.integer(exprDT$cell_ID)

    expMatrix <- Matrix::sparseMatrix(
        i = exprDT$gene_idx,
        j = exprDT$cell_ID,
        x = exprDT$count
    )

    colnames(expMatrix) <- cell_locations$cell_ID
    rownames(expMatrix) <- geneDT$gene
    rm(exprDT)
    if (isTRUE(verbose)) wrap_msg("finished expression matrix")

    # 4. create minimal giotto object
    if (isTRUE(verbose)) wrap_msg("4. create giotto object... \n")
    stereo <- createGiottoObject(
        expression = expMatrix,
        spatial_locs = cell_locations,
        verbose = FALSE,
        h5_file = h5_file
    )
    if (isTRUE(verbose)) wrap_msg("finished giotto object... \n")

    wrap_msg("gefToGiotto() finished \n")
    return(stereo)
}





## Seurat object ####



#' @title Deprecated
#' @name giottoToSeurat
#' @description Deprecated. Please use either [giottoToSeuratV4] or
#' [giottoToSeuratV5]
#' @param gobject Giotto object
#' @param obj_use Giotto object (deprecated, use gobject)
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param ... additional params to pass to \code{\link{get_spatial_locations}}
#' @returns Seurat object
#' @export
giottoToSeurat <- function(
        gobject,
        spat_unit = NULL,
        obj_use = NULL,
        ...) {
    stop(wrap_txt(
        "Deprecated. Please use either giottoToSeuratV4() or giottoToSeuratV5()"
    ))
}




#' @title Convert Giotto to Seurat V4
#' @name giottoToSeuratV4
#' @description Converts Giotto object into a Seurat object. This functions
#' extracts specific sets of data belonging to specified spatial unit.
#' The default values are 'cell' and 'rna' respectively.
#' @param gobject Giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param ... additional params to pass to \code{\link{get_spatial_locations}}
#' @returns Seurat object
#' @keywords seurat interoperability
#' @export
giottoToSeuratV4 <- function(
        gobject,
        spat_unit = NULL,
        ...) {
    # data.table vars
    feat_type <- name <- dim_type <- nn_type <- NULL
    # set default spat_unit and feat_type to be extracted as a Seurat assay
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )
    # verify if optional package is installed
    package_check(pkg_name = "Seurat", repository = "CRAN")
    requireNamespace("Seurat")
    # check whether any raw data exist -- required for Seurat
    avail_expr <- list_expression(gobject = gobject, spat_unit = spat_unit)
    raw_exist <- avail_expr[, "raw" %in% name, by = feat_type]
    # raw_exist <- sapply(gobject@expression_feat,function(x)
    #   'raw' %in% names(gobject@expression[[x]]))
    if (nrow(raw_exist) > 0) {
        assays_all <- raw_exist[, feat_type]
        # assays_all <- names(raw_exist[1])
        # assays_all <- union(assays_all,gobject@expression_feat)
    } else {
        stop("Raw count data not found. Required for Seurat object.")
    }
    # create Seurat object when at least one raw data is available
    for (i in seq_along(assays_all)) {
        assay_use <- assays_all[i]
        expr_use <- lapply(
            avail_expr[feat_type == assay_use, name],
            function(x) {
                get_expression_values(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = assay_use,
                    values = x,
                    output = "exprObj"
                )
            }
        )
        # expr_use <- gobject@expression[[assay_use]]
        names(expr_use) <- unlist(lapply(expr_use, objName))
        slot_use <- names(expr_use)
        if (i == 1) {
            data_raw <- expr_use[["raw"]][]
            sobj <- Seurat::CreateSeuratObject(
                counts = data_raw,
                assay = assay_use
            )
            if ("normalized" %in% slot_use) {
                sobj <- Seurat::SetAssayData(sobj,
                    slot = "data",
                    new.data = expr_use[["normalized"]][],
                    assay = assay_use
                )
            }
            if ("scaled" %in% slot_use) {
                sobj <- Seurat::SetAssayData(sobj,
                    slot = "scale.data",
                    # does not accept 'dgeMatrix'
                    new.data = as.matrix(expr_use[["scaled"]][]),
                    assay = assay_use
                )
            }
        } else {
            if ("raw" %in% slot_use) {
                data_raw <- expr_use[["raw"]][]
                flag_raw <- 1
            } else {
                flag_raw <- 0
            }
            if ("normalized" %in% slot_use) {
                data_norm <- expr_use[["normalized"]][]
                flag_norm <- 1
            } else {
                flag_norm <- 0
            }
            if (flag_raw == 1) {
                assay_obj <- Seurat::CreateAssayObject(counts = data_raw)
            } else if (flag_raw == 0 & flag_norm == 1) {
                assay_obj <- Seurat::CreateAssayObject(data = data_norm)
            } else {
                stop(paste0(
                    "Raw and normalized data not found for assay ",
                    assay_use
                ))
            }
            sobj[[assay_use]] <- assay_obj
            if ("scaled" %in% slot_use) {
                data_scale <- as.matrix(expr_use[["scaled"]][])
                sobj <- Seurat::SetAssayData(sobj,
                    slot = "scale.data",
                    new.data = data_scale,
                    assay = assay_use
                )
            }
        }
        # add cell metadata
        meta_cells <- data.table::setDF(
            get_cell_metadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_cells) <- meta_cells$cell_ID
        meta_cells <- meta_cells[, -which(colnames(meta_cells) == "cell_ID")]
        if (ncol(meta_cells) > 0) {
            colnames(meta_cells) <- paste0(assay_use, "_", colnames(meta_cells))
        }
        sobj <- Seurat::AddMetaData(sobj,
            metadata = meta_cells[Seurat::Cells(sobj), ]
        )
        # add feature metadata
        meta_genes <- data.table::setDF(
            get_feature_metadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_genes) <- meta_genes$feat_ID
        sobj[[assay_use]]@meta.features <- cbind(
            sobj[[assay_use]]@meta.features, meta_genes
        )
        # dim reduction
        # note: Seurat requires assay name specification for each dim reduc
        avail_dr <- list_dim_reductions(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use,
        )
        if (!is.null(avail_dr)) {
            if (nrow(avail_dr) > 0) {
                dr_use <- avail_dr[, name]
                for (i in seq(nrow(avail_dr))) {
                    dr_name <- avail_dr[i, name]
                    dr_type <- avail_dr[i, dim_type]
                    dr_obj <- get_dimReduction(
                        gobject = gobject,
                        output = "dimObj",
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        reduction_method = dr_type,
                        name = dr_name
                    )
                    emb_use <- dr_obj[][Seurat::Cells(sobj), ]
                    if (sum(c("loadings", "eigenvalues") %in%
                        names(slot(dr_obj, "misc"))) == 2) {
                        loadings_use <- slot(dr_obj, "misc")$loadings
                        stdev_use <- slot(dr_obj, "misc")$eigenvalues
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            loadings = loadings_use,
                            key = paste0(dr_name, "_"),
                            stdev = stdev_use,
                            assay = assay_use
                        )
                    } else {
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            key = paste0(dr_name, "_"),
                            assay = assay_use
                        )
                    }
                }
            }
        }
        # network objects
        # expression network
        avail_nn <- list_nearest_networks(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use
        )
        if (!is.null(avail_nn)) {
            if (nrow(avail_nn) > 0) {
                for (i in seq(nrow(avail_nn))) {
                    nn_name <- avail_nn[i, name]
                    nn_type <- avail_nn[i, nn_type]
                    nn_use <- get_NearestNetwork(
                        gobject = gobject,
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        nn_network_to_use = nn_type,
                        network_name = nn_name,
                        output = "data.table"
                    )
                    idx1 <- match(nn_use$from, Seurat::Cells(sobj))
                    idx2 <- match(nn_use$to, Seurat::Cells(sobj))
                    edge_weight <- nn_use$weight
                    edge_dist <- nn_use$distance
                    nn_mtx <- Matrix::sparseMatrix(
                        i = idx1, j = idx2, x = edge_dist,
                        dims = c(ncol(sobj), ncol(sobj))
                    )
                    rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                    nn_name <- paste0(nn_name)
                    sGraph <- Seurat::as.Graph(nn_mtx)
                    sGraph@assay.used <- assay_use
                    sobj[[nn_name]] <- sGraph
                }
            }
        }
    }
    # spatial coordinates
    loc_use <- data.table::setDF(
        get_spatial_locations(
            gobject = gobject,
            spat_unit = spat_unit,
            output = "data.table",
            copy_obj = TRUE,
            ... # allow setting of spat_loc_name through additional params
        )
    )
    rownames(loc_use) <- loc_use$cell_ID
    sobj <- Seurat::AddMetaData(sobj, metadata = loc_use)
    # add spatial coordinates as new dim reduct object
    loc_2 <- loc_use[, c("sdimx", "sdimy")]
    colnames(loc_2) <- c("spatial_1", "spatial_2")
    sobj[["spatial"]] <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(loc_2),
        assay = names(sobj@assays)[1],
        key = "spatial_"
    )
    # spatial network
    avail_sn <- list_spatial_networks(gobject = gobject, spat_unit = spat_unit)
    if (!is.null(avail_sn)) {
        if (nrow(avail_sn) > 0) {
            sn_all <- avail_sn[, name]
            for (i in sn_all) {
                snt_use <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    name = i,
                    output = "networkDT"
                )
                idx1 <- match(snt_use$from, Seurat::Cells(sobj))
                idx2 <- match(snt_use$to, Seurat::Cells(sobj))
                edge_weight <- snt_use$weight
                edge_dist <- snt_use$distance
                nn_mtx <- Matrix::sparseMatrix(
                    i = idx1, j = idx2, x = edge_dist,
                    dims = c(ncol(sobj), ncol(sobj))
                )
                rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                nn_name <- paste0(i)
                spatGraph <- Seurat::as.Graph(nn_mtx)
                spatGraph@assay.used <- names(sobj@assays)[1]
                sobj[[nn_name]] <- spatGraph
            }
        }
    }
    return(sobj)
}


#' @title Convert Giotto to Seurat V5
#' @name giottoToSeuratV5
#' @description Converts Giotto object into a Seurat object. This functions
#' extracts specific sets of data belonging to specified spatial unit.
#' The default values are 'cell' and 'rna' respectively.
#' @param gobject Giotto object
#' @param spat_unit spatial unit (e.g. 'cell')
#' @param res_type type of 10x image output resolution
#' @param ... additional params to pass to \code{\link{get_spatial_locations}}
#' @returns Seurat object
#' @keywords seurat interoperability
#' @export
giottoToSeuratV5 <- function(
        gobject,
        spat_unit = NULL,
        res_type = c("hires", "lowres", "fullres"),
        ...) {
    # data.table vars
    feat_type <- name <- dim_type <- nn_type <- NULL

    res_type <- match.arg(res_type, choices = c("hires", "lowres", "fullres"))
    
    # set default spat_unit and feat_type to be extracted as a Seurat assay
    spat_unit <- set_default_spat_unit(
        gobject = gobject,
        spat_unit = spat_unit
    )

    # verify if optional package is installed
    package_check(pkg_name = "Seurat", repository = "CRAN")
    loadNamespace("Seurat")

    # check whether any raw data exist -- required for Seurat
    avail_expr <- list_expression(gobject = gobject, spat_unit = spat_unit)
    raw_exist <- avail_expr[, "raw" %in% name, by = feat_type]
    # raw_exist <- sapply(gobject@expression_feat,function(x)
    # 'raw' %in% names(gobject@expression[[x]]))
    if (nrow(raw_exist) > 0) {
        assays_all <- raw_exist[, feat_type]
        # assays_all <- names(raw_exist[1])
        # assays_all <- union(assays_all,gobject@expression_feat)
    } else {
        stop("Raw count data not found. Required for Seurat object.")
    }

    # create Seurat object when at least one raw data is available
    for (i in seq_along(assays_all)) {
        assay_use <- assays_all[i]
        expr_use <- lapply(
            avail_expr[feat_type == assay_use, name],
            function(x) {
                get_expression_values(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = assay_use,
                    values = x,
                    output = "exprObj"
                )
            }
        )
        # expr_use <- gobject@expression[[assay_use]]
        names(expr_use) <- unlist(lapply(expr_use, objName))
        slot_use <- names(expr_use)
        if (i == 1) {
            data_raw <- expr_use[["raw"]][]
            sobj <- Seurat::CreateSeuratObject(
                counts = data_raw,
                assay = assay_use
            )
            if ("normalized" %in% slot_use) {
                sobj <- Seurat::SetAssayData(sobj,
                    slot = "data",
                    new.data = expr_use[["normalized"]][],
                    assay = assay_use
                )
            }
            if ("scaled" %in% slot_use) {
                sobj <- Seurat::SetAssayData(sobj,
                    slot = "scale.data",
                    # does not accept 'dgeMatrix'
                    new.data = as.matrix(expr_use[["scaled"]][]),
                    assay = assay_use
                )
            }
        } else {
            if ("raw" %in% slot_use) {
                data_raw <- expr_use[["raw"]][]
                flag_raw <- 1
            } else {
                flag_raw <- 0
            }
            if ("normalized" %in% slot_use) {
                data_norm <- expr_use[["normalized"]][]
                flag_norm <- 1
            } else {
                flag_norm <- 0
            }
            if (flag_raw == 1) {
                assay_obj <- Seurat::CreateAssayObject(counts = data_raw)
            } else if (flag_raw == 0 & flag_norm == 1) {
                assay_obj <- Seurat::CreateAssayObject(data = data_norm)
            } else {
                stop(paste0(
                    "Raw and normalized data not found for assay ",
                    assay_use
                ))
            }
            sobj[[assay_use]] <- assay_obj
            if ("scaled" %in% slot_use) {
                data_scale <- as.matrix(expr_use[["scaled"]][])
                sobj <- Seurat::SetAssayData(sobj,
                    slot = "scale.data",
                    new.data = data_scale,
                    assay = assay_use
                )
            }
        }

        # add cell metadata
        meta_cells <- data.table::setDF(
            get_cell_metadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_cells) <- meta_cells$cell_ID
        meta_cells <- meta_cells[, -which(colnames(meta_cells) == "cell_ID")]
        if (ncol(meta_cells) > 0) {
            colnames(meta_cells) <- paste0(
                assay_use, "_",
                colnames(meta_cells)
            )
        }
        sobj <- Seurat::AddMetaData(sobj,
            metadata = meta_cells[Seurat::Cells(sobj), ]
        )

        # add feature metadata
        meta_genes <- data.table::setDF(
            get_feature_metadata(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = assay_use,
                output = "data.table",
                copy_obj = TRUE
            )
        )
        rownames(meta_genes) <- meta_genes$feat_ID
        for (i in seq(sobj@assays)) {
            sobj@assays[[i]]@meta.data <- meta_genes
        }


        # dim reduction
        # note: Seurat requires assay name specification for each dim reduc
        avail_dr <- list_dim_reductions(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use,
        )
        if (!is.null(avail_dr)) {
            if (nrow(avail_dr) > 0) {
                dr_use <- avail_dr[, name]
                for (i in seq(nrow(avail_dr))) {
                    dr_name <- avail_dr[i, name]
                    dr_type <- avail_dr[i, dim_type]
                    dr_obj <- get_dimReduction(
                        gobject = gobject,
                        output = "dimObj",
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        reduction_method = dr_type,
                        name = dr_name
                    )
                    emb_use <- dr_obj[][Seurat::Cells(sobj), ]
                    if (sum(c("loadings", "eigenvalues") %in%
                        names(slot(dr_obj, "misc"))) == 2) {
                        loadings_use <- slot(dr_obj, "misc")$loadings
                        stdev_use <- slot(dr_obj, "misc")$eigenvalues
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            loadings = loadings_use,
                            key = paste0(dr_name, "_"),
                            stdev = stdev_use,
                            assay = assay_use
                        )
                    } else {
                        sobj[[dr_name]] <- Seurat::CreateDimReducObject(
                            embeddings = as.matrix(emb_use),
                            key = paste0(dr_name, "_"),
                            assay = assay_use
                        )
                    }
                }
            }
        }



        # network objects
        # expression network
        avail_nn <- list_nearest_networks(
            gobject = gobject,
            spat_unit = spat_unit,
            feat_type = assay_use
        )
        if (!is.null(avail_nn)) {
            if (nrow(avail_nn) > 0) {
                for (i in seq(nrow(avail_nn))) {
                    nn_name <- avail_nn[i, name]
                    nn_type <- avail_nn[i, nn_type]
                    nn_use <- get_NearestNetwork(
                        gobject = gobject,
                        spat_unit = spat_unit,
                        feat_type = assay_use,
                        nn_network_to_use = nn_type,
                        network_name = nn_name,
                        output = "data.table"
                    )
                    idx1 <- match(nn_use$from, Seurat::Cells(sobj))
                    idx2 <- match(nn_use$to, Seurat::Cells(sobj))
                    edge_weight <- nn_use$weight
                    edge_dist <- nn_use$distance
                    nn_mtx <- Matrix::sparseMatrix(
                        i = idx1, j = idx2, x = edge_dist,
                        dims = c(ncol(sobj), ncol(sobj))
                    )
                    rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                    nn_name <- paste0(nn_name)
                    sGraph <- Seurat::as.Graph(nn_mtx)
                    sGraph@assay.used <- assay_use
                    sobj[[nn_name]] <- sGraph
                }
            }
        }
    }

    # spatial coordinates
    loc_use <- getSpatialLocations(
        gobject = gobject,
        spat_unit = spat_unit,
        output = "spatLocsObj",
        copy_obj = TRUE,
        ... # allow setting of spat_loc_name through additional params
    )
    
    # flip y vals
     	
      loc_use <- flip(loc_use)[] %>%
      data.table::setDF()
      
    rownames(loc_use) <- loc_use$cell_ID
    sobj <- Seurat::AddMetaData(sobj, metadata = loc_use)
    # add spatial coordinates as new dim reduct object
    loc_2 <- loc_use[, c("sdimx", "sdimy")]
    colnames(loc_2) <- c("spatial_1", "spatial_2")
    sobj[["spatial"]] <- Seurat::CreateDimReducObject(
        embeddings = as.matrix(loc_2),
        assay = names(sobj@assays)[1],
        key = "spatial_"
    )



    # spatial network
    avail_sn <- list_spatial_networks(gobject = gobject, spat_unit = spat_unit)
    if (!is.null(avail_nn)) {
        if (nrow(avail_sn) > 0) {
            sn_all <- avail_sn[, name]
            for (i in sn_all) {
                snt_use <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    name = i,
                    output = "networkDT"
                )
                idx1 <- match(snt_use$from, Seurat::Cells(sobj))
                idx2 <- match(snt_use$to, Seurat::Cells(sobj))
                edge_weight <- snt_use$weight
                edge_dist <- snt_use$distance
                nn_mtx <- Matrix::sparseMatrix(
                    i = idx1, j = idx2, x = edge_dist,
                    dims = c(ncol(sobj), ncol(sobj))
                )
                rownames(nn_mtx) <- colnames(nn_mtx) <- Seurat::Cells(sobj)
                nn_name <- paste0(i)
                spatGraph <- Seurat::as.Graph(nn_mtx)
                spatGraph@assay.used <- names(sobj@assays)[1]
                sobj[[nn_name]] <- spatGraph
            }
        }
    }

    all_x <- NULL
    all_y <- NULL

    gimgs <- getGiottoImage(gobject, name = ":all:")
    if (length(gimgs) > 0) {
      for (i in seq_along(gimgs)) {
        gimg <- gimgs[[i]]
        key <- objName(gimg)
        imagerow <- loc_use$sdimy
        imagecol <- loc_use$sdimx
        img_array <- as(gimg, "array")
        img_array <- img_array / 255	
        coord <- data.frame(
          imagerow = imagerow, imagecol = imagecol, 
                    row.names = loc_use$cell_ID
                    )
        scalef <- .estimate_scalefactors(
          gimg,
          res_type = res_type,
          spatlocs = loc_use
          )
        # There does not seem to be a way to tell seurat which image type
        # you are using. The lowres scalefactor seems to be the important
        # one in mapping the image
            scalefactors <- Seurat::scalefactors(
              spot = scalef$spot,
                fiducial = scalef$fiducial,
                hires = scalef$hires,
                lowres = scalef[res_type] # this looks like the main one
                # so instead of strictly supplying lowres scalef, we use
                # the scalef belonging to whichever image was used in Giotto
                # since we allow use non-lowres images
            )
            # see https://github.com/satijalab/seurat/issues/3595
            newV1 <- new(
                Class = "VisiumV1",
                image = img_array,
                scale.factors = scalefactors,
                coordinates = coord, 
                spot.radius = 
                  scalef$fiducial * scalef$lowres / max(dim(img_array)),
                key = paste0(key, "_")
            )

            sobj@images[[key]] <- newV1
        }
    }



    return(sobj)
}

#' @param x image object
#' @param res_type type of 10x image output resolution
#' @param spatlocs a data.frame of spatial locations coordinates
#' @noRd

.estimate_scalefactors <- function(
    x, 
    res_type = c("hires", "lowres", "fullres"),
    spatlocs
    ){
  res_type <- match.arg(res_type, 
                        choices = c("hires", "lowres", "fullres"))
        pxdims <- dim(x)[1:2]
        edims <- range(ext(x))
        scalef <- mean(pxdims / edims)
# assume that lowres and hires follow a general ratio
# may not be that important since the scalefactor should theoretically
# only matter for the image res that we are using
# this ratio is roughly 3.333334 based on Visium BreastCancerA1 dataset
        res_ratio <- 3.333334
        # fullres should have a scalef of roughly 1.
        # No way to guess hires or lowres scalefs so use arbitrary values.
        hres_scalef <- switch(res_type,
                              "hires" = scalef,
                               "lowres" = scalef * res_ratio, 
                               "fullres" = 0.08250825 # arbitrary
                   	
                )
        lres_scalef <- switch(res_type,
                              "hires" = scalef / res_ratio,
                               "lowres" = scalef,
                               "fullres" = 0.02475247 # arbitrary
                 )
# spot diameter and fid diameter are variable based on how spatial info was
# mapped to the image. Estimate this by getting the center to center
 # px distance vs fullsize px dims ratio.
 # ! fullsize px dims is the same as edims !
            coords <- data.table::as.data.table(spatlocs)
             # create a delaunay
             dnet <- createNetwork(
             as.matrix(coords[, c("sdimx", "sdimy")]), 
             type = "delaunay",
             method = "geometry",
             include_distance = TRUE, 
             as.igraph = FALSE, 
             include_weight = TRUE, 
             verbose = FALSE
             )
             	
 # expect center to center be most common edge distance
 # this gives CC dist as fullres px distance
  distances <- sort(unique(dnet$distance))
  cc_px <- distances[which.max(table(dnet$distance))]
 # assume constant ratios between diameters and cc_px
 fid_cc_ratio <- 1.045909
fid_diam <- cc_px * fid_cc_ratio
spot_cc_ratio <- 0.6474675
spot_diam <- cc_px * spot_cc_ratio
scalef_list <- list(
  spot = spot_diam,
  fiducial = fid_diam,
  hires = hres_scalef,
  lowres = lres_scalef
  )
return(scalef_list)
                                  
  }


#' @title Deprecated
#' @name seuratToGiotto
#' @description Deprecated. Please use either [seuratToGiottoV4] or
#' [seuratToGiottoV5]
#' @param sobject Input Seurat object to convert to Giotto object
#' @param spatial_assay Specify name of the spatial assay slot in Seurat.
#' Default is \code{"Spatial"}.
#' @param dim_reduction Specify which dimensional reduction computations to
#' fetch from
#' @param subcellular_assay Specify name of the subcellular assay in input
#' object. Default is \code{"Vizgen"}.
#' @returns giotto object
#' @export
seuratToGiotto <- function(
        sobject,
        spatial_assay = "Spatial",
        dim_reduction = c("pca", "umap"),
        subcellular_assay = "Vizgen") {
    stop(wrap_txt(
        "Deprecated. Please use either seuratToGiottoV4() or seuratToGiottoV5()"
    ))
}




#' @title Convert a Seurat V4 object to a Giotto object
#' @name seuratToGiottoV4
#'
#' @param sobject Input Seurat object to convert to Giotto object
#' @param spatial_assay Specify name of the spatial assay slot in Seurat.
#' Default is \code{"Spatial"}.
#' @param dim_reduction Specify which dimensional reduction computations to
#' fetch from input Seurat object. Default is \code{"c('pca', 'umap')"}.
#' @param subcellular_assay Specify name of the subcellular assay in input
#' @param sp_network sp_network
#' @param nn_network nn_network
#' @param verbose logical. Default to TRUE
#' object. Default is \code{"Vizgen"}.
#' @returns A Giotto object converted from Seurat object with all computations
#' stored in it.
#' @keywords seurat interoperability
#' @examples
#' m_expression <- Matrix::Matrix(rnorm(100), nrow = 10, sparse = TRUE)
#' s <- Seurat::CreateSeuratObject(counts = m_expression)
#'
#' seuratToGiottoV5(s, spatial_assay = "RNA")
#' @export
seuratToGiottoV4 <- function(
        sobject,
        spatial_assay = "Spatial",
        dim_reduction = c("pca", "umap"),
        subcellular_assay = "Vizgen",
        sp_network = NULL,
        nn_network = NULL,
        verbose = TRUE) {
    package_check("Seurat")
    if (is.null(Seurat::GetAssayData(
        object = sobject, slot = "counts",
        assay = spatial_assay
    ))) {
        wrap_msg("No raw expression values are provided in spatial_assay")
        return(sobject)
    } else {
        exp <- Seurat::GetAssayData(
            object = sobject, slot = "counts",
            assay = spatial_assay
        )
        if (!is.null(sobject@assays$SCT)) {
            normexp <- Seurat::GetAssayData(
                object = sobject, slot = "counts",
                assay = "SCT"
            )
        }
        if (!is.null(slot(sobject, "assays")[[spatial_assay]]@data)) {
            normexp <- Seurat::GetAssayData(
                object = sobject, slot = "data",
                assay = spatial_assay
            )
        }
        # Cell Metadata
        cell_metadata <- sobject@meta.data
        # Dimension Reduction
        if (sum(vapply(
            dim_reduction,
            function(x) length(sobject@reductions[[x]]),
            FUN.VALUE = integer(1L)
        )) == 0) {
            dimReduc_list <- NULL
        } else {
            dimReduc_list <- lapply(dim_reduction, function(x) {
                dim_coord <- as.matrix(Seurat::Embeddings(
                    object = sobject,
                    reduction = x
                ))
                dim_load <- as.matrix(Seurat::Loadings(
                    object = sobject,
                    reduction = x
                ))
                dim_eig <- Seurat::Stdev(sobject, reduction = x)
                if (length(dim_eig) > 0) {
                    dim_eig <- sapply(dim_eig, function(x) x^2)
                }
                colnames(dim_coord) <- paste0("Dim.", seq_len(ncol(dim_coord)))
                if (length(dim_load) > 0) {
                    colnames(dim_load) <- paste0(
                        "Dim.", seq_len(ncol(dim_load))
                    )
                }
                dimReduc <- create_dim_obj(
                    name = x,
                    reduction = "cells",
                    reduction_method = x,
                    spat_unit = "cell",
                    feat_type = "rna",
                    provenance = NULL,
                    coordinates = dim_coord,
                    misc = list(eigenvalues = dim_eig, loadings = dim_load)
                )
                return(dimReduc)
            })
            # names(dimReduc_list) <- dim_reduction
        }
        # Spatial Locations
        if (length(sobject@assays[[spatial_assay]]) == 0) {
            spat_loc <- NULL
        } else {
            # !requires image objects!
            if (!is.null(Seurat::Images(
                object = sobject,
                assay = spatial_assay
            ))) {
                spat_coord <- Seurat::GetTissueCoordinates(sobject)
                spat_coord <- cbind(
                    rownames(spat_coord),
                    data.frame(spat_coord, row.names = NULL)
                )
                colnames(spat_coord) <- c("cell_ID", "sdimy", "sdimx")
                spat_loc <- spat_coord
            } else {
                message("Images for RNA assay not found in the data.
                        Skipping image processing.")
                spat_loc <- NULL
            }
        }


        # Subcellular
        name <- names(sobject@images)
        if (length(sobject@assays[[subcellular_assay]]) == 1) {
            spat_coord <- Seurat::GetTissueCoordinates(sobject)
            colnames(spat_coord) <- c("sdimx", "sdimy", "cell_ID")
            exp <- exp[, c(intersect(spat_coord$cell_ID, colnames(exp)))]
            spat_loc <- spat_coord
        }
        if (!length(sobject@images) == 0) {
            if ("molecules" %in% methods::slotNames(
                sobject@images[[name]]
            ) == TRUE) {
                if (!length(sobject@images[[name]][["molecules"]]) == 0) {
                    assay <- names(sobject@assays)
                    featnames <- rownames(sobject@assays[[assay]]@meta.features)
                    mol_spatlocs <- data.table::data.table()
                    for (x in featnames) {
                        df <- (Seurat::FetchData(
                            sobject[[name]][["molecules"]],
                            vars = x
                        ))
                        mol_spatlocs <- rbind(mol_spatlocs, df)
                    }
                    gpoints <- createGiottoPoints(mol_spatlocs,
                        feat_type = "rna"
                    )
                }
            }
        }
    }

    # Networks

    # spatial_network
    weight <- distance <- NULL

    if (!is.null(sp_network)) {
        for (i in seq(sp_network)) {
            if (verbose) message("Copying spatial networks")
            stempgraph <- Seurat::as.Graph(sobject@graphs[[sp_network[i]]])
            m <- as.matrix(stempgraph)
            sobjIgraph <- igraph::graph.adjacency(m, weighted = TRUE)

            if (verbose) message("Calculating distances")
            for (j in seq_along(sobject@graphs[[sp_network[i]]]@x)) {
                dist_matrix <- sobject@graphs[[sp_network[i]]]@x[[j]]
                sobjIgraph <- igraph::set_edge_attr(
                    graph = sobjIgraph, name = "distance",
                    value = dist_matrix, index = j
                )
            }
            e_attr <- igraph::edge_attr_names(sobjIgraph)
            if ("weight" %in% e_attr) {
                igDT <- data.table::setDT(igraph::as_data_frame(sobjIgraph))
                igDT[, weight := 1 / (1 + distance)]
                data.table::setcolorder(igDT,
                    neworder = c("from", "to", "weight", "distance")
                )
                sobjIgraph <- igraph::graph_from_data_frame(igDT)
            }
            DT <- data.table()

            edges <- igraph::ends(sobjIgraph,
                es = igraph::E(sobjIgraph),
                names = TRUE
            )
            DT$from <- edges[, 1]
            DT$to <- edges[, 2]
            ed_attr <- igraph::edge.attributes(sobjIgraph)
            DT$weight <- ed_attr[1]
            DT$distance <- ed_attr[2]
            spatNetObj <- create_spat_net_obj(
                networkDT = DT
            )
            gobject <- set_spatialNetwork(
                gobject = gobject,
                spatial_network = spatNetObj,
                name = sp_network[i]
            )
        }
    }
    ## nn_network
    if (!is.null(nn_network)) {
        nnNetObj_list <- lapply(seq_along(nn_network), function(i) {
            stempgraph <- Seurat::as.Graph(sobject@graphs[[nn_network[i]]])
            m <- as.matrix(stempgraph)
            sobjIgraph <- igraph::graph.adjacency(m, weighted = TRUE)

            if (verbose) message("Calculating distances")
            for (j in seq_along(sobject@graphs[[nn_network[i]]]@x)) {
                dist_matrix <- sobject@graphs[[nn_network[i]]]@x[[j]]

                sobjIgraph <- igraph::set_edge_attr(
                    graph = sobjIgraph, name = "distance",
                    value = dist_matrix, index = j
                )
            }
            e_attr <- igraph::edge_attr_names(sobjIgraph)
            if ("weight" %in% e_attr) {
                igDT <- data.table::setDT(igraph::as_data_frame(sobjIgraph))
                igDT[, weight := 1 / (1 + distance)]
                data.table::setcolorder(igDT,
                    neworder = c("from", "to", "weight", "distance")
                )
                sobjIgraph <- igraph::graph_from_data_frame(igDT)
            }

            if (verbose) message("Copying nearest neighbour networks")
            nnNetObj <- GiottoClass::createNearestNetObj(
                name = names(sobject@graphs)[i],
                network = sobjIgraph
            )
            return(nnNetObj)
        })

        for (i in seq_along(nnNetObj_list)) {
            gobject <- GiottoClass::set_NearestNetwork(
                gobject = gobject,
                nn_network = nnNetObj_list[[i]]
            )
        }
    }
    gobject <- createGiottoObject(exp,
        spatial_locs = spat_loc,
        dimension_reduction = dimReduc_list
    )
    if (exists("normexp") == TRUE) {
        exprObj <- create_expr_obj(
            name = "normalized",
            exprMat = normexp,
            spat_unit = "cell",
            feat_type = "rna",
            provenance = "cell"
        )
        gobject <- set_expression_values(
            gobject = gobject,
            values = exprObj, set_defaults = FALSE
        )
        # gobject@expression$cell$rna$normalized = normexp
    }
    gobject <- addCellMetadata(gobject = gobject, new_metadata = cell_metadata)
    if (exists("gpoints") == TRUE) {
        gobject <- addGiottoPoints(
            gobject = gobject,
            gpoints = list(gpoints)
        )
    }
    return(gobject)
}

#' @title Convert a Seurat V5 object to a Giotto object
#' @name seuratToGiottoV5
#'
#' @param sobject Input Seurat object to convert to Giotto object
#' @param spatial_assay Specify name of the spatial assay slot in Seurat.
#' Default is \code{"Spatial"}.
#' @param dim_reduction Specify which dimensional reduction computations to
#' fetch from input Seurat object. Default is \code{"c('pca', 'umap')"}.
#' @param subcellular_assay Specify name of the subcellular assay in input
#' @param sp_network sp_network
#' @param nn_network nn_network
#' @param verbose logical. Default to TRUE
#' object. Default is \code{"Vizgen"}.
#' @returns A Giotto object converted from Seurat object with all computations
#' stored in it.
#' @keywords seurat interoperability
#' @export
seuratToGiottoV5 <- function(
        sobject,
        spatial_assay = "Spatial",
        dim_reduction = c("pca", "umap"),
        subcellular_assay = "Vizgen",
        sp_network = NULL,
        nn_network = NULL,
        verbose = TRUE) {
    package_check("Seurat")
 
   # NSE vars
    sdimy <- NULL
    
    if (is.null(Seurat::GetAssayData(
        object = sobject, slot = "counts",
        assay = spatial_assay
    ))) {
        wrap_msg("No raw expression values are provided in spatial_assay")
        return(sobject)
    } else {
        exp <- Seurat::GetAssayData(
            object = sobject, slot = "counts",
            assay = spatial_assay
        )
        if (!is.null(sobject@assays$SCT)) {
            normexp <- Seurat::GetAssayData(
                object = sobject, slot = "counts",
                assay = "SCT"
            )
        }

        if ("data" %in% slotNames(sobject@assays[[spatial_assay]])) {
            if (!is.null(slot(sobject, "assays")[[spatial_assay]]@data)) {
                normexp <- Seurat::GetAssayData(
                    object = sobject,
                    slot = "data", assay = spatial_assay
                )
            }
        }
        if ("layers" %in% slotNames(sobject@assays[[spatial_assay]])) {
            if (!is.null(slot(sobject, "assays")[[spatial_assay]]@layers)) {
                normexp <- SeuratObject::LayerData(
                    object = sobject,
                    assay = spatial_assay
                )
            }
        }

        # Cell Metadata
        cell_metadata <- sobject@meta.data
        cell_metadata <- data.table::as.data.table(
        cell_metadata, keep.rownames = TRUE)
        
        # Feat Metadata
        feat_metadata <- sobject[[]]
        feat_metadata <- data.table::as.data.table(
        feat_metadata, keep.rownames = TRUE)

        # rownames of both kept as `rn`
        # Dimension Reduction
        if (sum(vapply(
            dim_reduction,
            function(x) length(sobject@reductions[[x]]),
            FUN.VALUE = integer(1L)
        )) == 0) {
            dimReduc_list <- NULL
        } else {
            dimReduc_list <- lapply(dim_reduction, function(x) {
                dim_coord <- as.matrix(Seurat::Embeddings(
                    object = sobject,
                    reduction = x
                ))
                dim_load <- as.matrix(Seurat::Loadings(
                    object = sobject,
                    reduction = x
                ))
                dim_eig <- Seurat::Stdev(sobject, reduction = x)
                if (length(dim_eig) > 0) {
                    dim_eig <- sapply(dim_eig, function(x) x^2)
                }
                colnames(dim_coord) <- paste0("Dim.", seq_len(ncol(dim_coord)))
                if (length(dim_load) > 0) {
                    colnames(dim_load) <- paste0(
                        "Dim.", seq_len(ncol(dim_load))
                    )
                }
                dimReduc <- create_dim_obj(
                    name = x,
                    reduction = "cells",
                    reduction_method = x,
                    spat_unit = "cell",
                    feat_type = "rna",
                    provenance = NULL,
                    coordinates = dim_coord,
                    misc = list(eigenvalues = dim_eig, loadings = dim_load)
                )

                return(dimReduc)
            })
            # names(dimReduc_list) <- dim_reduction
        }

        # Spatial Locations
        if (length(sobject@assays[[spatial_assay]]) == 0) {
            spat_loc <- NULL
        } else {
            # !requires image objects!
            if (!is.null(Seurat::Images(
                object = sobject,
                assay = spatial_assay
            ))) {
                spat_coord <- Seurat::GetTissueCoordinates(
                  sobject,
                  scale = NULL,
                  cols = c(
                  "imagerow",
                  "imagecol"))

                if (!("cell" %in% spat_coord)) {
                    spat_coord$cell_ID <- rownames(spat_coord)
                    colnames(spat_coord) <- c("sdimy", "sdimx", "cell_ID")
                } else {
                    colnames(spat_coord) <- c("sdimy", "sdimx", "cell_ID")
                }

                spat_loc <- data.table::as.data.table(spat_coord)
 
                # seurat has coords following imaging conventions
                # flip them for Giotto
                spat_loc[, sdimy := -sdimy]
                data.table::setcolorder(
                  spat_loc,
                  neworder = c("sdimx",
                               "sdimy", 
                               "cell_ID"
                               )
                  )
              
            } else {
                message("Images for RNA assay not found in the data.
                        Skipping image processing.")
                spat_loc <- NULL
            }
        }
        # Subcellular
        name <- names(sobject@images)
        # if (!is.null(subcellular_assay)){
        if (length(sobject@assays[[subcellular_assay]]) == 1) {
            spat_coord <- Seurat::GetTissueCoordinates(sobject)
            colnames(spat_coord) <- c("sdimx", "sdimy", "cell_ID")
            exp <- exp[, c(intersect(spat_coord$cell_ID, colnames(exp)))]
            spat_loc <- spat_coord
        }
        # }
        if (!length(sobject@images) == 0) {
            for (i in names(sobject@images)) {
                if ("molecules" %in% names(sobject@images[[i]]) == TRUE) {
                    if (!length(sobject@images[[i]][["molecules"]]) == 0) {
                        assay <- names(sobject@assays)
                        featnames <- rownames(sobject)
                        mol_spatlocs <- data.table::data.table()

                        for (x in featnames) {
                            df <- (Seurat::FetchData(
                                sobject[[i]][["molecules"]],
                                vars = x
                            ))
                            mol_spatlocs <- rbind(mol_spatlocs, df)
                        }
                        gpoints <- createGiottoPoints(mol_spatlocs,
                            feat_type = "rna"
                        )
                        if ("centroids" %in% names(sobject@images[[i]])) {
                            centroids_coords <-
                                sobject@images[[i]]$centroids@coords
                            centroids_coords <- vect(centroids_coords)
                            gpolygon <- create_giotto_polygon_object(
                                name = "cell", spatVector = centroids_coords
                            )
                        }
                        if ("segmentation" %in% names(sobject@images[[i]])) {
                            polygon_list <- list()

                            for (j in seq(sobject@images[[
                                i
                            ]]@boundaries$segmentation@polygons)) {
                                polygon_info <- sobject@images[[
                                    i
                                ]]@boundaries$segmentation@polygons[[j]]

                                # Get coordinates from segmentation
                                seg_coords <- polygon_info@Polygons[[1]]@coords

                                # Fetch cell_Id from polygon information
                                cell_ID <- polygon_info@ID

                                # Convert it to SpatVector
                                seg_coords <- vect(seg_coords)

                                # Create giotto_polygon_object
                                gpolygon <- create_giotto_polygon_object(
                                    name = "cell",
                                    spatVector = centroids_coords,
                                    spatVectorCentroids = seg_coords
                                )

                                # Add the cell_ID to the list of polygon names
                                polygon_list[[cell_ID]] <- gpolygon
                            }
                        }
                    }
                }
            }
        }
    }

    # Find SueratImages, extract them, and pass to create image

    for (i in names(sobject@images)) {
      simg <- sobject[[i]]
        # check if image slot has image in it
        if ("image" %in% slotNames(simg)) {
          img_array <- slot(simg, "image")
          if (!is.null(img_array)) { 
            scalef <- Seurat::ScaleFactors(simg)
                gImg <- createGiottoLargeImage(
                    raster_object = terra::rast(img_array) * 255,
                    name = i,
                    scale_factor = 1 / scalef$lowres
                )
            }
        }
    }


    gobject <- createGiottoObject(exp,
        spatial_locs = spat_loc,
        dimension_reduction = dimReduc_list
    )
    if (exists("normexp") == TRUE) {
        exprObj <- create_expr_obj(
            name = "normalized",
            exprMat = normexp,
            spat_unit = "cell",
            feat_type = "rna",
            provenance = "cell"
        )
        gobject <- set_expression_values(
            gobject = gobject,
            values = exprObj, set_defaults = FALSE
        )
        # gobject@expression$cell$rna$normalized = normexp
    }

    # Networks

    # spatial_network
    weight <- distance <- NULL

    if (!is.null(sp_network)) {
        for (i in seq(sp_network)) {
            if (verbose) message("Copying spatial networks")
            stempgraph <- Seurat::as.Graph(sobject@graphs[[sp_network[i]]])
            m <- as.matrix(stempgraph)
            sobjIgraph <- igraph::graph.adjacency(m, weighted = TRUE)

            if (verbose) message("Calculating distances")
            for (j in seq_along(sobject@graphs[[sp_network[i]]]@x)) {
                dist_matrix <- sobject@graphs[[sp_network[i]]]@x[[j]]
                sobjIgraph <- igraph::set_edge_attr(
                    graph = sobjIgraph, name = "distance",
                    value = dist_matrix, index = j
                )
            }
            e_attr <- igraph::edge_attr_names(sobjIgraph)
            if ("weight" %in% e_attr) {
                igDT <- data.table::setDT(igraph::as_data_frame(sobjIgraph))
                igDT[, weight := 1 / (1 + distance)]
                data.table::setcolorder(igDT,
                    neworder = c("from", "to", "weight", "distance")
                )
                sobjIgraph <- igraph::graph_from_data_frame(igDT)
            }
            DT <- data.table()

            edges <- igraph::ends(sobjIgraph,
                es = igraph::E(sobjIgraph),
                names = TRUE
            )
            DT$from <- edges[, 1]
            DT$to <- edges[, 2]
            ed_attr <- igraph::edge.attributes(sobjIgraph)
            DT$weight <- ed_attr[1]
            DT$distance <- ed_attr[2]
            spatNetObj <- create_spat_net_obj(
                networkDT = DT
            )
            gobject <- set_spatialNetwork(
                gobject = gobject,
                spatial_network = spatNetObj,
                name = sp_network[i]
            )
        }
    }

    ## nn_network
    if (!is.null(nn_network)) {
        nnNetObj_list <- lapply(seq_along(nn_network), function(i) {
            stempgraph <- Seurat::as.Graph(sobject@graphs[[nn_network[i]]])
            m <- as.matrix(stempgraph)
            sobjIgraph <- igraph::graph.adjacency(m, weighted = TRUE)

            if (verbose) message("Calculating distances")
            for (j in seq_along(sobject@graphs[[nn_network[i]]]@x)) {
                dist_matrix <- sobject@graphs[[nn_network[i]]]@x[[j]]

                sobjIgraph <- igraph::set_edge_attr(
                    graph = sobjIgraph, name = "distance",
                    value = dist_matrix, index = j
                )
            }
            e_attr <- igraph::edge_attr_names(sobjIgraph)
            if ("weight" %in% e_attr) {
                igDT <- data.table::setDT(igraph::as_data_frame(sobjIgraph))
                igDT[, weight := 1 / (1 + distance)]
                data.table::setcolorder(igDT,
                    neworder = c("from", "to", "weight", "distance")
                )
                sobjIgraph <- igraph::graph_from_data_frame(igDT)
            }

            if (verbose) message("Copying nearest neighbour networks")
            nnNetObj <- GiottoClass::createNearestNetObj(
                name = names(sobject@graphs)[i],
                network = sobjIgraph
            )
            return(nnNetObj)
        })

        for (i in seq_along(nnNetObj_list)) {
            gobject <- GiottoClass::set_NearestNetwork(
                gobject = gobject,
                nn_network = nnNetObj_list[[i]]
            )
        }
    }
    gobject <- addCellMetadata(
      gobject = gobject, new_metadata = cell_metadata, 
      by_column = TRUE, column_cell_ID = "rn")
    
    gobject <- addFeatMetadata(
      gobject = gobject, new_metadata = feat_metadata, 
      by_column = TRUE, column_feat_ID = "rn")

    if (exists("gpoints")) {
        gobject <- addGiottoPoints(
            gobject = gobject,
            gpoints = list(gpoints)
        )
    }

    if (exists("gpolygon")) {
        gobject <- addGiottoPolygons(
            gobject = gobject,
            gpolygons = polygon_list
        )
    }

    if (exists("gImg")) {
        gobject <- addGiottoLargeImage(
            gobject = gobject,
            largeImages = list(gImg)
        )
    }
    return(gobject)
}






## SpatialExperiment object ####

#' Utility function to convert a Giotto object to a SpatialExperiment object.
#'
#' @param giottoObj Input Giotto object to convert to a SpatialExperiment object
#' @param verbose A boolean value specifying if progress messages should be
#' displayed or not. Default \code{TRUE}.
#'
#' @returns A SpatialExperiment object that contains data from the input Giotto
#' object.
#' @examples
#' \dontrun{
#' mini_gobject <- GiottoData::loadGiottoMini("vizgen")
#' giottoToSpatialExperiment(mini_gobject)
#' }
#' @export
giottoToSpatialExperiment <- function(giottoObj, verbose = TRUE) {
    spat_unit <- NULL

    # Load required packages
    # package_check(pkg_name = "SummarizedExperiment", repository = 'Bioc')
    # SP should load this?
    # package_check(pkg_name = "SingleCellExperiment", repository = 'Bioc')
    # SP should load this?
    package_check(pkg_name = "SpatialExperiment", repository = "Bioc")
    package_check(pkg_name = "S4Vectors", repository = "Bioc")

    # SpatialExperiment objects list (one object for each spatial unit)
    speList <- list()

    # Expression Matrices
    giottoExpr <- list_expression(giottoObj)

    # Iterate over spatial units
    spatialUnits <- unique(giottoExpr$spat_unit) # a function to get spat units?
    for (su in seq(spatialUnits)) {
        if (verbose) {
            message("Processing spatial unit: '", spatialUnits[su], "'")
        }
        # Check if expression matrices exist in input object
        if (!is.null(giottoExpr)) {
            if (verbose) {
                message(
                    "Copying expression matrix: '", giottoExpr[1]$name,
                    "' for spatial unit: '", spatialUnits[su], "'"
                )
            }
            exprMat <- get_expression_values(
                gobject = giottoObj,
                spat_unit = spatialUnits[su],
                feat_type = giottoExpr[1]$feat_type,
                values = giottoExpr[1]$name,
                output = "matrix"
            )
            names(rownames(exprMat)) <- NULL
            names(colnames(exprMat)) <- NULL
            exprMat <- list(exprMat)
            names(exprMat)[1] <- giottoExpr[1]$name
            # Creating SPE object with first expression matrix
            spe <- SpatialExperiment::SpatialExperiment(assays = exprMat)
            SummarizedExperiment::assayNames(spe) <- paste0(
                giottoExpr[1]$name, "_",
                giottoExpr[1]$feat_type, "_",
                spatialUnits[su]
            )
            giottoExpr <- giottoExpr[-1, ]
        } else {
            stop("The input Giotto object must contain atleast one expression
                matrix.")
        }

        # Copying remaining expression matrices if they exist
        if (nrow(giottoExpr[spat_unit == spatialUnits[su]]) > 0) {
            for (i in seq(nrow(giottoExpr))) {
                if (verbose) {
                    message(
                        "Copying expression matrix: '",
                        giottoExpr[i]$name, "' for spatial unit: '",
                        spatialUnits[su], "'"
                    )
                }
                # SPE does not have specific slots for different units,
                # instead joining multiple unit names to identify them
                SummarizedExperiment::assay(
                    spe,
                    paste0(
                        giottoExpr[i]$name, "_",
                        giottoExpr[i]$feat_type, "_",
                        spatialUnits[su]
                    ),
                    withDimnames = FALSE
                ) <- get_expression_values(
                    gobject = giottoObj,
                    spat_unit = spatialUnits[su],
                    feat_type = giottoExpr[i]$feat_type,
                    values = giottoExpr[i]$name,
                    output = "matrix"
                )
            }
        }

        # Cell Metadata to ColData
        pData <- pDataDT(gobject = giottoObj, spat_unit = spatialUnits[su])
        if (nrow(pData) > 0) {
            if (verbose) {
                message(
                    "Copying phenotype data for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            SummarizedExperiment::colData(spe) <- S4Vectors::DataFrame(
                pData,
                row.names = pData$cell_ID
            )
        } else {
            if (verbose) {
                message("No phenotype data found in input Giotto object")
            }
        }

        # Feature Metadata to RowData
        fData <- fDataDT(gobject = giottoObj, spat_unit = spatialUnits[su])
        if (nrow(fData) > 0) {
            if (verbose) {
                message(
                    "Copying feature metadata for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            SummarizedExperiment::rowData(spe) <- fData
        } else {
            if (verbose) {
                message("No feature metadata found in input Giotto object")
            }
        }

        # Spatial Locations to Spatial Coordinates
        spatialLocs <- get_spatial_locations(
            gobject = giottoObj,
            spat_unit = spatialUnits[su],
            output = "data.table"
        )
        if (!is.null(spatialLocs)) {
            if (verbose) {
                message(
                    "Copying spatial locations for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            SpatialExperiment::spatialCoords(spe) <- data.matrix(
              spatialLocs[, seq_len(2)]
            )
        } else {
            if (verbose) {
                message("No spatial locations found in the input Giotto object")
            }
        }

        # DimReductions
        giottoReductions <- list_dim_reductions(
            gobject = giottoObj,
            spat_unit = spatialUnits[su]
        )
        if (!is.null(giottoReductions)) {
            if (verbose) {
                message(
                    "Copying reduced dimensions for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            for (i in seq(nrow(giottoReductions))) {
                SingleCellExperiment::reducedDim(
                    spe,
                    giottoReductions[i]$name
                ) <- get_dimReduction(
                    gobject = giottoObj,
                    reduction = "cells",
                    spat_unit = spatialUnits[su],
                    feat_type = giottoReductions[i]$feat_type,
                    reduction_method = giottoReductions[i]$dim_type,
                    name = giottoReductions[i]$name,
                    output = "data.table"
                )
            }
        } else {
            if (verbose) {
                message("No reduced dimensions found in the input Giotto
                        object")
            }
        }


        # NN Graph
        giottoNearestNetworks <- list_nearest_networks(
            gobject = giottoObj,
            spat_unit = spatialUnits[su]
        )
        if (!is.null(giottoNearestNetworks)) {
            if (verbose) {
                message(
                    "Copying nearest networks for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            for (i in seq(nrow(giottoNearestNetworks))) {
                nn_network <- get_NearestNetwork(
                    gobject = giottoObj,
                    spat_unit = spatialUnits[su],
                    nn_network_to_use = giottoNearestNetworks[i]$type,
                    network_name = giottoNearestNetworks[i]$name,
                    output = "data.table"
                )

                # SPE stores in colpairs, with col indices instead of colnames
                cell1 <- match(nn_network$from, pData$cell_ID)
                cell2 <- match(nn_network$to, pData$cell_ID)

                SingleCellExperiment::colPair(
                    spe,
                    giottoNearestNetworks[i]$name
                ) <- S4Vectors::SelfHits(
                    cell1, cell2,
                    nnode = ncol(spe),
                    nn_network[, seq(-1, -2)] # removing from and to
                )
            }
        } else {
            if (verbose) {
                message("No nearest networks found in the input Giotto object")
            }
        }

        # Spatial Networks
        giottoSpatialNetworks <- list_spatial_networks(
            gobject = giottoObj,
            spat_unit = spatialUnits[su]
        )
        if (!is.null(giottoSpatialNetworks)) {
            if (verbose) {
                message(
                    "Copying spatial networks for spatial unit: '",
                    spatialUnits[su], "'"
                )
            }
            for (i in seq(nrow(giottoSpatialNetworks))) {
                sp_network <- get_spatialNetwork(
                    gobject = giottoObj,
                    spat_unit = spatialUnits[su],
                    name = giottoSpatialNetworks[i]$name,
                    output = "networkDT"
                )

                # spe stores in colpairs, with col indices instead of colnames
                cell1 <- match(sp_network$from, pData$cell_ID)
                cell2 <- match(sp_network$to, pData$cell_ID)

                SingleCellExperiment::colPair(
                    spe,
                    giottoSpatialNetworks[i]$name
                ) <- S4Vectors::SelfHits(
                    cell1, cell2,
                    nnode = ncol(spe),
                    sp_network[, seq(-1, -2)] # removing from and to
                )
            }
        } else {
            if (verbose) {
                message("No spatial networks found in the input Giotto object")
            }
        }

        # SpatialImages
        giottoImages <- list_images(gobject = giottoObj)
        if (!is.null(giottoImages)) {
            for (i in seq(nrow(giottoImages))) {
                img <- get_giottoImage(
                    gobject = giottoObj,
                    image_type = giottoImages[i]$img_type,
                    name = giottoImages[i]$name
                )

                if (!is.null(img@file_path)) {
                    if (verbose) message("Copying spatial image: ", img@name)

                    tryCatch(
                        expr = {
                            spe <- SpatialExperiment::addImg(spe,
                                sample_id = spe$sample_id[i],
                                image_id = img@name,
                                imageSource = img@file_path,
                                scaleFactor = mean(img@scale_factor),
                                load = TRUE
                            )
                        },
                        error = function(e) {
                            message(
                                "Error copying spatial image: ",
                                img@name,
                                ". Please check if the image path is
                                    correct and the image exists at that path."
                            )
                        }
                    )
                } else {
                    if (verbose) {
                        message(
                            "\t - Skipping image with NULL file path: ",
                            img@name
                        )
                    }
                }
                S4Vectors::metadata(spe)[[img@name]] <- img
            }
        } else {
            if (verbose) {
                message("No spatial images found in the input Giotto object")
            }
        }

        if (verbose) message("")

        # Add spe for current spatial unit to speList
        speList[[su]] <- spe
    }

    # return list of spe objects
    return(speList)
}


#' Utility function to convert a SpatialExperiment object to a Giotto object
#'
#' @param spe Input SpatialExperiment object to convert to a Giotto object.
#' @param python_path Specify the path to python.
#' @param nn_network Specify the name of the nearest neighbour network(s)
#' in the input SpatialExperiment object. Default \code{NULL} will use
#' all existing networks.
#' @param sp_network Specify the name of the spatial network(s) in the input
#' SpatialExperiment object. Default \code{NULL} will use all existing
#' networks. This can be a vector of multiple network names.
#' @param verbose A boolean value specifying if progress messages should
#' be displayed or not. Default \code{TRUE}.
#' @import data.table
#' @returns Giotto object
#' @examples
#' \dontrun{
#' library(SpatialExperiment)
#' example(read10xVisium, echo = FALSE)
#' spatialExperimentToGiotto(spe)
#' }
#' @export
spatialExperimentToGiotto <- function(
        spe,
        python_path,
        nn_network = NULL,
        sp_network = NULL,
        verbose = TRUE) {
    # Create giotto instructions and set python path
    instrs <- createGiottoInstructions(python_path = python_path)

    # Create Giotto object with first matrix
    exprMats <- SummarizedExperiment::assays(spe)
    exprMatsNames <- SummarizedExperiment::assayNames(spe)
    firstMatrix <- exprMats[[1]]

    # check unique colnames
    if (length(unique(colnames(firstMatrix))) !=
        length(colnames(firstMatrix))) {
        colnames(firstMatrix) <- make.names(colnames(firstMatrix),
            unique = TRUE
        )
    }

    if (verbose) {
        message("Creating Giotto object with ", exprMatsNames[1], " matrix")
    }
    suppressWarnings(suppressMessages(giottoObj <- createGiottoObject(
        expression = firstMatrix,
        instructions = instrs
    )))
    exprMats[[1]] <- NULL
    exprMatsNames <- exprMatsNames[-1]

    # Copying remaining matrices
    if (length(exprMats) > 0) {
        for (i in seq(exprMats)) {
            if (verbose) {
                message("Copying expression matrix: ", exprMatsNames[i])
            }
            exprObj <- create_expr_obj(
                name = exprMatsNames[i],
                exprMat = exprMats[[i]]
            )
            giottoObj <- set_expression_values(
                gobject = giottoObj,
                values = exprObj
            )
        }
    }

    # Phenotype Data
    pData <- SummarizedExperiment::colData(spe)

    # To handle cell_ID duplication issues
    if ("cell_ID" %in% colnames(pData)) {
        pData$`cell_ID` <- NULL
    }
    if (nrow(pData) > 0) {
        if (verbose) message("Copying phenotype data")
        giottoObj <- addCellMetadata(
            gobject = giottoObj,
            new_metadata = as.data.table(pData)
        )
    }

    # Feature Metadata
    fData <- SummarizedExperiment::rowData(spe)
    if (nrow(fData) > 0) {
        if (verbose) message("Copying feature metadata")
        giottoObj <- addFeatMetadata(
            gobject = giottoObj,
            new_metadata = as.data.table(fData)
        )
    }

    # Reduced Dimensions
    redDims <- SingleCellExperiment::reducedDims(spe)
    redDimsNames <- SingleCellExperiment::reducedDimNames(spe)
    if (length(redDims) > 0) {
        for (i in seq(length(redDims))) {
            if (verbose) message("Copying reduced dimensions")
            dimRedObj <- create_dim_obj(
                name = redDimsNames[i],
                coordinates = redDims[[i]],
                reduction_method = redDimsNames[i]
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            giottoObj <- set_dimReduction(
                gobject = giottoObj,
                dimObject = dimRedObj
            )
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        }
    }

    # Spatial Locations
    spatialLocs <- SpatialExperiment::spatialCoords(spe)
    if (ncol(spatialLocs) > 0) {
        if (verbose) message("Copying spatial locations")
        spatialLocsDT <- data.table::data.table(
            sdimx = spatialLocs[, 1],
            sdimy = spatialLocs[, 2],
            cell_ID = colnames(spe)
        )
        spatLocsObj <- create_spat_locs_obj(
            name = "spatLocs",
            coordinates = spatialLocsDT
        )
        giottoObj <- set_spatial_locations(
            gobject = giottoObj,
            spatlocs = spatLocsObj
        )
    }

    # Spatial Images
    spatialImages <- SpatialExperiment::imgData(spe)
    if (nrow(spatialImages) > 0) {
        for (i in seq(nrow(spatialImages))) {
            if (verbose) message("Copying spatial images")
            spImg <- SpatialExperiment::getImg(
                spe,
                spatialImages[i, "sample_id"],
                spatialImages[i, "image_id"]
            )
            mObject <- magick::image_read(grDevices::as.raster(spImg))
            giottoImage <- createGiottoImage(
                gobject = giottoObj,
                mg_object = mObject,
                scale_factor = spatialImages[i, "scaleFactor"]
            )
            giottoObj <- addGiottoImage(
                gobject = giottoObj,
                images = list(giottoImage)
            )
        }
    }

    # Networks
    networks <- SingleCellExperiment::colPairs(spe)
    # Spatial Networks
    if (!is.null(sp_network)) {
        if (all(sp_network %in% names(networks))) {
            for (i in seq(sp_network)) {
                if (verbose) message("Copying spatial networks")
                networkDT <- as.data.table(networks[[sp_network[i]]])
                networkDT$to <- colnames(spe)[networkDT$to]
                networkDT$from <- colnames(spe)[networkDT$from]
                spatNetObj <- create_spat_net_obj(
                    networkDT = networkDT
                )
                giottoObj <- set_spatialNetwork(
                    gobject = giottoObj,
                    spatial_network = spatNetObj,
                    name = sp_network[i]
                )
                networks[[sp_network[i]]] <- NULL
            }
        }
    }

    # Nearest Neighbour Networks
    if (!is.null(nn_network)) {
        if (nn_network %in% names(networks)) {
            for (i in seq(nn_network)) {
                if (verbose) message("Copying nearest neighbour networks")
                nnNetObj <- create_nn_net_obj(
                    name = nn_network[i],
                    igraph = networks[[nn_network[i]]]
                )
                giottoObj <- set_NearestNetwork(
                    gobject = giottoObj,
                    nn_network = nnNetObj
                )
                networks[[nn_network[i]]] <- NULL
            }
        }
    }

    # if not specified, storing remaining as NN
    if (length(networks) > 0) {
        for (i in seq(networks)) {
            if (verbose) message("Copying additional networks")
            nnNetObj <- create_nn_net_obj(
                name = names(networks)[i],
                igraph = networks[[i]]
            )
            giottoObj <- set_NearestNetwork(
                gobject = giottoObj,
                nn_network = nnNetObj
            )
        }
    }

    return(giottoObj)
}





# SpatialExperiment methods ####

# register methods if SpatialExperiment is present
if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
    setMethod(
        SingleCellExperiment::colData,
        signature("giotto"),
        function(x, spat_unit = NULL, feat_type = NULL, ...) {
            pDataDT(
                gobject = x,
                spat_unit = spat_unit,
                feat_type = feat_type, ...
            ) %>%
                S4Vectors::DataFrame()
        }
    )

    setMethod(
        # TODO generate DataFrame of ENSG and gene symbols
        SingleCellExperiment::rowData,
        signature("giotto"),
        function(x, spat_unit = NULL, feat_type = NULL, ...) {
            fDataDT(
                gobject = x,
                spat_unit = spat_unit,
                feat_type = feat_type, ...
            ) %>%
                S4Vectors::DataFrame()
        }
    )
}






# older Giotto versions ####

#' Convert a master Giotto object to suite
#'
#' @param gobject A Giotto object created with master version
#' @param expression_feat available features (e.g. rna, protein, ...)
#'
#' @returns A Giotto object compatible with suite version
#' @export
#'
giottoMasterToSuite <- function(
        gobject,
        expression_feat = "rna") {
    master_object <- gobject

    spatial_locs <- cell_metadata <- feat_metadata <- instructions <- NULL

    if (!is.null(master_object@spatial_locs)) {
        spatial_locs <- master_object@spatial_locs
    }

    if (!is.null(master_object@cell_metadata)) {
        cell_metadata <- master_object@cell_metadata
    }

    if (!is.null(master_object@gene_metadata)) {
        feat_metadata <- master_object@gene_metadata
        colnames(feat_metadata)[1] <- "feat_ID"
    }

    # instructions
    if (!is.null(master_object@instructions)) {
        instructions <- master_object@instructions
    }

    # create Giotto object
    gobject <- createGiottoObject(
        expression = master_object@raw_exprs,
        expression_feat = expression_feat,
        spatial_locs = spatial_locs,
        cell_metadata = cell_metadata,
        feat_metadata = feat_metadata,
        instructions = instructions
    )

    # add normalized expression matrix
    if (!is.null(master_object@norm_expr)) {
        x <- createExprObj(master_object@norm_expr,
            name = "normalized",
            feat_type = expression_feat
        )
        gobject <- setExpression(gobject,
            x = x,
            spat_unit = "cell",
            feat_type = expression_feat,
            name = "normalized"
        )
    }

    # add scaled expression matrix
    if (!is.null(master_object@norm_scaled_expr)) {
        x <- createExprObj(master_object@norm_scaled_expr,
            name = "scaled",
            feat_type = expression_feat
        )
        gobject <- setExpression(gobject,
            x = x,
            spat_unit = "cell",
            feat_type = expression_feat,
            name = "scaled"
        )
    }

    # dimension_reduction
    if (!is.null(master_object@dimension_reduction)) {
        dimension_reduction_master <- master_object@dimension_reduction

        for (i in names(master_object@dimension_reduction$cells)) {
            j <- names(dimension_reduction_master$cells[[i]])
            dimension_reduction <- createDimObj(
                coordinates = dimension_reduction_master$cells[[i]][[
                    j
                ]]$coordinates,
                name = dimension_reduction_master$cells[[i]][[j]]$name,
                feat_type = expression_feat,
                method = dimension_reduction_master$cells[[i]][[
                    j
                ]]$reduction_method,
                misc = dimension_reduction_master$cells[[i]][[j]]$misc
            )
            gobject <- setDimReduction(gobject,
                dimension_reduction,
                spat_unit = "cell",
                feat_type = expression_feat,
                name = i,
                reduction_method = dimension_reduction_master$cells[[i]][[
                    j
                ]]$reduction_method
            )
        }
    }

    # nn_network
    if (!is.null(master_object@nn_network)) {
        nn_network_master <- master_object@nn_network

        for (i in names(nn_network_master)) {
            for (j in names(nn_network_master[[i]])) {
                k <- names(nn_network_master[[i]][[j]])
                nn_network <- createNearestNetObj(
                    name = i,
                    network = nn_network_master[[i]][[j]][[k]],
                    feat_type = expression_feat
                )

                gobject <- setNearestNetwork(gobject,
                    nn_network,
                    spat_unit = "cell",
                    feat_type = expression_feat,
                    nn_type = i,
                    name = j
                )
            }
        }
    }

    # spatial_network
    if (!is.null(master_object@spatial_network)) {
        spatial_network_master <- master_object@spatial_network$spatial_network
        spatial_network <- createSpatNetObj(
            network = spatial_network_master$networkDT,
            name = spatial_network_master$name,
            networkDT_before_filter =
                spatial_network_master$networkDT_before_filter,
            method = spatial_network_master$method,
            parameters = spatial_network_master$parameters,
            outputObj = spatial_network_master$outputObj,
            cellShapeObj = spatial_network_master$cellShapeObj,
            crossSectionObjects = spatial_network_master$crossSectionObjects,
            misc = spatial_network_master$misc
        )
        gobject <- setSpatialNetwork(gobject,
            spatial_network,
            spat_unit = "cell",
            name = spatial_network_master$name
        )
    }

    # spatial_enrichment
    if (!is.null(master_object@spatial_enrichment)) {
        spatial_enrichment_master <- master_object@spatial_enrichment

        for (i in names(spatial_enrichment_master)) {
            spatial_enrichment <- createSpatEnrObj(
                spatial_enrichment_master[[i]],
                name = i,
                feat_type = expression_feat,
                method = i
            )
        }

        gobject <- set_spatial_enrichment(gobject,
            spatial_enrichment,
            spat_unit = "cell",
            feat_type = expression_feat,
            enrichm_name = i
        )
    }

    # spatial_grid
    if (!is.null(master_object@spatial_grid)) {
        spatial_grid_master <- master_object@spatial_grid

        spatial_grid <- new("spatialGridObj",
            name = spatial_grid_master$spatial_grid$name,
            method = spatial_grid_master$spatial_grid$method,
            parameters = spatial_grid_master$spatial_grid$parameters,
            gridDT = spatial_grid_master$spatial_grid$gridDT,
            feat_type = expression_feat,
            misc = spatial_grid_master$spatial_grid$misc
        )

        gobject <- setSpatialGrid(gobject,
            spatial_grid,
            spat_unit = "cell",
            feat_type = expression_feat,
            name = spatial_grid_master$spatial_grid$name
        )
    }

    return(gobject)
}



## SpatialData object to Giotto ####

#' @title Convert SpatialData to Giotto
#' @name spatialdataToGiotto
#' @description Converts a saved SpatialData object into a Giotto object
#'
#' @param spatialdata_path path to SpatialData object
#' @param n_key_added equivalent of "key_added" argument from scanpy.pp.neighbors().
#'                    If multiple spatial networks are in the anndata object, a list of key_added
#'                    terms may be provided.
#'                    If converting an anndata object from giottoToAnnData, a .txt file may be
#'                    provided, which was generated in that function,
#'                          i.e. {spat_unit}_{feat_type}_nn_network_keys_added.txt
#'                    Cannot be "spatial". This becomes the name of the nearest network in the gobject.
#' @param spatial_n_key_added equivalent of "key_added" argument from squidpy.gr.spatial_neighbors.
#'                            If multiple spatial networks are in the anndata object, a list of key_added
#'                            terms may be provided.
#'                            If converting an anndata object from giottoToAnnData, a .txt file may be
#'                            provided, which was generated in that function,
#'                                i.e. {spat_unit}_{feat_type}_spatial_network_keys_added.txt
#'                            Cannot be the same as n_key_added.
#' @param delaunay_spat_net binary parameter for spatial network. If TRUE, the spatial network is a delaunay network.
#' @param spat_unit desired spatial unit for conversion, default NULL
#' @param feat_type desired feature type for conversion, default NULL
#' @param python_path path to python executable within a conda/miniconda environment
#' @param env_name name of environment containing python_path executable
#'
#' @return Giotto object
#' @details Function in beta. Converts a structured SpatialData file into a Giotto object.
#'    The returned Giotto Object will take default insructions with the
#'    exception of the python path, which may be customized.
#'    See \code{\link{changeGiottoInstructions}} to modify instructions after creation.
#' @export

spatialdataToGiotto <- function(
        spatialdata_path = NULL,
        n_key_added = NULL,
        spatial_n_key_added = NULL,
        delaunay_spat_net = TRUE,
        spat_unit = NULL,
        feat_type = NULL,
        python_path = NULL,
        env_name = NULL) {
    
    # File check
    if (is.null(spatialdata_path)) {
        stop("Please provide a path to SpatialData object for conversion.\n")
    }
    if (!file.exists(spatialdata_path)) {
        stop("The provided path to SpatialData object does not exist.\n")
    }

    # Initialize reticulate
    instrs <- createGiottoInstructions(
        python_path = python_path,
        show_plot = FALSE,
        return_plot = FALSE,
        save_plot = TRUE,
        save_dir = NULL,
        plot_format = NULL,
        dpi = NULL,
        units = NULL,
        height = NULL,
        width = NULL,
        is_docker = FALSE,
        plot_count = 0,
        fiji_path = NULL,
        no_python_warn = FALSE
    )

    # Check spatialdata dependencies
    spatialdata_installed <- checkPythonPackage(package_name = "spatialdata", env_to_use = env_name)

    # Import sd2g, a python module for parsing SpatialData
    sd2g_path <- system.file("python", "sd2g.py", package = "GiottoClass")
    reticulate::source_python(sd2g_path)
    sdata <- read_spatialdata_from_path(spatialdata_path)

    # Extract expression matrix
    expr_df <- extract_expression(sdata)
    cID <- extract_cell_IDs(sdata)
    fID <- extract_feat_IDs(sdata)

    # Extract spatial locations
    spatial_df <- extract_spatial(sdata)
    sp <- parse_obsm_for_spat_locs(sdata)

    # Set up metadata
    cm <- extract_cell_metadata(sdata)
    cm <- as.data.table(cm)
    if ("leiden" %in% names(cm)) {
        cm$leiden <- as.numeric(cm$leiden)
    }

    fm <- extract_feat_metadata(sdata)
    fm <- as.data.table(fm)

    # Create baseline Giotto object
    gobject <- createGiottoObject(
        expression = expr_df,
        spatial_locs = spatial_df,
        instructions = instrs
    )

    # Attach hires image
    raster <- terra::rast(extract_image(sdata))
    giotto_image <- createGiottoLargeImage(raster)
    gobject <- addGiottoLargeImage(gobject = gobject, largeImages = c(giotto_image))

    # Attach metadata
    cm <- readCellMetadata(cm)
    gobject <- setCellMetadata(gobject, x = cm)
    fm <- readFeatMetadata(fm)
    gobject <- setFeatureMetadata(gobject, x = fm)

    spat_unit <- activeSpatUnit(gobject)
    feat_type <- activeFeatType(gobject)

    # Add PCA
    p <- extract_pca(sdata)
    if (!is.null(p)) {
        pca <- p$pca
        evs <- p$eigenvalues
        loads <- p$loadings
        # Add pca to giottoObject
        dobj <- createDimObj(
            coordinates = pca,
            name = "pca",
            spat_unit = spat_unit,
            feat_type = feat_type,
            method = "pca",
            reduction = "cells",
            provenance = NULL,
            misc = list(
                eigenvalues = evs,
                loadings = loads
            ),
            my_rownames = colnames(expr_df)
        )
        gobject <- set_dimReduction(gobject = gobject, dimObject = dobj)
    }

    # Add UMAP
    u <- extract_umap(sdata)
    if (!is.null(u)) {
        # Add UMAP to giottoObject
        dobj <- createDimObj(
            coordinates = u,
            name = "umap",
            spat_unit = spat_unit,
            feat_type = feat_type,
            method = "umap",
            reduction = "cells",
            provenance = NULL,
            misc = NULL,
            my_rownames = colnames(expr_df)
        )
        gobject <- set_dimReduction(gobject = gobject, dimObject = dobj)
    }

    # Add tSNE
    t <- extract_tsne(sdata)
    if (!is.null(t)) {
        # Add TSNE to giottoObject
        dobj <- createDimObj(
            coordinates = t,
            name = "tsne",
            spat_unit = spat_unit,
            feat_type = feat_type,
            method = "tsne",
            reduction = "cells",
            provenance = NULL,
            misc = NULL,
            my_rownames = colnames(expr_df)
        )
        gobject <- set_dimReduction(gobject = gobject, dimObject = dobj)
    }

    ## Nearest Network

    weights_sd <- NULL
    num_NN_nets <- length(n_key_added)

    if (is.null(n_key_added) && !is.null(extract_NN_connectivities(sdata, key_added = n_key_added))) {
        num_NN_nets <- 1
    }

    for (i in num_NN_nets) {
        if (inherits(n_key_added, "list")) {
            n_key_added_it <- n_key_added[[i]]
        } else {
            n_key_added_it <- n_key_added
        }

        weights_sd <- extract_NN_connectivities(sdata, key_added = n_key_added_it)
        # adw = methods::as(weights_ad, "TsparseMatrix")
        if (!is.null(weights_sd)) {
            distances_sd <- extract_NN_distances(sdata, key_added = n_key_added_it)

            nn_dt <- align_network_data(distances = weights_sd, weights = distances_sd)

            # pre-allocate DT variables
            from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <- uniq_ID <- NULL
            nn_dt <- data.table::data.table(nn_dt)

            nn_dt[, from_cell_ID := cID[from]]
            nn_dt[, to_cell_ID := cID[to]]
            nn_dt[, uniq_ID := paste0(from, to)]
            nn_dt[order(uniq_ID)]
            nn_dt[, uniq_ID := NULL]
            vert <- unique(x = c(nn_dt$from_cell_ID, nn_dt$to_cell_ID))
            nn_network_igraph <- igraph::graph_from_data_frame(nn_dt[, .(from_cell_ID, to_cell_ID, weight, distance)], directed = TRUE, vertices = vert)

            nn_info <- extract_NN_info(adata = adata, key_added = n_key_added_it)

            net_type <- "kNN" # anndata default
            if (("sNN" %in% n_key_added_it) & !is.null(n_key_added_it)) {
                net_type <- "sNN"
                net_name <- paste0(n_key_added_it, ".", nn_info["method"])
            } else if (!("sNN" %in% n_key_added_it) & !is.null(n_key_added_it)) {
                net_name <- paste0(n_key_added_it, ".", nn_info["method"])
            } else {
                net_name <- paste0(net_type, ".", nn_info["method"])
            }

            netObj <- createNearestNetObj(
                name = net_name,
                network = nn_network_igraph,
                spat_unit = spat_unit,
                feat_type = feat_type
            )

            gobject <- set_NearestNetwork(
                gobject = gobject,
                nn_network = netObj,
                spat_unit = spat_unit,
                feat_type = feat_type,
                nn_network_to_use = net_type,
                network_name = net_name,
                set_defaults = FALSE
            )
        }
    }

    ## Spatial Network

    s_weights_sd <- NULL
    num_SN_nets <- length(spatial_n_key_added)

    # Check for the case where NULL is provided, since the
    # anndata object takes the default value for SN

    if (is.null(spatial_n_key_added) && !is.null(extract_SN_connectivities(sdata, key_added = spatial_n_key_added))) {
        num_SN_nets <- 1
    }

    for (i in 1:num_SN_nets) {
        if (inherits(spatial_n_key_added, "list")) {
            spatial_n_key_added_it <- spatial_n_key_added[[i]]
        } else {
            spatial_n_key_added_it <- spatial_n_key_added
        }

        s_weights_sd <- extract_SN_connectivities(sdata, key_added = spatial_n_key_added_it)
        if (!is.null(s_weights_sd)) {
            s_distances_sd <- extract_SN_distances(sdata, key_added = spatial_n_key_added_it)
            ij_matrix <- methods::as(s_distances_sd, "TsparseMatrix")
            from_idx <- ij_matrix@i + 1 # zero index!!!
            to_idx <- ij_matrix@j + 1 # zero index!!!

            # pre-allocate DT variables
            from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <- uniq_ID <- NULL
            sn_dt <- data.table::data.table(
                from = from_idx,
                to = to_idx,
                weight = s_weights_sd@x,
                distance = s_distances_sd@x
            )

            sn_dt[, from_cell_ID := cID[from]]
            sn_dt[, to_cell_ID := cID[to]]

            sdimx <- "sdimx"
            sdimy <- "sdimy"
            xbegin_name <- paste0(sdimx, "_begin")
            ybegin_name <- paste0(sdimy, "_begin")
            xend_name <- paste0(sdimx, "_end")
            yend_name <- paste0(sdimy, "_end")

            network_DT <- data.table::data.table(
                from = sn_dt$from_cell_ID,
                to = sn_dt$to_cell_ID,
                xbegin_name = sp[sn_dt$from, sdimx],
                ybegin_name = sp[sn_dt$from, sdimy],
                xend_name = sp[sn_dt$to, sdimx],
                yend_name = sp[sn_dt$to, sdimy],
                weight = s_weights_sd@x,
                distance = s_distances_sd@x
            )
            data.table::setnames(network_DT,
                old = c("xbegin_name", "ybegin_name", "xend_name", "yend_name"),
                new = c(xbegin_name, ybegin_name, xend_name, yend_name)
            )
            data.table::setorder(network_DT, from, to)

            dist_mean <- get_distance(network_DT, method = "mean")
            dist_median <- get_distance(network_DT, method = "median")
            cellShapeObj <- list(
                "meanCellDistance" = dist_mean,
                "medianCellDistance" = dist_median
            )

            # TODO filter network?
            # TODO 3D handling?
            if (delaunay_spat_net) {
                spatObj <- create_spat_net_obj(
                    name = "sNN",
                    method = "delaunay",
                    networkDT = network_DT,
                    cellShapeObj = cellShapeObj
                )
            } else {
                spatObj <- create_spat_net_obj(
                    name = "sNN",
                    method = "non-delaunay",
                    networkDT = network_DT,
                    cellShapeObj = cellShapeObj
                )
            }

            gobject <- set_spatialNetwork(
                gobject = gobject,
                spatial_network = spatObj,
                name = "sNN"
            )
        }
    }
    return(gobject)
}



##  Giotto to SpatialData####

#' @title Convert Giotto to SpatialData
#' @name giottoToSpatialData
#' @description Converts a Giotto object to SpatialData object
#'
#' @param gobject giotto object to be converted
#' @param spat_unit spatial unit which will be used in conversion
#' @param feat_type feature type which will be used in conversion
#' @param spot_radius radius of the spots
#' @param python_path path to python executable within a conda/miniconda environment
#' @param env_name name of environment containing python_path executable
#' @param save_directory directory in which the SpatialData object will be saved
#'
#' @return SpatialData object saved on disk.
#' @details Function in beta. Converts and saves a Giotto object in SpatialData format on disk.
#' @export

giottoToSpatialData <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        spot_radius = NULL,
        python_path = NULL,
        env_name = NULL,
        save_directory = NULL) {
        
    # Initialize reticulate
    instrs <- createGiottoInstructions(python_path = python_path)

    # Check spatialdata dependencies
    spatialdata_installed <- checkPythonPackage(package_name = "spatialdata", env_to_use = env_name)

    # Import sd2g, a python module for parsing SpatialData
    g2sd_path <- system.file("python", "g2sd.py", package = "GiottoClass")
    reticulate::source_python(g2sd_path)

    # Get metadata
    spat_unit <- activeSpatUnit(gobject)
    feat_type <- activeFeatType(gobject)

    # Create a temporary folder to hold anndata
    temp <- "temp_anndata/"

    # First, convert Giotto object to AnnData using an existing function
    giottoToAnnData(
        gobject = gobject,
        spat_unit = spat_unit,
        feat_type = feat_type,
        python_path = python_path,
        env_name = env_name,
        save_directory = temp
    )

    # Extract GiottoImage
    gimg <- getGiottoImage(gobject, image_type = "largeImage")

    # Temporarily save the image to disk
    writeGiottoLargeImage(
        giottoLargeImage = gimg,
        gobject = gobject,
        largeImage_name = "largeImage",
        filename = "temp_image.png",
        dataType = NULL,
        max_intensity = NULL,
        overwrite = TRUE,
        verbose = TRUE
    )
    
    spat_locs <- getSpatialLocations(gobject, output="data.table")

    # Create SpatialData object
    createSpatialData(temp, spat_locs, spot_radius, save_directory)

    # Delete temporary files and folders
    unlink("temp_image.png")
    unlink("temp_image.png.aux.xml")
    unlink(temp, recursive = TRUE)

    # Successful Conversion
    cat("Giotto object has been converted and saved to SpatialData object at: ", save_directory, "\n")

}
