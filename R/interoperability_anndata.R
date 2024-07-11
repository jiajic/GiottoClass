
## anndata object ####


#' @title Check Scanpy Installation
#' @name check_py_for_scanpy
#' @import reticulate
#' @returns character
#' @description checks current python environment for scanpy 1.9.0
#' @keywords internal
check_py_for_scanpy <- function() {
    # test if scanpy is found
    module_test <- reticulate::py_module_available("scanpy")
    py_path <- reticulate::py_config()$python
    genv_in_use <- grepl(pattern = "giotto_env", x = py_path)
    
    if (module_test == FALSE && !genv_in_use) {
        stop(wrap_txt("scanpy python module is not installed:
            install in the environment or python path with:

            'pip install scanpy==1.9.0'

            Alternatively, install in the active python
            environment with reticulate:

            reticulate::py_install(packages = 'scanpy==1.9.0',
                                   pip = TRUE)
            \n
            ", errWidth = TRUE))
    } else if (module_test == FALSE && genv_in_use) {
        cat("Python module scanpy is required for conversion.
          Installing scanpy now in the Giotto Miniconda Environment.\n")
        
        conda_path <- reticulate::miniconda_path()
        py_ver <- reticulate::py_config()$version_string
        py_ver <- strsplit(py_ver, "|", fixed = TRUE)[[1]][1]
        py_ver <- gsub(" ", "", py_ver, fixed = TRUE)
        conda_full_path <- paste0(conda_path, "/", "bin/conda")
        full_envname <- paste0(conda_path, "/envs/giotto_env")
        
        reticulate::py_install(
            packages = "scanpy==1.9.0",
            envname = full_envname,
            method = "conda",
            conda = conda_full_path,
            pip = TRUE,
            python_version = py_ver
        )
    } else {
        cat("Required Python module scanpy has been previously installed.
            Proceeding with conversion.\n")
    }
}


#' @title Convert anndata to Giotto
#' @name anndataToGiotto
#' @description Converts a spatial anndata (e.g. scanpy) .h5ad file into a
#' Giotto object
#'
#' @param anndata_path path to the .h5ad file
#' @param n_key_added equivalent of "key_added" argument from scanpy.pp.
#' neighbors(). If multiple spatial networks are in the anndata object, a list
#' of key_added terms may be provided. If converting an anndata object from
#' giottoToAnnData, a .txt file may be provided, which was generated in that
#' function, i.e. \{spat_unit\}_\{feat_type\}_nn_network_keys_added.txt. Cannot
#' be "spatial". This becomes the name of the nearest network in the gobject.
#' @param spatial_n_key_added equivalent of "key_added" argument from
#' squidpy.gr.spatial_neighbors. If multiple spatial networks are in the
#' anndata object, a list of key_added terms may be provided. If converting an
#' anndata object from giottoToAnnData, a .txt file may be provided, which was
#' generated in that function,
#' i.e. {spat_unit}_{feat_type}_spatial_network_keys_added.txt
#' Cannot be the same as n_key_added.
#' @param delaunay_spat_net binary parameter for spatial network. If TRUE, the
#' spatial network is a deluanay network.
#' @param spat_unit desired spatial unit to use for conversion, default NULL
#' @param feat_type desired feature type to use for conversion, default NULL
#' @param h5_file name to create and on-disk HDF5 file
#' @param python_path path to python executable within a conda/miniconda
#' environment
#' @param env_name name of environment containing python_path executable
#'
#' @details Function in beta. Converts a .h5ad file into a Giotto object.
#' The returned Giotto Object will take default insructions with the
#' exception of the python path, which may be customized.
#' See \code{\link{changeGiottoInstructions}} to modify instructions after
#' creation.
#' @returns Giotto object
#' @export
anndataToGiotto <- function(
        anndata_path = NULL,
        n_key_added = NULL,
        spatial_n_key_added = NULL,
        delaunay_spat_net = TRUE,
        spat_unit = NULL,
        feat_type = NULL,
        h5_file = NULL,
        python_path = NULL,
        env_name = "giotto_env") {
    # Preliminary file checks and guard clauses
    if (is.null(anndata_path)) {
        stop("Please provide a path to an AnnData .h5ad file for conversion.\n")
    }
    
    if (!file.exists(anndata_path)) {
        stop("The provided path to the AnnData .h5ad file does not exist.\n")
    }
    if (!is.null(n_key_added) && !is.null(spatial_n_key_added)) {
        for (n in n_key_added) {
            for (s in spatial_n_key_added) {
                if (n == s) {
                    stop("Arguments n_key_added and spatial_n_key_added may not
                        take the same value.")
                }
            }
        }
    }
    
    # Required step to properly initialize reticulate
    set_giotto_python_path(python_path = python_path)
    
    checkPythonPackage("anndata", env_to_use = env_name)
    checkPythonPackage("scanpy", env_to_use = env_name)
    # should trigger a stop() downstream if not installed
    
    # Import ad2g, a python module for parsing anndata
    ad2g_path <- system.file("python", "ad2g.py", package = "GiottoClass")
    reticulate::source_python(ad2g_path)
    adata <- read_anndata_from_path(anndata_path)
    
    ### Set up expression matrix
    X <- extract_expression(adata)
    cID <- extract_cell_IDs(adata)
    fID <- extract_feat_IDs(adata)
    X@Dimnames[[1]] <- fID
    X@Dimnames[[2]] <- cID
    # Expression matrix X ready
    
    ### Set up spatial info
    sp <- parse_obsm_for_spat_locs(adata)
    # Spatial locations sp ready
    
    ### Set up metadata
    cmeta <- extract_cell_metadata(adata)
    cmeta <- as.data.table(cmeta)
    if ("leiden" %in% names(cmeta)) {
        cmeta$leiden <- as.numeric(cmeta$leiden)
    }
    
    fm <- extract_feat_metadata(adata)
    fm <- as.data.table(fm)
    # Metadata ready
    
    ### Create Minimal giottoObject
    gobject <- createGiottoObject(
        expression = X,
        spatial_locs = sp,
        h5_file = h5_file
    )
    
    ### Add metadata
    cmeta <- readCellMetadata(cmeta)
    gobject <- setCellMetadata(gobject,
                               x = cmeta
    )
    fm <- readFeatMetadata(fm)
    gobject <- setFeatureMetadata(gobject,
                                  x = fm
    )
    
    spat_unit <- set_default_spat_unit(gobject,
                                       spat_unit = spat_unit
    )
    feat_type <- set_default_feat_type(gobject,
                                       feat_type = feat_type,
                                       spat_unit = spat_unit
    )
    
    ### Set up PCA
    p <- extract_pca(adata)
    if (!is.null(p)) {
        pca <- p$pca
        evs <- p$eigenvalues
        loads <- p$loadings
        # Add PCA to giottoObject
        dobj <- create_dim_obj(
            name = "pca.ad",
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = NULL,
            reduction = "cells",
            reduction_method = "pca",
            coordinates = pca,
            misc = list(
                eigenvalues = evs,
                loadings = loads
            ),
            my_rownames = colnames(X)
        )
        
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dobj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }
    
    ### Set up UMAP
    u <- extract_umap(adata)
    if (!is.null(u)) {
        # Add UMAP to giottoObject
        dobj <- create_dim_obj(
            name = "umap.ad",
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = NULL,
            reduction = "cells",
            reduction_method = "umap",
            coordinates = u,
            misc = NULL,
            my_rownames = colnames(X)
        )
        
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dobj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }
    ### Set up TSNE
    t <- extract_tsne(adata)
    if (!is.null(t)) {
        # Add TSNE to giottoObject
        dobj <- create_dim_obj(
            name = "tsne.ad",
            spat_unit = spat_unit,
            feat_type = feat_type,
            provenance = NULL,
            reduction = "cells",
            reduction_method = "tsne",
            coordinates = t,
            misc = NULL,
            my_rownames = colnames(X)
        )
        
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
        gobject <- set_dimReduction(gobject = gobject, dimObject = dobj)
        ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
    }
    
    ### NN Network
    
    # Need to create nnNetObj or igraph object to use with setter for NN
    
    weights_ad <- NULL
    num_NN_nets <- length(n_key_added)
    
    if (is.null(n_key_added) &&
        !is.null(extract_NN_connectivities(adata, key_added = n_key_added))) {
        num_NN_nets <- 1
    }
    
    for (i in num_NN_nets) {
        if (inherits(n_key_added, "list")) {
            n_key_added_it <- n_key_added[[i]]
        } else {
            n_key_added_it <- n_key_added
        }
        
        weights_ad <- extract_NN_connectivities(adata,
                                                key_added = n_key_added_it
        )
        # adw = methods::as(weights_ad, "TsparseMatrix")
        if (!is.null(weights_ad)) {
            distances_ad <- extract_NN_distances(adata,
                                                 key_added = n_key_added_it
            )
            
            nn_dt <- align_network_data(
                distances = weights_ad,
                weights = distances_ad
            )
            
            # pre-allocate DT variables
            from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <-
                uniq_ID <- NULL
            nn_dt <- data.table::data.table(nn_dt)
            
            nn_dt[, from_cell_ID := cID[from]]
            nn_dt[, to_cell_ID := cID[to]]
            nn_dt[, uniq_ID := paste0(from, to)]
            nn_dt[order(uniq_ID)]
            nn_dt[, uniq_ID := NULL]
            vert <- unique(x = c(nn_dt$from_cell_ID, nn_dt$to_cell_ID))
            nn_network_igraph <- igraph::graph_from_data_frame(
                nn_dt[, .(from_cell_ID, to_cell_ID, weight, distance)],
                directed = TRUE, vertices = vert
            )
            
            nn_info <- extract_NN_info(
                adata = adata,
                key_added = n_key_added_it
            )
            
            net_type <- "kNN" # anndata default
            if (("sNN" %in% n_key_added_it) & !is.null(n_key_added_it)) {
                net_type <- "sNN"
                net_name <- paste0(n_key_added_it, ".", nn_info["method"])
            } else if (!("sNN" %in% n_key_added_it) &
                       !is.null(n_key_added_it)) {
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
    s_weights_ad <- NULL
    num_SN_nets <- length(spatial_n_key_added)
    
    # Check for the case where NULL is provided, since the
    # anndata object takes the default value for SN
    
    if (is.null(spatial_n_key_added) &&
        !is.null(extract_SN_connectivities(adata,
                                           key_added = spatial_n_key_added
        ))) {
        num_SN_nets <- 1
    }
    
    for (i in seq_len(num_SN_nets)) {
        if (inherits(spatial_n_key_added, "list")) {
            spatial_n_key_added_it <- spatial_n_key_added[[i]]
        } else {
            spatial_n_key_added_it <- spatial_n_key_added
        }
        
        s_weights_ad <- extract_SN_connectivities(
            adata,
            key_added = spatial_n_key_added_it
        )
        if (!is.null(s_weights_ad)) {
            s_distances_ad <- extract_SN_distances(
                adata,
                key_added = spatial_n_key_added_it
            )
            ij_matrix <- methods::as(s_distances_ad, "TsparseMatrix")
            from_idx <- ij_matrix@i + 1 # zero index!!!
            to_idx <- ij_matrix@j + 1 # zero index!!!
            
            # pre-allocate DT variables
            from <- to <- weight <- distance <- from_cell_ID <- to_cell_ID <-
                uniq_ID <- NULL
            sn_dt <- data.table::data.table(
                from = from_idx,
                to = to_idx,
                weight = s_weights_ad@x,
                distance = s_distances_ad@x
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
                weight = s_weights_ad@x,
                distance = s_distances_ad@x
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
                    name = "Spat_Net_from_AnnData",
                    method = "delaunay",
                    networkDT = network_DT,
                    cellShapeObj = cellShapeObj
                )
            } else {
                spatObj <- create_spat_net_obj(
                    name = "Spat_Net_from_AnnData",
                    method = "non-delaunay",
                    networkDT = network_DT,
                    cellShapeObj = cellShapeObj
                )
            }
            
            gobject <- set_spatialNetwork(
                gobject = gobject,
                spatial_network = spatObj,
                name = "Spat_Net_from_AnnData"
            )
        }
    }
    
    ### Layers
    lay_names <- extract_layer_names(adata)
    if (!is.null(lay_names)) {
        for (l_n in lay_names) {
            lay <- extract_layered_data(adata, layer_name = l_n)
            if ("data.frame" %in% class(lay)) {
                names(lay) <- fID
                row.names(lay) <- cID
            } else {
                lay@Dimnames[[1]] <- fID
                lay@Dimnames[[2]] <- cID
            }
            layExprObj <- createExprObj(lay, name = l_n)
            gobject <- set_expression_values(
                gobject = gobject,
                spat_unit = spat_unit,
                feat_type = feat_type,
                name = l_n,
                values = layExprObj
            )
        }
    }
    
    gobject <- update_giotto_params(
        gobject = gobject,
        description = "_AnnData_Conversion"
    )
    
    wrap_msg("\nAnnData object successfully converted to Giotto.\n")
    return(gobject)
}


#' @title Convert Giotto to anndata
#' @name giottoToAnnData
#' @description Converts a Giotto object to a spatial anndata (e.g. scanpy)
#' .h5ad file
#' @param gobject giotto object to be converted
#' @param spat_unit spatial unit which will be used in conversion.
#' @param feat_type feature type which will be used in conversion.
#' @param python_path path to python executable within a conda/miniconda
#' environment
#' @param env_name name of environment containing python_path executable
#' @param save_directory directory in which the file will be saved.
#' @details Function in beta. Converts a Giotto object into .h5ad file(s).
#'
#' If there are multiple spatial units and/or feature types, but only
#' one spatial unit and/or feature type is specified, then only the
#' specified spatial unit and/or feature type will be used. If NULL,
#' by default, all spatial units will be used in conversion.
#'
#' If multiple spatial units or feature types are specified, multiple
#' AnnData object will be created and returned.
#'
#' This function will create .txt files which will record any `key_added`
#' parameters for networks. They are named after the corresponding spatial unit
#' and feature type pair.
#'
#' The save_directory will be created if it does not already exist.
#' The default save_directory is the working directory.
#' @returns vector containing .h5ad file path(s)
#' @export
giottoToAnnData <- function(
        gobject = NULL,
        spat_unit = NULL,
        feat_type = NULL,
        python_path = NULL,
        env_name = "giotto_env",
        save_directory = NULL) {
    # Check gobject
    invalid_obj <- !("giotto" %in% class(gobject))
    if (is.null(gobject) || invalid_obj) {
        stop(wrap_msg("Please provide a valid Giotto Object for conversion."))
    }
    
    checkPythonPackage("anndata", env_to_use = env_name)
    checkPythonPackage("scanpy", env_to_use = env_name)
    
    # Python module import
    g2ad_path <- system.file("python", "g2ad.py", package = "GiottoClass")
    reticulate::source_python(g2ad_path)
    if (!is.null(save_directory)) dir_guard(save_directory)
    
    # Check directory, make it if it doesn't exist
    if (is.null(save_directory)) {
        save_directory <- paste0(getwd(), "/")
    } else if (!dir.exists(save_directory)) {
        warning(wrap_msg("Provided save directory not found. Creating save
                        directory at location:"))
        cat(save_directory)
        dir.create(save_directory, recursive = TRUE)
        if (dir.exists(save_directory)) {
            cat("Created directory", save_directory)
        } else {
            stop(wrap_msg("Unable to create directory. Please change the
                        provided path and try again."))
        }
    } else {
        wrap_msg(
            "Directory", save_directory,
            "found. The converted Giotto object will be saved here as an
                .h5ad file."
        )
    }
    
    # Expresion
    expr_dt <- list_expression(gobject)
    # ID spat_unit and feat_type if not already provided.
    if (is.null(spat_unit) && is.null(feat_type)) {
        spat_unit <- unique(expr_dt$spat_unit)
        feat_type <- unique(expr_dt$feat_type)
    } else if (is.null(spat_unit) && !is.null(feat_type)) {
        spat_unit <- unique(expr_dt$spat_unit)
    } else if (!is.null(spat_unit) && is.null(feat_type)) {
        feat_type <- unique(expr_dt$feat_type)
    }
    
    for (su in spat_unit) {
        wrap_msg("Spatial unit(s)", su, "will be used in conversion.")
    }
    for (ft in feat_type) {
        wrap_msg("Feature type(s)", ft, "will be used in conversion.")
    }
    
    # Iterate through spat_unit and feat_type to pull out expression data.
    # By default, the raw expression in the slot of the first spatial unit
    # and first feature type (if multiple are present) will be transferred to
    # the AnnData.anndata.X slot
    # Any other expression data will be inserted into AnnData.anndata.layers
    # By default, layer names are formed by
    # "'spatial_unit'_'feature_type'_'value'"
    
    adata <- NULL # scope
    su_ft_length <- 0
    
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            su_ft_length <- su_ft_length + 1
        }
    }
    
    
    adata_list <- lapply(seq_len(su_ft_length), function(i) adata)
    adata_pos <- 1
    
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            expr_names <- list_expression_names(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft
            )
            
            for (en in expr_names) {
                if (en == "raw") {
                    raw_x <- get_expression_values(
                        gobject = gobject,
                        values = en,
                        spat_unit = su,
                        feat_type = ft,
                        output = "matrix"
                    )
                    
                    adata <- ad_obj(x = raw_x)
                } else {
                    ad_layer_name <- paste0(su, "_", ft, "_", en)
                    
                    x <- get_expression_values(
                        gobject = gobject,
                        values = en,
                        spat_unit = su,
                        feat_type = ft,
                        output = "matrix"
                    )
                    
                    if ("dgeMatrix" %in% class(x)) x <- methods::as(x, "array")
                    
                    adata <- set_adg_layer_data(
                        adata = adata,
                        lay = x,
                        lay_name = ad_layer_name
                    )
                }
            }
            adata_list[[adata_pos]] <- adata
            adata_pos <- adata_pos + 1
            adata <- NULL
        }
    }
    # Reset indexing variable
    adata_pos <- 1
    
    # Spatial Locations
    for (su in spat_unit) {
        for (ft_ in names(gobject@expression[[su]])) {
            sl <- get_spatial_locations(
                gobject = gobject,
                output = "data.table",
                spat_unit = su
            )
            n_col_sl <- dim(sl)[2]
            
            # preallocate data.table params
            sdimx <- sdimy <- sdimz <- NULL
            
            if (n_col_sl == 3) {
                sl <- sl[, .(sdimx, sdimy)]
            } else {
                sl <- sl[, .(sdimx, sdimy, sdimz)]
            }
            adata <- adata_list[[adata_pos]]
            adata <- set_adg_spat_locs(
                adata = adata,
                spat_locs = sl
            )
            adata_pos <- adata_pos + 1
        }
    }
    # Reset indexing variable
    adata_pos <- 1
    
    # Spatial Info
    
    # Cell Metadata
    # Feat Metadata
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            cmeta <- get_cell_metadata(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft,
                output = "data.table",
                set_defaults = FALSE
            )
            
            fm <- get_feature_metadata(
                gobject = gobject,
                spat_unit = su,
                feat_type = ft,
                output = "data.table",
                set_defaults = FALSE
            )
            
            adata_list[[adata_pos]] <- set_adg_metadata(
                adata = adata_list[[adata_pos]],
                cell_meta = cmeta,
                feat_meta = fm
            )
            
            adata_pos <- adata_pos + 1
        }
    }
    # Reset indexing variable
    adata_pos <- 1
    
    # Dimension Reductions
    
    # error hanldling wrapper to get_dimReduction
    try_get_dimReduction <- function(
        gobject,
        spat_unit,
        feat_type,
        reduction,
        reduction_method,
        name,
        output,
        set_defaults) {
        tryCatch(
            {
                dim_red <- get_dimReduction(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    reduction = reduction,
                    reduction_method = reduction_method,
                    name = name,
                    output = output,
                    set_defaults = set_defaults
                )
                return(dim_red)
            },
            error = function(e) {
                return(NULL)
            }
        )
    }
    
    ## PCA
    
    # pca on feats not supported by anndata because of dimensionality
    # agreement reqs
    reduction_options <- names(gobject@dimension_reduction)
    dim_red <- NULL
    
    for (ro in reduction_options) {
        if (ro != "cells") {
            warning("AnnData does not support storing PCA by features.
                    Skipping PCA data conversion.")
            break
        }
        for (su in spat_unit) {
            for (ft in names(gobject@expression[[su]])) {
                name <- "pca"
                if (ft != "rna") name <- paste0(ft, ".pca")
                dim_red <- try_get_dimReduction(
                    gobject = gobject,
                    spat_unit = su,
                    feat_type = ft,
                    reduction = ro,
                    reduction_method = "pca",
                    name = name,
                    output = "dimObj",
                    set_defaults = FALSE
                )
                if (is.null(dim_red)) {
                    adata_pos <- adata_pos + 1
                    next
                }
                pca_coord <- dim_red[]
                pca_loadings <- data.table(dim_red@misc$loadings)
                feats_used <- dimnames(dim_red@misc$loadings)[[1]]
                evs <- dim_red@misc$eigenvalues
                
                adata_list[[adata_pos]] <- set_adg_pca(
                    adata = adata_list[[adata_pos]],
                    pca_coord = pca_coord,
                    loadings = pca_loadings,
                    eigenv = evs,
                    feats_used = feats_used
                )
                adata_pos <- adata_pos + 1
            }
        }
        adata_pos <- 1
    }
    
    # Reset indexing variable
    adata_pos <- 1
    
    ## UMAP
    for (ro in reduction_options) {
        for (su in spat_unit) {
            for (ft in names(gobject@expression[[su]])) {
                name <- "umap"
                if (ft != "rna") name <- paste0(ft, ".umap")
                dim_red <- try_get_dimReduction(
                    gobject = gobject,
                    spat_unit = su,
                    feat_type = ft,
                    reduction = ro,
                    reduction_method = "umap",
                    name = name,
                    output = "dimObj",
                    set_defaults = FALSE
                )
                
                if (is.null(dim_red)) {
                    adata_pos <- adata_pos + 1
                    next
                }
                umap_data <- dim_red[]
                adata_list[[adata_pos]] <- set_adg_umap(
                    adata = adata_list[[adata_pos]],
                    umap_data = umap_data
                )
                adata_pos <- adata_pos + 1
            }
        }
        # Reset indexing variable
        adata_pos <- 1
    }
    
    # Reset indexing variable
    adata_pos <- 1
    
    ## T-SNE
    for (ro in reduction_options) {
        for (su in spat_unit) {
            for (ft in names(gobject@expression[[su]])) {
                name <- "tsne"
                if (ft != "rna") name <- paste0(ft, ".tsne")
                dim_red <- try_get_dimReduction(
                    gobject = gobject,
                    spat_unit = su,
                    feat_type = ft,
                    reduction = ro,
                    reduction_method = "tsne",
                    name = name,
                    output = "dimObj",
                    set_defaults = FALSE
                )
                
                if (is.null(dim_red)) {
                    adata_pos <- adata_pos + 1
                    next
                }
                tsne_data <- dim_red[]
                adata_list[[adata_pos]] <- set_adg_tsne(
                    adata = adata_list[[adata_pos]],
                    tsne_data = tsne_data
                )
                adata_pos <- adata_pos + 1
            }
        }
        # Reset indexing variable
        adata_pos <- 1
    }
    
    # Reset indexing variable
    adata_pos <- 1
    
    # Nearest Neighbor Network
    
    # error hanldling wrapper to get_NearestNetwork
    try_get_NN <- function(
        gobject,
        spat_unit,
        feat_type,
        nn_network_to_use,
        network_name,
        output,
        set_defaults) {
        tryCatch(
            {
                nearest_net <- get_NearestNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    feat_type = feat_type,
                    nn_network_to_use = nn_network_to_use,
                    network_name = network_name,
                    output = output,
                    set_defaults = set_defaults
                )
                return(nearest_net)
            },
            error = function(e) {
                return(NULL)
            }
        )
    }
    
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            nn_network_to_use <- c("sNN", "kNN")
            for (nn_net_tu in nn_network_to_use) {
                network_name <- list_nearest_networks_names(
                    gobject = gobject,
                    spat_unit = su,
                    feat_type = ft,
                    nn_type = nn_net_tu
                )
                if (is.null(network_name)) {
                    next
                }
                for (n_name in network_name) {
                    gob_NN <- try_get_NN(
                        gobject = gobject,
                        spat_unit = su,
                        feat_type = ft,
                        nn_network_to_use = nn_net_tu,
                        network_name = n_name,
                        output = "nnNetObj",
                        set_defaults = FALSE
                    )
                    
                    pidx <- grep("nn_network", names(gobject@parameters))
                    for (p in pidx) {
                        if (gobject@parameters[[p]]["type"] == nn_net_tu) {
                            kval <- gobject@parameters[[p]]["k"]
                            dim_red_used <- gobject@parameters[[p]][
                                "dim_red_to_use"
                            ]
                        }
                    }
                    
                    df_gob_NN <- igraph::as_data_frame(gob_NN[])
                    
                    adata_list[[adata_pos]] <- set_adg_nn(
                        adata = adata_list[[adata_pos]],
                        df_NN = df_gob_NN,
                        net_name = n_name,
                        n_neighbors = kval,
                        dim_red_used = dim_red_used
                    )
                }
                
                fname_nn <- paste0(su, "_", ft, "_nn_network_keys_added.txt")
                network_name <- network_name[!grepl("kNN.", network_name)]
                append_n <- FALSE
                if (length(network_name) != 0) {
                    if (nn_net_tu == "kNN") append_n <- TRUE
                    write(network_name, fname_nn, append = append_n)
                }
            }
            adata_pos <- adata_pos + 1
        }
    }
    
    # Reset indexing variable
    adata_pos <- 1
    
    try_get_SN <- function(
        gobject,
        spat_unit,
        name,
        output,
        set_defaults,
        verbose) {
        tryCatch(
            {
                spatial_net <- get_spatialNetwork(
                    gobject = gobject,
                    spat_unit = spat_unit,
                    name = name,
                    output = output,
                    set_defaults = set_defaults,
                    verbose = verbose
                )
                return(spatial_net)
            },
            error = function(e) {
                return(NULL)
            }
        )
    }
    
    
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            # Spatial networks do not have a feature type slot.
            # Iterate through anyways to properly assign to anndata objects
            network_name <- list_spatial_networks_names(
                gobject = gobject,
                spat_unit = su
            )
            if (is.null(network_name)) {
                next
            }
            for (sn_name in network_name) {
                gob_SN <- try_get_SN(
                    gobject = gobject,
                    spat_unit = su,
                    name = sn_name,
                    output = "networkDT",
                    set_defaults = FALSE,
                    verbose = TRUE
                )
                
                pidx <- grep("spatial_network", names(gobject@parameters))
                for (p in pidx) {
                    if (gobject@parameters[[p]]["name of spatial network"] == sn_name) {
                        current_param <- gobject@parameters[[p]]
                        kval <- current_param["k neighbours"]
                        maxdist <- current_param["maximum distance threshold"]
                        dimused <- current_param["dimensions used"]
                    }
                }
                
                adata_list[[adata_pos]] <- set_adg_sn(
                    adata = adata_list[[adata_pos]],
                    df_SN = gob_SN,
                    net_name = sn_name,
                    n_neighbors = kval,
                    max_distance = maxdist,
                    dim_used = dimused
                )
            }
            
            fname_sn <- paste0(su, "_", ft, "_spatial_network_keys_added.txt")
            if (length(network_name) != 0) write(network_name, fname_sn)
        }
        
        adata_pos <- adata_pos + 1
    }
    
    # Reset indexing variable
    adata_pos <- 1
    
    # Write AnnData object to .h5ad file
    # Verify it exists, and return upon success
    fname_list <- lapply(seq_len(su_ft_length), function(i) NULL)
    wrap_msg("\n")
    for (su in spat_unit) {
        for (ft in names(gobject@expression[[su]])) {
            adata <- adata_list[[adata_pos]]
            path_adata <- write_ad_h5ad(
                adata = adata,
                save_directory = save_directory,
                spat_unit = su,
                feat_type = ft
            )
            if (!is.null(path_adata)) {
                wrap_msg(
                    "Spatial unit", su, "and feature type",
                    ft, "converted to:"
                )
                wrap_msg(path_adata)
                fname_list[[adata_pos]] <- path_adata
                adata_pos <- adata_pos + 1
            } else {
                wrap_msg(
                    "Unable to convert spatial unit feature type pair",
                    su, ft
                )
                stop(wrap_msg("File writing error. Please try again."))
            }
        }
    }
    
    wrap_msg("\nGiotto object successfully converted to .h5ad file(s)\n")
    
    return(fname_list)
}
