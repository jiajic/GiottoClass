#' @include generics.R
NULL

#' @title Row and column names
#' @name row-plus-colnames-generic
#' @aliases colnames rownames
#' @description Retrieve or set the row or column names of an object
#' @param x object
#' @return A character vector of row or col names
#' @examples
#' g <- GiottoData::loadSubObjectMini("exprObj")
#'
#' colnames(g)
NULL

#' @title Dimnames of an object
#' @name dimnames
#' @description
#' Retrieve or set the dimnames of an object
#' @param x object
#' @returns list of 2 character vectors. First is rownames. Second is colnames
#' @examples
#' g <- GiottoData::loadSubObjectMini("exprObj")
#'
#' dimnames(g)
NULL

#' @title Attribute names of an object
#' @name names
#' @description
#' Retrieve or set the names of an object
#' @param x object
#' @returns 
#' For `names`, a character vector
#' @examples
#' g <- GiottoData::loadSubObjectMini("exprObj")
#'
#' names(g)
NULL

#' @rdname row-plus-colnames-generic
#' @export
setMethod("colnames", signature(x = "exprObj"), function(x) colnames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod("colnames", signature(x = "metaData"), function(x) colnames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod("colnames", signature(x = "spatEnrObj"), function(x) colnames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod("colnames", signature(x = "spatLocsObj"), function(x) colnames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod("colnames", signature(x = "dimObj"), function(x) colnames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod("colnames", signature(x = "spatialNetworkObj"), 
          function(x) colnames(x[]))



#' @rdname row-plus-colnames-generic
#' @export
setMethod("rownames", signature(x = "exprObj"), function(x) rownames(x[]))

#' @rdname row-plus-colnames-generic
#' @export
setMethod("rownames", signature(x = "dimObj"), function(x) rownames(x[]))





#' @rdname dimnames
#' @export
setMethod("dimnames", signature(x = "exprObj"), function(x) dimnames(x[]))

#' @rdname dimnames
#' @export
setMethod("dimnames", signature(x = "dimObj"), function(x) dimnames(x[]))




#' @rdname names
#' @export
setMethod("names", signature("dimObj"), function(x) colnames(x[]))

#' @rdname names
#' @export
setMethod("names", signature("terraVectData"), function(x) names(x[]))

#' @rdname names
#' @export
setMethod("names", signature("metaData"), function(x) names(x[]))

#' @rdname names
#' @export
setMethod("names", signature("nnNetObj"), 
          function(x) igraph::edge_attr_names(x[]))

#' @rdname names
#' @export
setMethod("names", signature("spatialNetworkObj"), function(x) names(x[]))

#' @rdname names
#' @export
setMethod("names", signature("spatEnrObj"), function(x) names(x[]))
