% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interoperability.R
\name{seuratToGiottoV4}
\alias{seuratToGiottoV4}
\title{Convert a Seurat V4 object to a Giotto object}
\usage{
seuratToGiottoV4(
  sobject,
  spatial_assay = "Spatial",
  dim_reduction = c("pca", "umap"),
  subcellular_assay = "Vizgen",
  sp_network = NULL,
  nn_network = NULL,
  verbose = TRUE
)
}
\arguments{
\item{sobject}{Input Seurat object to convert to Giotto object}

\item{spatial_assay}{Specify name of the spatial assay slot in Seurat.
Default is \code{"Spatial"}.}

\item{dim_reduction}{Specify which dimensional reduction computations to
fetch from input Seurat object. Default is \code{"c('pca', 'umap')"}.}

\item{subcellular_assay}{Specify name of the subcellular assay in input}

\item{sp_network}{sp_network}

\item{nn_network}{nn_network}

\item{verbose}{logical. Default to TRUE
object. Default is \code{"Vizgen"}.}
}
\value{
A Giotto object converted from Seurat object with all computations
stored in it.
}
\description{
Convert a Seurat V4 object to a Giotto object
}
\examples{
m_expression <- Matrix::Matrix(rnorm(100), nrow = 10, sparse = TRUE)
s <- Seurat::CreateSeuratObject(counts = m_expression)

seuratToGiottoV5(s, spatial_assay = "RNA")
}
\keyword{interoperability}
\keyword{seurat}
