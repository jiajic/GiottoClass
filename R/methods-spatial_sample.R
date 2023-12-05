

# density mapping
# mapping to how many clusters are expected per bin

#' @export
setMethod(
  'sample_giotto', signature('giotto'),
  function(
    x,
    spat_unit = NULL,
    feat_type = NULL,
    method = 'random',
    subset = FALSE,
    bins,
    bin_shape,
    ...
  )
  {

  }
)


#' @export
setMethod(
  'sample_giotto',
  signature('giottoPolygon'),
  function(
    x,
    strata,
    ...
  ) {
    # pass to SpatVector points method
    sample_giotto(centroids(x))
  }
)


#' @export
setMethod(
  'sample_giotto',
  signature(x = 'spatLocsObj'),
  function(
    x,
    ...
  ) {
    sv <- terra::vect(
      x[][, c('sdimx', 'sdimy')],
      geom = c('sdimx', 'sdimy')
    )

    # pass to SpatVector points method
    sample_giotto(sv)
  }
)


#' @param sample_ratio Proportion of the dataset to sample. Used to calculate
#' `size` param if it is not directly given
#' @param size a non-negative integer giving the number of items to sample.
#' @export
setMethod(
  'sample_giotto',
  signature(x = 'SpatVector'),
  function(
    x,
    sample_ratio = 0.1,
    size = ceiling(nrow(x) * sample_ratio),
    strata,
    density = TRUE,
    ...
  ) {
    # only points will work for this
    checkmate::assert_true(terra::geomtype(x) == 'points')
    size = ceiling(size)


  }
)


gpoints <- GiottoData::loadSubObjectMini('giottoPoints')

r <- terra::rast(nrows = 20, ncol = 20, extent = ext(gpoints))
r <- terra::rasterize(
  gpoint@spatVector,
  r,
  fun = 'count',
  update = TRUE,
  na.rm = TRUE
)









