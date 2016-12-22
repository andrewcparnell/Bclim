#' Plots of Bclim slice clouds
#'
#' Create bivariate climate plots of individual slices. For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the author's personal webpage (https://maths.ucd.ie/parnell)?
#'
#' @param x The output of a run from \code{\link{slice_clouds}}
#' @param slice The chosen slice to plot
#' @param dims A vector of length 2. dim=1 corresponds to GDD5, dim=2 to MTCO, and dim=3 to AET/PET
#' @param n The resolution of the resulting plot. A higher value of n will yield finer plots but might require some colour adjustment
#' @param ... Other arguments to the plot function such as axis labels, titles, and colours
#'
#' @details This function creates a bivariate density plot of two climate dimensions (two of GDD5, MTCO and AET/PET) using the MASS library function \code{\link{kde2d}}
#'
#' @return Just a plot
#' 
#' @seealso The main Bclim functions are \code{\link{slice_clouds}} and \code{\link{climate_histories}}
#' @export
plot.slice_clouds = function(x,slice = 1, dims=1:2, n=50,...) {

if(length(dims)!=2) stop('Only plots of 2 climate dimensions are currently supported')

  # Get current slice and climate dimensions
  curr_slice = x$slice_clouds[,slice,dims]

  # Create density and standardise
  dens = MASS::kde2d(curr_slice[,1],curr_slice[,2], n=n)
  dens$z = dens$z/sum(dens$z)
  graphics::filled.contour(dens,...)

}

