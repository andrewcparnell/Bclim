plot.slice_clouds = function(x,slice = 1, dims=1:2, n=50,...) {

if(length(dims)!=2) stop('Only plots of 2 climate dimensions are currently supported')

  # Get current slice and climate dimensions
  curr_slice = x$slice_clouds[,slice,dims]

  # Create density and standardise
  dens = MASS::kde2d(curr_slice[,1],curr_slice[,2], n=n)
  dens$z = dens$z/sum(dens$z)
  graphics::filled.contour(dens,...)

}

