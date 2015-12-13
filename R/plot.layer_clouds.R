plot.layer_clouds = function(x,layer = 1, dims=1:2, n=50,...) {
  
if(length(dims)!=2) stop('Only plots of 2 climate dimensions are currently supported')  
  
  # Get current layer and climate dimensions
  curr_layer = x$layer_clouds[,layer,dims]

  # Create density and standardise
  dens = MASS::kde2d(curr_layer[,1],curr_layer[,2], n=n)
  dens$z = dens$z/sum(dens$z)
  filled.contour(dens,...)

}
    
    