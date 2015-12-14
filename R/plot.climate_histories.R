plot.climate_histories = function(x,dim=1,layer_clouds=TRUE,chron=NULL,climate_ribbon=TRUE,most_representative=1,conf=c(0.95,0.75,0.5), col_clouds = rgb(0,0,1,0.2), col_ribbon=rgb(1,0,0,0.4),col_representative=rgb(0,1,0),present_left=TRUE,...) {

# First create base plot
x_range = range(x$time_grid)
if(!present_left) xrange = rev(xrange)
y_range = range(x$layers$layer_clouds[,,dim])

plot(1,1,type="n",xlim=x_range,ylim=y_range,xaxt='n',yaxt='n',...)
axis(side=1,at=pretty(x$time_grid,n=10))
axis(side=2,at=pretty(x$layers$layer_clouds[,,dim],n=10))
grid()

# Get some objects to make coding neater
n_layers = x$layers$n_layers
n_samples = x$layers$n_samples

# Second add in layer clouds if required
if(layer_clouds) {
  if(is.null(chron)) stop("A chronology is required for plotting layer clouds")
  for(i in 1:n_layers) {
    # Get current MDPs and a suitable number of chronologies to match
    curr_MDP = x$layers$layer_clouds[,i,dim]
    curr_chron = sample(chron[,i],size=n_samples,replace=TRUE)

    if(var(curr_chron)>0) {
      dens = MASS::kde2d(curr_chron,curr_MDP)

      # Standardise and find 95% limit
      dens$z = dens$z/sum(dens$z)
      limits_diff = (max(dens$z)-min(dens$z))/100

      # Loop through confidence levels to plot
      for(k in 1:length(conf)) {
        limits = min(dens$z)
        prop = 1
        while(prop>conf[k]) {
          limits = limits+limits_diff
          prop = sum(dens$z[dens$z>limits])
        }

        con_lines = contourLines(dens$x,dens$y,dens$z,levels=limits)
        for(j in 1:length(con_lines)) polygon(con_lines[[j]]$x,con_lines[[j]]$y,col=col_clouds,border=col_clouds)
      # End of loop through confidence levels
      }
    # End of if statement for variance of chronology
    } else {
      # If there's no variance so it's a vertical stripe
      for(k in 1:length(conf)) {
        lines(c(curr_chron[1],curr_chron[1]),quantile(curr_MDP,probs=c((1-conf[k])/2,conf[k]+(1-conf[k])/2)),col=col_clouds)
      }
    # End of if statement for stripe line
    }
  # End of loop through layers
  }
# End of layer_clouds if statement
}

# Third add in climate ribbon if required
if(climate_ribbon) {
  for(k in 1:length(conf)) {
    curr_hists = x$histories[,,dim]
    clim_quantile = apply(curr_hists,2,'quantile',probs=c((1-conf[k])/2,conf[k]+(1-conf[k])/2))
    polygon(c(x$time_grid,rev(x$time_grid)),c(clim_quantile[1,],rev(clim_quantile[2,])),col=col_ribbon,border=col_ribbon)
  }
}

# Finally add most representative plots if required
if(most_representative>0) {
  if(most_representative%%2==0) stop("For legal reasons number of most representative lines must be odd")
  # Function to create the most representative line
  rep_line = function(z,num_lines) {
    med = apply(z,2,median)
    diffs = rowSums(sweep(z,2,med,'-')^2)
    middle_val = round(length(diffs)/2,0)
    choose_vals = middle_val:(middle_val+num_lines-1)-(num_lines-1)/2
    return(z[order(diffs)[choose_vals],])
  }
  plot_lines = matrix(rep_line(x$histories[,,dim],num_lines=most_representative),nrow=most_representative)
  for(i in 1:most_representative) {
    lines(x$time_grid,plot_lines[i,],col=col_representative,lwd=2)
  }
}

}