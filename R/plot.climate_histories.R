plot.climate_histories = function(x,dim=1,slice_clouds=TRUE,chron=NULL,climate_ribbon=TRUE,most_representative=1,conf=c(0.95,0.75,0.5), col_clouds = grDevices::rgb(0,0,1,0.2), col_ribbon=grDevices::rgb(1,0,0,0.4),col_representative=grDevices::rgb(0,1,0),present_left=TRUE,...) {

# Get extra arguments if provided
ex = list(...)

clim_names = dimnames(x$slice_clouds$slice_clouds)[[3]]
if(is.null(ex$xlab)) ex$xlab = 'Age (k cal years BP)'
if(is.null(ex$ylab)) ex$ylab = clim_names[dim]

# First create base plot
x_range = range(x$time_grid)
if(!present_left) x_range = rev(x_range)
y_range = range(x$slice_clouds$slice_clouds[,,dim])

if(is.null(ex$xlim)) ex$xlim = x_range
if(is.null(ex$ylim)) ex$ylim = y_range
ex$xaxt = 'n'
ex$yaxt = 'n'
ex$x = ex$y = 1
ex$type = 'n'
if(is.null(ex$las)) ex$las = 1

args = utils::modifyList(ex, list(...))
do.call("plot", args)

graphics::axis(side=1,at=pretty(x$time_grid,n=10))
graphics::axis(side=2,at=pretty(x$slice_clouds$slice_clouds[,,dim],n=10), las = ex$las)
graphics::grid()

# Get some objects to make coding neater
n_slices = x$slice_clouds$n_slices
n_samples = x$slice_clouds$n_samples

# Second add in slice clouds if required
if(slice_clouds) {
  if(is.null(chron)) stop("A chronology is required for plotting slice clouds")
  for(i in 1:n_slices) {
    # Get current MDPs and a suitable number of chronologies to match
    curr_MDP = x$slice_clouds$slice_clouds[,i,dim]
    if(n_samples>nrow(chron)) {
      curr_chron = sample(chron[,i],size=n_samples,replace=TRUE)
    } else {
      curr_chron = chron[1:n_samples,i]
    }
    
    if(stats::var(curr_chron)>0) {
      dens = MASS::kde2d(curr_chron,curr_MDP)
      
      # Standardise and find 95% limit
      dens$z = dens$z/sum(dens$z)
      
      # Turn dens into a vector and find the highest values
      de = as.vector(dens$z)
      do = order(de)
      cu = cumsum(de[do])
      
      for(k in 1:length(conf)) {
        # Find which ones are above the threshold
        good_cu = which(cu>1-conf[k])
        good_clim_vec = sort((1:length(de))[do][good_cu])
        
        good_clim = 1*matrix((1:length(de)) %in% good_clim_vec, 
                           ncol=ncol(dens$z), 
                           nrow=nrow(dens$z))
        curr_dens = dens
        curr_dens$z = good_clim
        cont = grDevices::contourLines(curr_dens, nlevels = 1)
        graphics::polygon(cont[[1]]$x,cont[[1]]$y,col=col_clouds,border=col_clouds)
      }
      
      limits_diff = (max(dens$z)-min(dens$z))/100

      # Loop through confidence levels to plot
      for(k in 1:length(conf)) {
        limits = min(dens$z)
        prop = 1
        while(prop>conf[k]) {
          limits = limits+limits_diff
          prop = sum(dens$z[dens$z>limits])
        }
      }
      
    # End of if statement for variance of chronology
    } else {
      # If there's no variance so it's a vertical stripe
      for(k in 1:length(conf)) {
        graphics::lines(c(curr_chron[1],curr_chron[1]),stats::quantile(curr_MDP,probs=c((1-conf[k])/2,conf[k]+(1-conf[k])/2)),col=col_clouds)
      }
    # End of if statement for stripe line
    }
  # End of loop through slices
  }
# End of slice_clouds if statement
}

# Third add in climate ribbon if required
if(climate_ribbon) {
  for(k in 1:length(conf)) {
    curr_hists = x$histories[,,dim]
    clim_quantile = apply(curr_hists,2,stats::quantile,probs=c((1-conf[k])/2,conf[k]+(1-conf[k])/2))
    graphics::polygon(c(x$time_grid,rev(x$time_grid)),c(clim_quantile[1,],rev(clim_quantile[2,])),col=col_ribbon,border=col_ribbon)
  }
}

# Finally add most representative plots if required
if(most_representative>0) {
  if(most_representative%%2==0) stop("For legal reasons number of most representative lines must be odd")
  # Function to create the most representative line
  rep_line = function(z,num_lines) {
    med = apply(z,2,stats::median)
    diffs = rowSums(sweep(z,2,med,'-')^2)
    middle_val = round(length(diffs)/2,0)
    choose_vals = middle_val:(middle_val+num_lines-1)-(num_lines-1)/2
    return(z[order(diffs)[choose_vals],])
  }
  plot_lines = matrix(rep_line(x$histories[,,dim],num_lines=most_representative),nrow=most_representative)
  for(i in 1:most_representative) {
    graphics::lines(x$time_grid,plot_lines[i,],col=col_representative,lwd=2)
  }
}

}