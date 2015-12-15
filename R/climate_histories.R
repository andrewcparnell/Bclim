climate_histories = function(layer_clouds,
                             chronology,
                             time_grid,
                             n_mix=10,
                             mix_warnings=FALSE,
                             n_chron=2000,
                             control_mcmc=list(iterations=100000,
                                               burnin=20000,
                                               thinby=40,
                                               report=100),
                             control_chains=list(v_mh_sd=2,
                                                 phi1_mh_sd=1,
                                                 phi2_mh_sd=10,
                                                 v_start=statmod::rinvgauss(layer_clouds$n_layers-1,2,1),
                                                 Z_start=sample(1:n_mix,
                                                                layer_clouds$n_layers,
                                                               replace=TRUE),
                                                 phi1_start=rep(3,layer_clouds$n_dimensions),
                                                 phi2_start=rep(20,layer_clouds$n_dimensions)),
                             control_priors=list(phi1_dl_mean=rep(1.275,layer_clouds$n_dimensions),
                                                 phi1_dl_sd=rep(0.076,layer_clouds$n_dimensions),
                                                 phi2_dl_mean=rep(4.231,layer_clouds$n_dimensions),
                                                 phi2dl_sd=rep(0.271,layer_clouds$n_dimensions))) {

################# USEFUL FUNCTIONS #################
  
NIGB = function(delta, IG.sum, tb.points, c.start, c.end){
  total.points = length(tb.points)
  NIG.bridge = rep(NaN, nrow=total.points)
  NIG.bridge[1] = c.start
  NIG.bridge[total.points] = c.end
  IG.increment = rep(NaN, total.points-2)
  z = IG.sum
  eps=1E-5
  
  l = 2
  while(l < total.points){
    if((tb.points[l]-tb.points[l-1])<eps) {
      IG.increment[l-1] = 0
    } else {
      
      q = rnorm(1)^2  # Generate chi-square random variate
      
      # Reparameterise
      mu = (tb.points[total.points]-tb.points[l]) / (tb.points[l]-tb.points[l-1])
      lambda = (delta^2 * (tb.points[total.points]-tb.points[l])^2) / z
      
      if(z==0) {
        print(c(delta,IG.sum,tb.points,c.start,c.end))
        stop("Problem in Inverse Gaussian Bridge - IG sum is zero")
      }
      
      # Compute the roots of the chi-square random variate
      s1 = mu + (mu^2*q)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*q + mu^2*q^2)
      if(lambda < eps) { s1 = mu }
      s2 = mu^2 / s1
      
      # Acceptance/rejection of root
      s = ifelse(runif(1) < mu*(1+s1)/((1+mu)*(mu+s1)), s1, s2)
      
      IG.increment[l-1] = z / (1 + s)  # Compute the IG incrrement
      
      if(any(IG.increment<0, na.rm=TRUE)) stop("Inverse Gaussian bridge producing negative variances")
    } # End of if statement
    
    # Rescale the sum of left-over distance of the IG bridge
    z = z - IG.increment[l-1]
    
    #Compute the current value of the NIG bridge
    NIG.bridge[l] = (c.start*z + c.end*(IG.sum-z)) / IG.sum + rnorm(1) * (IG.sum-z)*z / IG.sum
    l = l + 1
  }
  return(list(IGB = c(IG.increment, (IG.sum - sum(IG.increment))), NIGB = t(NIG.bridge)))
}
  
# Extrapolation function
NIGextrap = function(curr.clim, curr.chron, tg.select, phi1, phi2, future=FALSE) {
  t.diff = abs(diff(sort(c(tg.select, curr.chron))))
  mu = phi1*t.diff # Re-parameterise on to the R version of the IG distribution
  lambda = phi1*phi2*t.diff^2
  v.out = statmod::rinvgauss(length(t.diff),mu,lambda)
  if(future==TRUE) {
    return(list(NIG = rev(cumsum(rev(rnorm(length(t.diff), mean=0, sd=sqrt(v.out))))) + curr.clim, IG=v.out))
  } else {
    return(list(NIG=cumsum(rnorm(length(t.diff), mean=0, sd=sqrt(v.out)))+ curr.clim, IG=v.out))
  }
  
}  
  
################# MIXTURE ESTIMATION #################

# Calculate n.samp = number of samples, n = number of layers, m = number of climate dimensions
n_samples = layer_clouds$n_samples
n_layers = layer_clouds$n_layers
n_dimensions = layer_clouds$n_dimensions

scale_mean = rep(0,n_dimensions)
scale_var = rep(1,n_dimensions)
MDP = layer_clouds$layer_clouds
for(i in 1:n_dimensions) {
  scale_mean[i] = mean(layer_clouds$layer_clouds[,,i])
  scale_var[i] = median(diag(var(layer_clouds$layer_clouds[,,i])))
  MDP[,,i] = (layer_clouds$layer_clouds[,,i]-scale_mean[i])/sqrt(scale_var[i])
}

# Set up mixture components
mean_mat = array(NA,dim=c(n_layers,n_dimensions,n_mix))
prec_mat = array(NA,dim=c(n_layers,n_dimensions,n_mix))
prop_mat = matrix(NA,nrow=n_layers,ncol=n_mix)

ans.all = list()
cat('Mixture estimation:\n')
for(i in 1:n_layers) {
  cat("\r")
  cat(format(round(100*i/n_layers,2), nsmall = 2),"% completed",sep='')
  ans.all[[i]] = mclust::Mclust(MDP[,i,],G=n_mix,modelNames="EII",warn=mix_warnings)
}
cat('\n')

for(i in 1:n_layers) {
  mean_mat[i,,] = ans.all[[i]]$parameters$mean
  for(g in 1:n_mix) {
    prec_mat[i,,g] = 1/diag(ans.all[[i]]$parameters$variance$sigma[,,g])
  }
  prop_mat[i,] = ans.all[[i]]$parameters$pro
  # End of loop through n layers
}

################# MCMC #################

# Create output matrices
remaining = (control_mcmc$iterations-control_mcmc$burnin)/control_mcmc$thinby
if(remaining!=as.integer(remaining))
  stop("Iterations minus burnin divided by thinby must be an integer")

vout = rep(0,length=n_dimensions*(n_layers-1)*remaining)
zout = rep(0,length=n_dimensions*n_layers*remaining)
chronout = rep(0,length=n_layers*remaining)
cout = rep(0,length=n_dimensions*(n_layers)*remaining)
phi1out = phi2out = rep(0,length=n_dimensions*remaining)

# Re-dim the precisions matrix
Bclimprec = prec_mat[,1,]

# Write out the chronologies to a temporary file
chron_loc = paste0(getwd(),'/chron.txt')
if(any(chronology>1000)) warning("Ensure chronologies are provided in thousands of years.")
write.table(chronology[1:n_chron,],file=chron_loc,row.names=FALSE,col.names=FALSE,quote=FALSE)

# Run C code
out = .C("BclimMCMC3D",
          as.integer(n_mix),
          as.integer(n_layers),
          as.integer(n_dimensions),
          as.integer(n_chron),
          as.double(prop_mat),
          as.double(mean_mat),
          as.double(Bclimprec),
          as.character(chron_loc),
          as.integer(control_mcmc$iterations),
          as.integer(control_mcmc$burnin),
          as.integer(control_mcmc$thinby),
          as.integer(control_mcmc$report),
          as.double(control_chains$v_mh_sd),
          as.double(control_chains$v_start),
          as.integer(control_chains$Z_start),
          as.double(control_chains$phi1_start),
          as.double(control_chains$phi2_start),
          as.double(vout),
          as.integer(zout),
          as.double(chronout),
          as.double(cout),
          as.double(phi1out),
          as.double(phi2out),
          as.double(control_priors$phi1_dl_mean),
          as.double(control_priors$phi1_dl_sd),
          as.double(control_priors$phi2_dl_mean),
          as.double(control_priors$phi2_dl_sd),
          as.double(control_chains$phi1_mh_sd),
          as.double(control_chains$phi2_mh_sd)
)

vout  = array(NA,dim=c(remaining,n_layers-1,n_dimensions))
cout  = array(NA,dim=c(remaining,n_layers,n_dimensions))
for(i in 1:remaining) {
  for(j in 1:n_dimensions) {
    vout[i,,j] = out[[18]][seq(1,n_layers-1)+(j-1)*(n_layers-1)+(i-1)*(n_layers-1)*n_dimensions]
    cout[i,,j] = out[[21]][seq(1,n_layers)+(j-1)*(n_layers)+(i-1)*(n_layers)*n_dimensions]
  }
}
chronout = matrix(out[[20]],ncol=n_layers,nrow=remaining,byrow=TRUE)
zout = matrix(out[[19]],ncol=n_layers,nrow=remaining,byrow=TRUE)
phi1out = matrix(out[[22]],ncol=n_dimensions,nrow=remaining,byrow=TRUE)
phi2out = matrix(out[[23]],ncol=n_dimensions,nrow=remaining,byrow=TRUE)

###### Stage 2 - Interpolation

chron = chronout  # n.s-by-n where n.s is the numer of iterations
phi1 = phi1out #n.s-by-m
phi2 = phi2out
gamma = sqrt(phi2 / phi1) #n.s-by-m
delta = sqrt(phi1 * phi2)
clim = cout   # n.s-by-n-by-m
v = vout  # n.s-by-(n-1)-by-m squared volatility
n.s = dim(clim)[1]
n = dim(clim)[2]
m = dim(clim)[3]
n.g = length(time_grid)
eps = 1e-5

#	Create some arrays to hold the interpolated climates and squared volatilities v
# Note that the storage here is longer as we have to interpolate onto a finer grid which includes the chronologies
clim.interp = array(NaN, dim = c(n.s, n.g, m))
v.interp.all = array(NA, dim = c(n.s , n.g-1+n, m))
v.interp = array(NaN, dim = c(n.s, n.g-1, m))

for (j in 1:m) {  #loop through clim dim
  cat("Interpolation for climate dimension ", j,":\n",sep='')

  for (i in 1:n.s) { #loop through each sample
    if(i%%10==0) {
      cat("\r")
      cat("",format(round(100*i/n.s,2), nsmall = 2),"% completed",sep='')
      if(i<n.s) {
        cat("\r")
      } else {
        cat("\n")
      }
    }

    # Get small vectors of the current values
    curr.clim = clim[i,,j] # length n
    curr.v = v[i,,j] # length n-1
    curr.chron = chron[i,] # length n
    time_grid.all = sort(c(curr.chron, time_grid)) # length n + n.g

    # If there are any chronology times which exactly match the time_grid then fill in those climates and variances as zero
    if(any(diff(time_grid.all)<eps)) {
      curr.chron = curr.chron+eps  # shift the chronology by a little bit
      time_grid.all = sort(c(curr.chron, time_grid))
    }

    # if extrapolation into the future is required
    if(time_grid[1] < min(curr.chron)) {
      # Get all the time_grid that are less than the smallest value of the chronology
      t.select.index = which(time_grid < min(curr.chron))
      t.all.select.index = which(time_grid.all < min(curr.chron))
      t.select = time_grid[t.select.index]

      clim.temp = NIGextrap(curr.clim[1], curr.chron[1], t.select, phi1[i,j], phi2[i,j], future=TRUE)
      clim.interp[i, t.select.index, j] = clim.temp$NIG
      v.interp.all[i, t.all.select.index, j] = rev(cumsum(clim.temp$IG))
    }

    # if extrapolation into the past is required
    if(time_grid[length(time_grid)] > max(curr.chron)) {
      t.select.index = which(time_grid>max(curr.chron))
      t.all.select.index = which(time_grid.all > max(curr.chron))
      t.select = time_grid[t.select.index]

      # Perform extrapolation
      clim.temp = NIGextrap(curr.clim[length(curr.clim)], curr.chron[length(curr.clim)],
                            t.select, phi1[i,j], phi2[i,j])
      clim.interp[i, t.select.index, j] = clim.temp$NIG
      v.interp.all[i, t.all.select.index-1, j] = cumsum(clim.temp$IG)
    }

    # Now look at how many gaps there are
    for(k in 1:(n-1)) {
      # Find out if in this section there are any time_grid points inside
      t.select.index = which(time_grid >= curr.chron[k] & time_grid < curr.chron[k+1])
      t.all.select.index = which(time_grid.all >= curr.chron[k] & time_grid.all < curr.chron[k+1])

      # If the length of t.select.index is positive then there are points inside and we should use the NIG bridge
      if(length(t.select.index) > 0) {
        # Select which bits of the grid are inbetween the current section of the chronology
        t.select = time_grid[t.select.index]

        # Now interpolate using the NIGB code
        clim.temp = NIGB(delta[i,j], curr.v[k], c(curr.chron[k],t.select,curr.chron[k+1]),
                          curr.clim[k], curr.clim[k+1])

        clim.interp[i,t.select.index,j] = clim.temp$NIGB[-c(1,length(clim.temp$NIGB))]
        v.interp.all[i,t.all.select.index,j] = clim.temp$IGB
      }

      # If there's no grid points just store the v value from the original data
      if(length(t.select.index)==0 & length(t.all.select.index)>0) {
        v.interp.all[i, t.all.select.index, j] = curr.v[k]
      }
    }

    # Sum over the interpolated v.interp values so as to get the correct interpolated v values
    time_grid.select.lower = match(time_grid, time_grid.all)
    time_grid.select.upper = c((time_grid.select.lower-1)[2:length(time_grid.select.lower)], length(time_grid.all))
    for(l in 1:(n.g-1)) {
      v.interp[i,l,j] = sum(v.interp.all[i, time_grid.select.lower[l]:time_grid.select.upper[l], j])
    }

    if(any(is.na(clim.interp[i,,j]))) browser()# stop("Some interpolated climates have been given NAs")
    if(any(is.na(v.interp[i,,j]))) browser()# stop("Some interpolated v values have been given NAs")

  } #end of sample loops
} #end climate dim loop

cat("Completed! \n")
clim.interp.resc = sweep(sweep(clim.interp,3,sqrt(scale_var),'*'),3,scale_mean,'+')
out = list(histories = clim.interp.resc, time_grid=time_grid, layer_clouds=layer_clouds)
class(out) = 'climate_histories'

return(out)

}