\name{climate_histories}
\alias{climate_histories}
\title{
Create Bclim climate_histories
}
\description{
Runs a number of algorithms to create climate histories for a given set of layer clouds (from \code{\link{layer_clouds}} and a set of chronologies
}
\usage{
climate_histories(layer_clouds,
chronology,
time_grid,
n_mix=10,
mix_warnings=FALSE,
n_chron=2000,
keep_parameters=TRUE,
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
                    phi2dl_sd=rep(0.271,layer_clouds$n_dimensions)))
}
\arguments{
  \item{layer_clouds}{
An object of class \code{layer_clouds} obtained from \code{\link{layer_clouds}}
}
  \item{chronology}{
A set of chronologies given as a matrix. These should be provided in thousands of years before present. See details below
}
  \item{time_grid}{
The time grid on which to create the climate histories
}
  \item{n_mix}{
The number of mixture components for the Mclust mixture algorithm.
}
  \item{mix_warnings}{
Whether to display warnings related to the mixture algorithm
}
  \item{n_chron}{
The number of chronologies to use
}
  \item{keep_parameters}{
Whether to keep latent parameters or not. Useful for convergence checking so default is TRUE
}
  \item{control_mcmc}{
A list containing elements that control the MCMC, including the number of iterations, the size of the burn-in period, the amount to thinby, and how often for the algorithm to report its progress.
}
  \item{control_chains}{
A list containing elements that control the starting values of the parameters (v_start, Z_start, phi1_start and phi2_start) and the Metropolis-Hastings proposal standard deviation for v, phi1 and phi2.
}
  \item{control_priors}{
A list containing the prior parameters for the volatilities, given by phi1 and phi2, both of which should be the log-mean and log-sd of the log-normal distribution. The values provided here are for the GISP2 ice core for the period 0 to 10k years BP.
}
}
\details{
This function takes the layer_clouds produced by \code{\link{layer_clouds}} uses a set of algorithms to produce cliamte histories on the provided time grid. The full details are in the paper referenced below. The options listed above allow quite a detailed level of control over the behaviour of the algorithm, and convergence should be checked using suitable means (see e.g. the R package boa or coda).

One of the key inputs to this function is a chronology. This should be a matrix of n_chron by n_layers containing sample chronologies as produced by, e.g. the R package Bchron. These are used by the \code{climate_histories} function to take account of chronological uncertainty. In the (unlikely) event that there is no chronological uncertainty, the rows of the chronologies can be identical.

}
\value{
A list object with the following elements
  \item{v.store }{Samples of the posterior estimated volatilities}
  \item{chron.store }{Samples of the used chronologies}
  \item{c.store }{Samples of the posterior estimated climates}
  \item{z.store }{Samples of the posterior mixture indices}
  \item{phi1 }{Values used for the IG prior on v for each climate dimension}
  \item{phi2 }{Values used for the IG prior on v for each climate dimension}
  \item{chron.loc }{A character string giving the location of the chronology file}
  \item{nchron }{The number of chronologies in the chronology file}
  \item{parameters }{A list containing further latent parameter values for convergence checking (only if \code{keep_parameters} is TRUE)}
}
\references{
Parnell, A. C., et al. (2015), Bayesian inference for palaeoclimate with time uncertainty and stochastic volatility. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 115–138.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie> and Thinh Doan
}

\seealso{
\code{\link{layer_clouds}} for producing the input for this function. See \code{\link{plot.climate_histories}} and  \code{\link{summary.climate_histories}} for plotting and summary details

For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the authors personal webpage (http://mathsci.ucd.ie/~parnell_a/Bclim.html)?
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
