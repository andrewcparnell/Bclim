\name{BclimMCMC}
\alias{BclimMCMC}

\title{
Bclim Markov chain Monte Carlo run
}
\description{
A Bclim Markov chain Monte Carlo (MCMC) run to determine the volatilities and climate from a set of marginal data posteriors approximated as mixtures.
}
\usage{
BclimMCMC(Bclimdata,chron.loc,nchron=10000,control.mcmc=list(iterations=100000,
burnin=20000,thinby=40,report=100),control.chains=list(v.mh.sd=2,phi1.mh.sd=1,
phi2.mh.sd=10,vstart=statmod::rinvgauss(Bclimdata$n-1,2,1),
Zstart=sample(1:Bclimdata$G,Bclimdata$n,replace=TRUE),phi1start=rep(3,Bclimdata$m),
phi2start=rep(20,Bclimdata$m)),control.priors=list(phi1dlmean=rep(1.275,Bclimdata$m),
phi1dlsd=rep(0.076,Bclimdata$m),phi2dlmean=rep(4.231,Bclimdata$m),
phi2dlsd=rep(0.271,Bclimdata$m)))
}
\arguments{
  \item{Bclimdata}{
A set of approximated marginal data posteriors (MDPs) taken from a run of \code{\link{BclimMixPar}} or \code{\link{BclimMixSer}}
}
  \item{chron.loc}{
A character string giving the location of the chronologies file
}
  \item{nchron}{
The number of chronologies in the chronologies file
}
  \item{control.mcmc}{
A list containing elements that control the MCMC, including the number of iterations, the size of the burn-in period, the amount to thinby, and how often for the algorithm to report its progress.
}
  \item{control.chains}{
A list containing elements that control the starting values of the parameters (vstart, Zstart, phi1start and phi2) and the Metropolis-Hastings proposal standard deviation for v, phi1 and phi2.
}
  \item{control.priors}{
A list containing the prior parameters for the volatilities, given by phi1 and phi2, both of which should be the log-mean and log-sd of the .
}
}
\details{
This function takes the MDPs, approximated as Gaussian mixtures from \code{\link{BclimMixPar}} or \code{\link{BclimMixSer}}, and produces estimated climate and climate volatilities using an MCMC algorithm. The full details are in the Arxiv paper referenced below. The options listed above allow quite a detailed level of control over the behaviour of the algorithm, and convergence should be checked using suitable means (see e.g. the R package boa or coda).
}
\value{
A list object wiht the following elements
  \item{v.store }{Samples of the posterior estimated volatilities}
  \item{chron.store }{Samples of the used chronologies}
  \item{c.store }{Samples of the posterior estimated climates}
  \item{z.store }{Samples of the posterior mixture indices}
  \item{phi1 }{Values used for the IG prior on v for each climate dimension}
  \item{phi2 }{Values used for the IG prior on v for each climate dimension}
  \item{chron.loc }{A character string giving the location of the chronology file}
  \item{nchron }{The number of chronologies in the chronology file}
}
\references{
See Arxiv paper at http://arxiv.org/abs/1206.5009.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}

\seealso{
The main Bclim function is \code{\link{BclimRun}}. See also the other 3 stages: \code{\link{BclimInterp}}, \code{\link{BclimMixSer}} (or \code{\link{BclimMixPar}}), and \code{\link{BclimLayer}}, 

}
\examples{
\dontrun{
# Set the working directory using setwd (not shown)

# Download and load in the response surfaces:
url1 <- 'http://mathsci.ucd.ie/~parnell_a/media/requireddata3D.RData'
download.file(url1,'required_data3D.RData')

# and now the pollen
url2 <- 'http://mathsci.ucd.ie/~parnell_a/media/SlugganPollen.txt'
download.file(url2,'SlugganPollen.txt')

# and finally the chronologies
url3 <- 'http://mathsci.ucd.ie/~parnell_a/media/Sluggan_2chrons.txt'
download.file(url3,'Slugganchrons.txt')

# Create variables which state the locations of the pollen and chronologies
pollen.loc <- paste(getwd(),'/SlugganPollen.txt',sep='')
chron.loc <- paste(getwd(),'/Slugganchrons.txt',sep='')

# Load in the response surfaces
load('required_data3D.RData')

## note that all of these functions have further options you can change with
step1 <- BclimLayer(pollen.loc,required.data3D=required.data3D)
step2 <- BclimMixSer(step1) 
# See also the parallelised version BclimMixPar if you have doMC and foreach installed
step3 <- BclimMCMC(step2,chron.loc) 
# You should probably do some convergence checking after this step
step4 <- BclimInterp(step2,step3) 
results <- BclimCompile(step1,step2,step3,step4,core.name="Sluggan Moss")

# Create a plot of MTCO (dim=2)
plotBclim(results,dim=2)

# Create a volatility plot
plotBclimVol(results,dim=2)
}
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
