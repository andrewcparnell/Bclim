\name{BclimCompile}
\alias{BclimCompile}
\title{
Compiles results from other Bclim stages
}
\description{
This function is designed to compile results from other Bclim stages. 
}
\usage{
BclimCompile(Layers, Mixtures, MCMC, Interpolations, core.name = "Core")
}
\arguments{
  \item{Layers}{
Output from the function \code{\link{BclimLayer}}
}
  \item{Mixtures}{
Output from the function \code{\link{BclimMixSer}} or \code{\link{BclimMixPar}}
}
  \item{MCMC}{
Output from the function \code{\link{BclimMCMC}}
}
  \item{Interpolations}{
Output from the function \code{\link{BclimInterp}}
}
  \item{core.name}{
A character string giving the name of the core
}
}
\details{
A Bclim run can be separated into 4 stages: 1 - Turning each individual pollen layer into climate marginal data posteriors, 2 - Approximating the marginal data posteriors as mixtures of Gaussians, 3 - Running an MCMC stage to estimate climates and climate volatilities, 4 - Interpolating the climates and volatilities onto a regular grid. This function takes the output from each of these functions and returns a Bclim list object that can be further manipulated.
}
\value{
A list of class Bclim is returned with the following components:
  \item{time.grid}{The time grid used for interpolation}
  \item{core.name}{The name of the core}
  \item{clim.interp}{The interpolated values of climate for the 3 climate dimensions}
  \item{vol.interp}{The interpolated volatilities for the 3 climate dimensions}
  \item{MDP}{The marginal data posteriors for each layer in the core}
  \item{ScMean}{The raw climate means (used for standardisation purposes)}
  \item{ScVar}{The raw climate variances (used for standardisation purposes)}
  \item{clim.dims}{The names of the different climate dimensions}
  \item{n}{The number of layers in the core}
  \item{m}{The number of climate dimensions (always 3)}
  \item{n.samp}{The number of samples created for the MDPs}
  \item{Chronsfile}{The location of the chronologies file}
  \item{nchron}{The number of chronologies}
  \item{chron.store}{The chronologies used for creating the posterior distributions of climate}
}
\references{
See Arxiv paper at http://arxiv.org/abs/1206.5009.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}

\seealso{
The output here can be used as an input to \code{\link{plotBclim}} and \code{\link{plotBclimVol}}. For advanced use of Bclim, see the functions referenced at the top of this file.
}
\examples{
\dontrun{
# Set the working directory using setwd (not shown)

# Download and load in the response surfaces:
url1 <- 'http://mathsci.ucd.ie/~parnell_a/required.data3D.RData'
download.file(url1,'required_data3D.RData')

# and now the pollen
url2 <- 'http://mathsci.ucd.ie/~parnell_a/SlugganPollen.txt'
download.file(url2,'SlugganPollen.txt')

# and finally the chronologies
url3 <- 'http://mathsci.ucd.ie/~parnell_a/Sluggan_2chrons.txt'
download.file(url3,'Slugganchrons.txt')

# Create variables which state the locations of the pollen and chronologies
pollen.loc <- paste(getwd(),'/SlugganPollen.txt',sep='')
chron.loc <- paste(getwd(),'/Slugganchrons.txt',sep='')

# Load in the response surfaces
load('required.data3D.RData')

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
