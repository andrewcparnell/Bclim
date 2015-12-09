\name{BclimMixSer}
\alias{BclimMixSer}
\title{
Serial version of Bclim mixture analysis
}
\description{
Function to approximate marginal data posteriors as mixtures of Gaussians
}
\usage{
BclimMixSer(MDP, G = 10, mixwarnings = FALSE)
}
\arguments{
  \item{MDP}{
A set of marginal data posteriors, as produced by \code{\link{BclimLayer}}
}
  \item{G}{
The number of Gaussian groups required for each layer to be partitioned into. The default of 10 is usually fine.
}
  \item{mixwarnings}{
Whether to suppress mixture warnings (default) or not.
}
}
\details{
This function approximates marginal data posteriors (MDPs) as mixtures of Gaussians. The mixture algorithm is taken from the Mclust package which is a required installation for this to run. This is the serial version, i.e. it only uses one processor, as opposed to the \code{\link{BclimMixPar}} parallel version which will run much faster but requires extra packages to be installed and a multi-core machine.
}
\value{
Outputs a list containing the following objects:
  \item{MDP }{A nsamples x n x m array (these values are described below) }
  \item{n }{The number of layers}
  \item{m }{The number of climate dimensions (always 3)}
  \item{n.samp }{The number of samples given in \code{\link{BclimLayer}}}
  \item{ScMean}{The raw climate means (used for standardisation purposes)}
  \item{ScVar}{The raw climate variances (used for standardisation purposes)}
  \item{G }{The number of mixture groups (as above)}
  \item{mu.mat }{An estimate of the Gaussian mixture mean components}
  \item{tau.mat }{An estimate of the Gaussian mixture precision components}
  \item{p.mat }{An estimate of the Gaussian mixture proportions} 
}
\references{
See Arxiv paper at http://arxiv.org/abs/1206.5009.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}

\seealso{
The output here can be used as an input to \code{\link{BclimMCMC}}. See the main \code{\link{BclimRun}} function for more details of the other stages you might need to run.
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
