\name{BclimLayer}
\alias{BclimLayer}
\title{
Function to approximate pollen layers as marginal data posteriors
}
\description{
This function takes a set of pollen data slice by slice and turns it into marginal data posteriors by means of a zero-inflated Negative Binomial model and some response surfaces.
}
\usage{
BclimLayer(pollen.loc, required.data3D, nsamples = 1000)
}
\arguments{
  \item{pollen.loc}{
A character string detailing the location of the 28 taxa pollen file. The 28 taxa should be in the following order: Abies Alnus Betula Carpinus Castanea Cedrus Corylus Ephedra Fagus Juniperus Larix Olea Ostrya Phillyrea Picea Pinus.D Pinus.H Pistacia Quercus.D Quercus.E Salix Tilia Ulmus Artemisia Chenopodiaceae Cyperaceae Ericales Gramineae
}
  \item{required.data3D}{
A list object of response surfaces. A suitable list object can be downloaded as in the example below. You need to load this in before this part of Bclim can be run.
}
  \item{nsamples}{
The number of samples of each marginal data posterior to take. You are unlikely to need to change this from the default value of 1000.
}
}
\details{
A marginal data posterior (MDP) is of the form pi(ci|yi) where pi is a probability distribution, ci is the 3D climate at depth di and yi is the 28-dimensional pollen vector. This function loops through each layer in the core to produce MDPs which represent the information about climate obtained only from that layer of pollen. 

This function could be parallelised to enable speedier computation.
}
\value{
The output is a nsamples x n x m array where nsamples is as above, n is the number of layers and m is the number of climate dimensions (always 3).
}
\references{
Fore more detail on the algorithm see:
Salter-Townshend, M. and J. Haslett (2012). Fast Inversion of a Flexible Regression Model for Multivariate, Zero-Inflated Pollen Counts. Environmetrics.
Sweeney, J. (2012). Advances in Bayesian Model Development and Inversion in Multivariate Inverse Inference Problems with application to palaeoclimate reconstruction. Ph. D. thesis, Trinity College Dublin.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie> and James Sweeney
}
\seealso{
The main Bclim function is \code{\link{BclimRun}}. See also the other 3 stages: \code{\link{BclimInterp}}, \code{\link{BclimMixSer}} (or \code{\link{BclimMixPar}}), and \code{\link{BclimMCMC}}, 
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

