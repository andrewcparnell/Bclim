\name{BclimRun}
\alias{BclimRun}
\title{
Run all stages of Bclim together
}
\description{
This function is intended for beggining/light users of Bclim who wish to create palaeoclimate reconstructions from their own pollen data. It collects together the functions \code{\link{BclimLayer}}, \code{\link{BclimMixSer}}, \code{\link{BclimMCMC}}, \code{\link{BclimInterp}} and \code{\link{BclimCompile}}.
}
\usage{
BclimRun(pollen.loc, chron.loc, core.name = "Core", time.grid = seq(0, 14, length = 100), 
  required.data3D = NULL, nchrons = 10000, parallel = FALSE, save.as.you.go = TRUE)
}
\arguments{
  \item{pollen.loc}{
A character string detailing the location of the 28 taxa pollen file. The 28 taxa should be in the following order: Abies Alnus Betula Carpinus Castanea Cedrus Corylus Ephedra Fagus Juniperus Larix Olea Ostrya Phillyrea Picea Pinus.D Pinus.H Pistacia Quercus.D Quercus.E Salix Tilia Ulmus Artemisia Chenopodiaceae Cyperaceae Ericales Gramineae
}
  \item{chron.loc}{
A character string detailing the location of the chronologies file. This file should contain a matrix with the same number of columns as the number of slices in the core (i.e. the number of rows of the pollen data) and at least 2000 rows. Such a file can be created via the Bchron R package from radiocarbon and non-radiocarbon data.
}
  \item{core.name}{
A character string giving the name of the core. Defaults to 'Core'.
}
  \item{time.grid}{
A sequence of values of the time variable over which climates and volatilities are to be predicted. Defaults to 0 to 14ka BP with centurial time steps.
}
  \item{required.data3D}{
A list object of response surfaces. A suitable list object can be downloaded as in the example below. You need to load this in before this part of Bclim can be run.
}
  \item{nchrons}{
The number of chronologies in the chronologies file. When using output from Bchron this can be safely left at 10,000.
}
  \item{parallel}{
Whether to run the model using parallel processing technology. To set this as true you must have installed the doMC and foreach R packages.
}
  \item{save.as.you.go}{
Whether to save output from the differing stages. You should probably keep this set to TRUE in case any of them fail for some reason.
}
}
\details{
A Bclim run can be separated into 4 stages: 1 - Turning each individual pollen layer into climate marginal data posteriors, 2 - Approximating the marginal data posteriors as mixtures of Gaussians, 3 - Running an MCMC stage to estimate climates and climate volatalities, 4 - Interpolating the climates and volatilities onto a regular grid. This function calls all of these stages in turn and will tidy them up for presentation or use in further analysis, such as for plotting using \code{\link{plotBclim}}. Climate is always 3 dimensional and represented as Growing Degree Days Above 5C (GDD5), The Mean Temperature of the Coldest Month (MTCO) and the ratio of Actual to Potential Evapotranspiration (AET/PET).
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
# First set the working directory using the setwd command (not shown)

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

# Now run Bclim the simple way
RunSluggan <- BclimRun(pollen.loc,chron.loc,core.name="Sluggan Moss",time.grid=seq(0,14,length=100),
  nchrons=10000,required.data3D=required.data3D,parallel=FALSE)

# Create some plots of climate
plotBclim(RunSluggan,dim=1)

# Create some plots of the volatilities
plotBclimVol(RunSluggan,dim=3)
}
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
