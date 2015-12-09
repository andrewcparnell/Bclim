\name{BclimInterp}
\alias{BclimInterp}
\title{
Interpolation function for producing palaeoclimate histories on a grid
}
\description{
The function takes output from previous stages of Bclim and produces interpolated values of palaeoclimate (in 3 dimensions) and palaeoclimate volatility.
}
\usage{
BclimInterp(Bclim.data, Bclim.res, tgrid = seq(0, 14, by = 0.1))
}
\arguments{
  \item{Bclim.data}{
Output from the function \code{\link{BclimMixSer}} or \code{\link{BclimMixPar}}
}
  \item{Bclim.res}{
Output from the function \code{\link{BclimMCMC}}
}
  \item{tgrid}{
A grid (not necessarily regular) on which to interpolate climate
}
}
\details{
A Bclim run can be separated into 4 stages: 1 - Turning each individual pollen layer into climate marginal data posteriors, 2 - Approximating the marginal data posteriors as mixtures of Gaussians, 3 - Running an MCMC stage to estimate climates and climate volatilities, 4 - Interpolating the climates and volatilities onto a regular grid. This function runs stage 4.

}
\value{
Returns a list containing the following elements:
  \item{clim.interp}{Interpolated values of 3D climate}
  \item{vol.interp}{Interpolated values of 3D climate volatility}
  \item{time.grid}{The time grid used for interpolation}
}
\references{
See Arxiv paper at http://arxiv.org/abs/1206.5009.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
The main Bclim function is \code{\link{BclimRun}}. See also the other 3 stages: \code{\link{BclimLayer}}, \code{\link{BclimMixSer}} (or \code{\link{BclimMixPar}}), and \code{\link{BclimMCMC}}, 
}
\examples{
## See help for BclimMCMC and webpage at http://mathsci.ucd.ie/~parnell_a/Bclim.html
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }

