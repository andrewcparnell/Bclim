\name{Bclim-package}
\alias{Bclim-package}
\alias{Bclim}
\docType{package}
\title{
Create palaeoclimate reconstructions from pollen cores
}
\description{
This package takes pollen core and chronological data to produce palaeoclimate estimates with quantified uncertainties. The package runs in 4 stages, collected together in the function BclimCompile. 
}
\details{
\tabular{ll}{
Package: \tab Bclim\cr
Type: \tab Package\cr
Version: \tab 2.3.1\cr
Date: \tab 2013-07-29\cr
License: \tab GPL >=2 \cr
}
See the function \code{\link{BclimRun}} for a full description of the workings of the package and an example script for getting it to work.
}
\author{
Andrew Parnell, Thinh Doan and James Sweeney

Maintainer: Andrew Parnell <andrew.parnell@ucd.ie>
}
\references{
See Arxiv paper at http://arxiv.org/abs/1206.5009, see also references in \code{\link{BclimLayer}}. The webpage http://mathsci.ucd.ie/~parnell_a/Bclim.html contains more detail and background information.
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
\keyword{ package }
\seealso{
\code{\link{BclimRun}}, \code{\link{BclimLayer}}, \code{\link{BclimMixSer}}, \code{\link{BclimMixPar}}, \code{\link{BclimMCMC}}, \code{\link{BclimInterp}}, \code{\link{BclimCompile}}, \code{\link{plotBclim}}, \code{\link{plotBclimVol}}
}
\examples{
## See help for BclimMCMC and webpage at http://mathsci.ucd.ie/~parnell_a/Bclim.html
}
