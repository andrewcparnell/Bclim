\name{Bclim-package}
\alias{Bclim-package}
\alias{Bclim}
\docType{package}
\title{
Create palaeoclimate reconstructions from pollen cores
}
\description{
This package takes pollen counts and chronological data to produce climate histories. The package runs in two stages, first by calling \code{\link{slice_clouds}} and subsequently via \code{\link{climate_histories}}.
}
\details{
\tabular{ll}{
Package: \tab Bclim\cr
Type: \tab Package\cr
Version: \tab 3.0.0\cr
Date: \tab 2015-12-14\cr
License: \tab GPL >=2 \cr
}
See the functions \code{\link{slice_clouds}} and \code{\link{climate_histories}} for a full description of the workings of the package. See the package vignette for worked examples
}
\author{
Andrew Parnell, Thinh Doan and James Sweeney

Maintainer: Andrew Parnell <andrew.parnell@ucd.ie>
}
\references{
See the pacakge vignette (at https://cran.r-project.org/web/packages/Bclim/index.html), the package web page (at http://mathsci.ucd.ie/~parnell_a/Bclim.html), and the following references:
Salter-Townshend, M. and J. Haslett (2012). Fast Inversion of a Flexible Regression Model for Multivariate, Zero-Inflated Pollen Counts. Environmetrics.
Sweeney, J. (2012). Advances in Bayesian Model Development and Inversion in Multivariate Inverse Inference Problems with application to palaeoclimate reconstruction. Ph. D. thesis, Trinity College Dublin.
Parnell, A. C., et al. (2015), Bayesian inference for palaeoclimate with time uncertainty and stochastic volatility. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 115–138.
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
\keyword{ package }
\seealso{
\code{\link{slice_clouds}}, \code{\link{climate_histories}}
}
