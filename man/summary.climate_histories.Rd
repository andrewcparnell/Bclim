\name{summary.climate_histories}
\alias{summary.climate_histories}
\title{
Summarises the output created by \code{\link{climate_histories}}
}
\description{
Produces estimated climate values for a chosen climate dimension for each of the values supplied to the \code{time_grid} argument to \code{\link{climate_histories}}
}
\usage{
\method{summary}{climate_histories}(object,
dim=1,
probs=c(0.05,0.95),...)
}
\arguments{
  \item{object}{
An object of class \code{climate_histories} produced by the function \code{\link{climate_histories}}
}
  \item{dim}{
The chosen climate dimension. This could be GDD5 (dim=1), MTCO (dim=2) or AET/PET (dim=3)
}
  \item{probs}{
The chosen values at which to compute time-wise quantiles. The default is a 90\% interval, i.e. from 5\% to 95\%
}
  \item{...}{
Not used
}

}
\details{
The output is a table of time-wise confidence/credibility intervals for the climate histories at each time point given on the time grid for the specified climate dimension. The results can be saved in an object if required.
}
\value{
A data frame with the following columns:
\item{time_grid }{The provided time grid points}
\item{quantiles }{The quantiles of the climate variable for the specified probabilities}
Note that this object is reported silently so will be discarded unless the function is called with an object as in the vignette.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
See \code{\link{climate_histories}} for creating objects suitable for this function

For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the authors personal webpage (http://mathsci.ucd.ie/~parnell_a/Bclim.html)?
}
