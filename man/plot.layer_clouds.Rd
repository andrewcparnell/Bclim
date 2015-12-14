\name{plot.layer_clouds}
\alias{plot.layer_clouds}
\title{
Plots of Bclim layer clouds
}
\description{
Create bivariate climate plots of individual layers
}
\usage{
\method{plot}{layer_clouds}(x,
layer = 1,
dims = 1:2,
n = 50,...)
}
\arguments{
  \item{x}{
The output of a run from \code{\link{layer_clouds}}
}
  \item{layer}{
The chosen layer to plot
}
  \item{dims}{
A vector of length 2. dim=1 corresponds to GDD5, dim=2 to MTCO, and dim=3 to AET/PET
}
  \item{n}{
The resolution of the resulting plot. A higher value of n will yield finer plots but might require some colour adjustment
}
  \item{\dots}{
Other arguments to the plot function such as axis labels, titles, and colours
}
}
\details{
This function creates a bivariate density plot of two climate dimensions (two of GDD5, MTCO and AET/PET) using the MASS library function \code{\link{kde2d}}
}
\value{
No output, just a plot
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
The main Bclim functions are \code{\link{layer_clouds}} and \code{\link{climate_histories}}

For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the authors personal webpage (http://mathsci.ucd.ie/~parnell_a/Bclim.html)?
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
