\name{plotBclimVol}
\alias{plotBclimVol}
\title{
Plotting of interpolated climate volatilities
}
\description{
This function plots volatilities over time for the selected climate dimension after a Bclim run.
}
\usage{
plotBclimVol(x, dim = 1, title = NULL, presentleft = TRUE, denscol = "red", 
  denstransp = 0.5, leg = TRUE, mean = TRUE, legloc = "topleft", ...)
}
\arguments{
  \item{x}{
A Bclim output object either created by \code{\link{BclimCompile}} or \code{\link{BclimRun}}.
}
  \item{dim}{
The chosen climate dimension. This could be GDD5 (dim=1), MTCO (dim=2) or AET/PET x 100 (dim=3)
}
  \item{title}{
The title for the plot
}
  \item{presentleft}{
Whether you wish to have the present on the left hand side (default) or on the right hand side
}
  \item{denscol}{
The colour of the 95\% volatility ranges
}
  \item{denstransp}{
The transparency (between 0 and 1) of the 95\% volatility ranges
}
  \item{leg}{
Whether to include a legend or not (default TRUE)
}
  \item{mean}{
Whether to include a line indicating the mean volatility
}
  \item{legloc}{
Where to place the legend (e.g. 'topleft', 'bottomright', etc)
}
  \item{\dots}{
Other arguments to the plot function
}
}
\value{
No output produced, just a plot.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie> and James Sweeney
}
\seealso{
The main Bclim function is \code{\link{BclimRun}}. See there for all the other information.
}
\examples{
## See main function for examples
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
