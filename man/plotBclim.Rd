\name{plotBclim}
\alias{plotBclim}
\title{
Plots of posterior Bclim climate histories 
}
\description{
Create plots of climate histories from a Bclim run
}
\usage{
plotBclim(x, dim = 1, title = NULL, presentleft = TRUE, blob = TRUE, MDPcol = "blue", 
  denscol = "red", MDPtransp = 0.1, denstransp = 0.5, leg = TRUE, legloc = "topleft", ...)
}
\arguments{
  \item{x}{
The output of a Bclim run from either \code{\link{BclimCompile}} or \code{\link{BclimRun}}.
}
  \item{dim}{
The chosen climate dimension. This could be GDD5 (dim=1), MTCO (dim=2) or AET/PET x 100 (dim=3)
}
  \item{title}{
The title for the plot
}
  \item{presentleft}{
ether you wish to have the present on the left hand side (default) or on the right hand side
}
  \item{blob}{
Whether to include the marginal data posteriors (MDPs) as blobs or not
}
  \item{MDPcol}{
The colour of the MDPs (only works if blob=TRUE)
}
  \item{denscol}{
The colour of the 50/75/95\% posterior climate ranges
}
  \item{MDPtransp}{
The transparency (between 0 and 1) of the MDPs
}
  \item{denstransp}{
The transparency (between 0 and 1) of the posterior climate ranges
}
  \item{leg}{
Whether to include a legend or not (default TRUE)
}
  \item{legloc}{
Where to place the legend (e.g. 'topleft', 'bottomright', etc)
}
  \item{\dots}{
Other arguments to the plot function
}
}
\details{
This function creates the default Bclim plots of climate posteriors and MDPs from a Bclim run. If you don't like these then remember that the Bclim object created for x above contains all the information you need to create any other plot you would like.
}
\value{
No output, just a plot
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie> and James Sweeney
}
\seealso{
The main Bclim function is \code{\link{BclimRun}}. See there for all the other information.
}
\examples{
## See help for BclimMCMC and webpage at http://mathsci.ucd.ie/~parnell_a/Bclim.html
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
