\name{plot.climate_histories}
\alias{plot.climate_histories}
\title{
Plots of posterior Bclim climate histories
}
\description{
Create plots of climate histories from a Bclim run
}
\usage{
\method{plot}{climate_histories}(x,
dim=1,
layer_clouds=TRUE,
chron=NULL,
climate_ribbon=TRUE,
most_representative=1,
conf=c(0.95,0.75,0.5),
col_clouds = rgb(0,0,1,0.2),
col_ribbon=rgb(1,0,0,0.4),
col_representative=rgb(0,1,0),
present_left=TRUE,...)
}
\arguments{
  \item{x}{
The output of a Bclim run from \code{\link{climate_histories}}
}
  \item{dim}{
The chosen climate dimension. This could be GDD5 (dim=1), MTCO (dim=2) or AET/PET (dim=3)
}
  \item{layer_clouds}{
Whether to ploy the individual layer clouds (default TRUE)
}
  \item{chron}{
A chronology file (see \code{\link{climate_histories}} for details). Only required if layer_clouds=TRUE
}
  \item{climate_ribbon}{
Whether to plot the climate ribbon, i.e. the time slices quantiles (default=TRUE)
}
  \item{most_representative}{
The number of representative climate histories to plot. See Details section below. Can be set to zero if none are required
}
  \item{conf}{
The confidence levels of the layer clouds and the climate histories. Default is 95\%, 75\% and 50\% shading
}
  \item{col_clouds}{
The colour of the climate clouds. Default is blue with 20\% transparency
}
  \item{col_ribbon}{
The colour of the climate ribbon. Default is red with 40\% transparency
}
  \item{col_representative}{
The colour of the representative climate histories. Default is green
}
  \item{present_left}{
Whether the present (i.e. 0 years before present) should be on the left or the right of the plot. Default is to put it on the left
}
  \item{\dots}{
Other arguments to the plot function, such as axis labels, titles etc
}
}
\details{
This function creates the default Bclim plots of climate histories and layer clouds from a Bclim run. Users can turn on or off the layer clouds and summaries of the the climate histories (the `climate ribbon'), and change the confidence level shown on the plots. The function also allows for a number of `representative histories' to be plotted. These are considered to be the climate histories that are the median distance away from the point-wise medians.
}
\value{
No output, just a plot
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>
}
\seealso{
The main Bclim functions are \code{\link{layer_clouds}} and \code{\link{climate_histories}}.

For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the authors personal webpage (http://mathsci.ucd.ie/~parnell_a/Bclim.html)?
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }
