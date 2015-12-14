\name{layer_clouds}
\alias{layer_clouds}
\title{
Function to approximate pollen layers as climate clouds
}
\description{
This function takes a set of pollen data and turns it layer-by-layer into climate estimates
}
\usage{
layer_clouds(pollen,
path_to_rs = 'http://mathsci.ucd.ie/~parnell_a/',
n_samples=1000)
}
\arguments{
  \item{pollen}{
A matrix or data frame of pollen counts (they can be normalised or not) which contains an unspecified number of rows and precisely 28 columns. These columns should represent counts of the following taxa in order: Abies Alnus Betula Carpinus Castanea Cedrus Corylus Ephedra Fagus Juniperus Larix Olea Ostrya Phillyrea Picea Pinus.D Pinus.H Pistacia Quercus.D Quercus.E Salix Tilia Ulmus Artemisia Chenopodiaceae Cyperaceae Ericales Gramineae
}
  \item{path_to_rs}{
A web address which links to the file \code{requireddata3D.RData} which contains response surfaces. The default should work fine
}
  \item{n_samples}{
The number of samples taken for each layer cloud. Default is 1000
}
}
\details{
A layer cloud is a multivariate probability distribution of the three climate dimensions (Growing Degree Days above 5C, GDD5; Mean Temperature of Coldest Month, MTCO; the ratio of actual to potential evapotranspiration, AET/PET) given the pollen information at that layer only. This function loops through each layer in the core to produce layer clouds which represent the information about climate obtained only from that layer of pollen. See references below for the technical details of this technique
}
\value{
A list object the the following elements
  \item{layer_clouds}{The layer clouds, an n_samples x n_layers x n_dimensions array}
  \item{n_samples}{The number of layers (i.e. the number of rows in the pollen file)}
  \item{n_dimensions}{The number of climate dimensions (currently always 3)}
}
\references{
Fore more detail on the algorithm see:
Salter-Townshend, M. and J. Haslett (2012). Fast Inversion of a Flexible Regression Model for Multivariate, Zero-Inflated Pollen Counts. Environmetrics.
Sweeney, J. (2012). Advances in Bayesian Model Development and Inversion in Multivariate Inverse Inference Problems with application to palaeoclimate reconstruction. Ph. D. thesis, Trinity College Dublin.
Parnell, A. C., et al. (2015), Bayesian inference for palaeoclimate with time uncertainty and stochastic volatility. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 115–138.
}
\author{
Andrew Parnell <andrew.parnell@ucd.ie>, James Sweeney <James.Sweeney@ucd.ie>, and Thinh Doan
}
\seealso{
\code{\link{climate_histories}}, \code{\link{plot.layer_clouds}}

For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the authors personal webpage (http://mathsci.ucd.ie/~parnell_a/Bclim.html)?
}
\keyword{ model }
\keyword{ multivariate }
\keyword{ smooth }

