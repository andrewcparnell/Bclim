#' Function to approximate pollen slices as climate clouds
#'
#' This function takes a set of pollen data and turns it slice-by-slice into climate estimates. For examples why not see the wonderful Bclim vignette (available at https://cran.r-project.org/web/packages/Bclim/index.html) and the author's personal webpage (https://maths.ucd.ie/parnell)?

#'
#' @param pollen A matrix or data frame of pollen counts (they can be normalised or not) which contains an unspecified number of rows and precisely 28 columns. These columns should represent counts of the following taxa in order: Abies Alnus Betula Carpinus Castanea Cedrus Corylus Ephedra Fagus Juniperus Larix Olea Ostrya Phillyrea Picea Pinus.D Pinus.H Pistacia Quercus.D Quercus.E Salix Tilia Ulmus Artemisia Chenopodiaceae Cyperaceae Ericales Gramineae
#' @param path_to_rs A web address which links to the file \code{requireddata3D.RData} which contains response surfaces. The default should work fine
#' @param n_samples The number of samples taken for each slice cloud. Default is 1000
#'
#' @details A slice cloud is a multivariate probability distribution of the three climate dimensions (Growing Degree Days above 5C, GDD5; Mean Temperature of Coldest Month, MTCO; the ratio of actual to potential evapotranspiration, AET/PET) given the pollen information at that slice only. This function loops through each slice in the core to produce slice clouds which represent the information about climate obtained only from that slice of pollen. See references below for the technical details of this technique
#'
#' @useDynLib PalaeoRecon3D
#'
#' @references Fore more detail on the algorithm see:
#' Salter-Townshend, M. and J. Haslett (2012). Fast Inversion of a Flexible Regression Model for Multivariate, Zero-Inflated Pollen Counts. Environmetrics.
#' Sweeney, J. (2012). Advances in Bayesian Model Development and Inversion in Multivariate Inverse Inference Problems with application to palaeoclimate reconstruction. Ph. D. thesis, Trinity College Dublin.
#' Parnell, A. C., et al. (2015), Bayesian inference for palaeoclimate with time uncertainty and stochastic volatility. Journal of the Royal Statistical Society: Series C (Applied Statistics), 64: 115–138.
#' 
#' @return A list object the the following elements
#' \itemize{
#'   \item{slice_clouds}{The slice clouds, an n_samples x n_slices x n_dimensions array}
#' \item{n_samples}{The number of slices (i.e. the number of rows in the pollen file)}
#' \item{n_dimensions}{The number of climate dimensions (currently always 3)
#'   }
#' }
#' @export
#'
#' @seealso \code{\link{climate_histories}}, \code{\link{plot.slice_clouds}}
slice_clouds = function(pollen, path_to_rs = 'https://maths.ucd.ie/parnell/', n_samples=1000) {

# Function to turn pollen data into marginal data posteriors and then fit them to mixtures of normal distributions

# Load in the response surfaces from the package
required.data3D = NULL # Fix so that R doesn't complain about 'visible bindings'
con = url(paste0(path_to_rs,'requireddata3D.RData'))
load(con)

############ PART 1 - CREATE MDPs - WRITTEN BY JAMES

# Number of slices
nslices = nrow(pollen)

# Format pollen data
#getting rowSums ~ 1000
myrowsums = rowSums(pollen)
myrowsums[myrowsums==0] = 1
Pollen = as.matrix(round(pollen/myrowsums*1000))
Mu = V = matrix(0,175616,28)

##### loading predictors & parameters
M = required.data3D$Mu
Var = required.data3D$V
quadpts = required.data3D$quadpts
Id = required.data3D$Id
Mu[Id,] = M
V[Id,] = Var
Mu = c(Mu)
V = c(V)
quadprobs = required.data3D$quadprobs
nquadpts = 13
alpha = required.data3D$alpha
delta = required.data3D$delta
Buffer = required.data3D$Buffer
Resultslices = rep(0,length=(175616*nslices))
Climate = required.data3D$Climate

##### Running C function
result = .C("PalaeoRecon3D",
          as.integer(nslices),
          as.double(Resultslices),
          as.integer(Pollen),
          as.double(Mu),
          as.double(V),
          as.double(alpha),
          as.double(delta),
          as.double(Buffer),
          as.integer(nquadpts),
          as.double(quadpts),
          as.double(quadprobs))
res = result[[2]]
res[res==0] = -Inf

##### return posterior dists in a matrix of dimension nslices x 2500
post = matrix(exp(res),nslices,175616,byrow=TRUE)
post=post/rowSums(post)

#####provide n samples from MTCO/GDD5 space including jittering
MDP = rep(0,n_samples*nslices*3)

for(i in 1:nslices) {
  Firstsamp = sample(1:175616, n_samples, replace = TRUE, prob = post[i,])
  jitter = cbind(stats::runif(n_samples,-1,1),stats::runif(n_samples,-1,1),stats::runif(n_samples,-1,1))*cbind(rep(72.53,n_samples),rep(.7,n_samples),rep(10.206,n_samples))

  jit_locations = Climate[Firstsamp,]+jitter
  MDP[(i-1)*n_samples + (1:n_samples)] = jit_locations[,1]
  MDP[nslices*n_samples+(i-1)*n_samples + (1:n_samples)] = jit_locations[,2]
  MDP[nslices*n_samples*2+(i-1)*n_samples + (1:n_samples)] = jit_locations[,3]
}

#restructuring samples array
dim(MDP) = c(n_samples,nslices,3)

#naming dimensions
dimnames(MDP) = list(NULL,NULL,c("GDD5","MTCO","AET/PET"))

# Divide AET/PET by 1000 so that it's a proportion again
MDP[,,3] = MDP[,,3]/1000

Bclimdata = list(slice_clouds=MDP,n_samples=n_samples,n_slices=nslices,n_dimensions=3)
class(Bclimdata) = 'slice_clouds'
return(Bclimdata)

}
