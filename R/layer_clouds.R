layer_clouds = function(pollen, path_to_rs = 'http://mathsci.ucd.ie/~parnell_a/', n_samples=1000) {

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
  jitter = cbind(runif(n_samples,-1,1),runif(n_samples,-1,1),runif(n_samples,-1,1))*cbind(rep(72.53,n_samples),rep(.7,n_samples),rep(10.206,n_samples))

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

Bclimdata = list(layer_clouds=MDP,n_samples=n_samples,n_layers=nslices,n_dimensions=3)
class(Bclimdata) = 'layer_clouds'
return(Bclimdata)

}
