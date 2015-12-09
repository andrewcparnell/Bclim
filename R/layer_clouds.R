layer_clouds = function(pollen, path_to_rs = system.file('data',package='Bclim'), nsamples=1000, G=10, mixwarnings=FALSE) {

# Function to turn pollen data into marginal data posteriors and then fit them to mixtures of normal distributions

# Load in the response surfaces from the package
load(paste(path_to_rs,'/data/required.data3D'))

############ PART 1 - CREATE MDPs - WRITTEN BY JAMES

# Number of slices
nslices = nrow(pollen.data)

# Format pollen data
#getting rowSums ~ 1000
myrowsums = rowSums(pollen.data)
myrowsums[myrowsums==0] = 1
Pollen = as.matrix(round(pollen.data/myrowsums*1000))
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
          as.double(quadprobs)
          ,PACKAGE="Bclim")
res = result[[2]]
res[res==0] = -Inf

##### return posterior dists in a matrix of dimension nslices x 2500
post = matrix(exp(res),nslices,175616,byrow=TRUE)
post=post/rowSums(post)

#####provide n samples from MTCO/GDD5 space including jittering
MDP = rep(0,nsamples*nslices*3)

for(i in 1:nslices) {
  Firstsamp = sample(1:175616, nsamples, replace = TRUE, prob = post[i,])
  jitter = cbind(runif(nsamples,-1,1),runif(nsamples,-1,1),runif(nsamples,-1,1))*cbind(rep(72.53,nsamples),rep(.7,nsamples),rep(10.206,nsamples))

  jit_locations = Climate[Firstsamp,]+jitter
  MDP[(i-1)*nsamples + (1:nsamples)] = jit_locations[,1]
  MDP[nslices*nsamples+(i-1)*nsamples + (1:nsamples)] = jit_locations[,2]
  MDP[nslices*nsamples*2+(i-1)*nsamples + (1:nsamples)] = jit_locations[,3]
}

#restructuring samples array
dim(MDP) = c(nsamples,nslices,3)

#naming dimensions
dimnames(MDP) = list(NULL,NULL,c("GDD5","MTCO","AET/PET"))

############ PART 2 - Approximate as mixtures - written by Andrew

# Calculate n.samp = number of samples, n = number of layers, m = number of climate dimensions
n.samp = dim(MDP)[1]
n = dim(MDP)[2]
m = 3

ScMean = rep(0,m)
ScVar = rep(1,m)
MDP2 = MDP
for(i in 1:m) {
  ScMean[i] = mean(MDP[,,i])
  ScVar[i] = median(diag(var(MDP[,,i])))
  MDP2[,,i] = (MDP[,,i]-ScMean[i])/sqrt(ScVar[i])
}

################# MIXTURE ESTIMATION #################

# Set up mixture components
mu.mat = array(NA,dim=c(n,m,G))
tau.mat = array(NA,dim=c(n,m,G))
p.mat = matrix(NA,nrow=n,ncol=G)

ans.all = list()
for(i in 1:n) {
  cat("\r")
  cat("Completed:",format(round(100*i/n,2), nsmall = 2),"%")
  ans.all[[i]] = mclust::Mclust(MDP2[,i,],G=G,modelNames="EII",warn=mixwarnings)
}

for(i in 1:n) {
  mu.mat[i,,] = ans.all[[i]]$parameters$mean
  for(g in 1:G) {
    tau.mat[i,,g] = 1/diag(ans.all[[i]]$parameters$variance$sigma[,,g])
  }
  p.mat[i,] = ans.all[[i]]$parameters$pro
  # End of loop through n layers
}


# Now output everything to a nice neat list
Bclimdata = list(MDP=MDP2,n=n,m=m,n.samp=n.samp,ScMean=ScMean,ScVar=ScVar,G=G,mu.mat=mu.mat,tau.mat=tau.mat,p.mat=p.mat)
class(Bclimdata) = 'layer_clouds'
cat("\n")
return(Bclimdata)

}
