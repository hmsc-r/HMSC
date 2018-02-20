# Defining a HMSC with 2 latent factors levels (non-spatial at the level of ovservations and 2D at the level of sites),
# probit data, and traits
rm(list=ls())

# download, install the package from GitHub and load it to the session
# library(devtools)
# install_github("gtikhonov/HMSC")
library(Hmsc)

ny = 100
ns = 41
nc = 3
np2 = 33


# read or generate some data
X = matrix(rnorm(ny*nc),ny,nc)
S = data.frame(s1=1:np2,s2=rnorm(np2))
Beta = matrix(rnorm(nc*ns),nc,ns)
Y = X %*% Beta
Pi = data.frame(L1=as.character(1:ny), L2=as.character(((1:ny-1) %% np2)+1))



# create 2 random levels and specify both data and priors for them
rL1 = HmscRandomLevel$new(N=ny, priors="default")
rL2 = HmscRandomLevel$new(data=S)
rL2$setPriors(nfMax=10, mu=3, a1=5, b1=1, a2=3, b2=1, alphapw=cbind(p=c(0,1,2,3),w=c(0.4,0.2,0.2,0.2)))

# create the main model and specify data, priors, parameters
m = Hmsc$new(Y=Y, X=X, dist="normal", rL=list(rL1,rL2), Pi=Pi)
# m = Hmsc$new()
# m$setData(Y=Y)

# m$setPriors()

m
# m$setMcmcParameters()

m$sampleMcmc(100, thin=10)
# m$getPosterior() # returns posterior

# postprocessing....

library(coda)
mcmc.list()
mcmc

postBeta = array(unlist(lapply(m$postList, function(a) a$Beta)),c(nc,ns,m$samples))
plot(Beta,apply(postBeta,c(1,2),mean))
