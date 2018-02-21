# Defining a HMSC with 2 latent factors levels (non-spatial at the level of ovservations and 2D at the level of sites),
# probit data, and traits
rm(list=ls())
set.seed(1)

# download, install the package from GitHub and load it to the session
# library(devtools)
# install_github("gtikhonov/HMSC")
library(Hmsc)
library(lattice)

ny = 1000
ns = 41
nc = 3
np2 = 70


# read or generate some data
X = matrix(rnorm(ny*nc),ny,nc)
S = data.frame(s1=1:np2,s2=rnorm(np2))
Beta = matrix(rnorm(nc*ns),nc,ns)

Pi = data.frame(L1=as.character(1:ny), L2=as.character(((1:ny-1) %% np2)+1))
nf = c(3,2)
np = apply(Pi,2,function(a) return(length(unique(a))))
Lambda1 = matrix(rnorm(nf[1]*ns),nf[1],ns)
Lambda2 = matrix(rnorm(nf[2]*ns),nf[2],ns)
Eta1 = matrix(rnorm(np[1]*nf[1]),np[1],nf[1])
# Eta1[,2] = 0
Eta2 = matrix(rnorm(np[2]*nf[2]),np[2],nf[2])
sigma = 1+0*rgamma(ns,1,1)

L = X %*% Beta + Eta1[as.numeric(as.character(Pi$L1)),]%*%Lambda1 + Eta2[as.numeric(as.character(Pi$L2)),]%*%Lambda2
Y = L + matrix(rnorm(ny*ns),ny,ns)*matrix(sqrt(sigma),ny,ns,byrow=TRUE)
# Y = matrix(as.numeric(Y>0),ny,ns)


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

m$sampleMcmc(100, thin=10, adaptNf=c(100,100)) #initPar=list(Eta=list(Eta1,Eta2))
# m$getPosterior() # returns posterior

# postprocessing....

library(coda)


postBeta = array(unlist(lapply(m$postList, function(a) a$Beta)),c(nc,ns,m$samples))
plot(Beta,apply(postBeta,c(1,2),mean),main="Beta")
abline(0,1,col="red")

# mcmcBeta = as.mcmc(matrix(as.vector(postBeta),m$samples,nc*ns,byrow=TRUE))
# plot(mcmcBeta)
# acfplot(mcmcBeta)

getOmega = function(a,r=1)
   return(crossprod(a$Lambda[[r]]))
postOmega1 = array(unlist(lapply(m$postList,getOmega)),c(ns,ns,m$samples))
postOmega1Mean = apply(postOmega1,c(1,2),mean)
plot(crossprod(Lambda1),postOmega1Mean,main="Omega1")
abline(0,1,col="red")

postOmega2 = array(unlist(lapply(m$postList,getOmega,2)),c(ns,ns,m$samples))
postOmega2Mean = apply(postOmega2,c(1,2),mean)
plot(crossprod(Lambda2),postOmega2Mean,main="Omega2")
abline(0,1,col="red")

# levelplot(postOmega1Mean)
# levelplot(crossprod(Lambda1))

postSigma = array(unlist(lapply(m$postList, function(a) a$sigma)),c(ns,m$samples))
plot(sigma,apply(postSigma,1,mean),,main="Sigma")
abline(0,1,col="red")
