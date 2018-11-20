# Defining a HMSC without latent factors levels
rm(list=ls())
set.seed(1)

####################################### installation ########################################################
# library(devtools)
# install_github("gtikhonov/HMSC")
library(Hmsc)

ny = 1000
ns = 21
nc = 3
np2 = 70
distr = "normal"


####################################### generating data #####################################################
X = matrix(rnorm(ny*nc),ny,nc)
S = data.frame(s1=1:np2,s2=rnorm(np2))
Beta = matrix(rnorm(nc*ns),nc,ns)

dfPi = NULL
sigma = 1+0*rgamma(ns,1,1)

L = X %*% Beta
Y = L + matrix(rnorm(ny*ns),ny,ns)*matrix(sqrt(sigma),ny,ns,byrow=TRUE)
if(distr=="probit")
   Y = matrix(as.numeric(Y>0),ny,ns)

####################################### create the main model and specify data, priors, parameters ##########
m = Hmsc(Y=Y, X=X, dist=distr, dfPi)
m = sampleMcmc(m, 100, thin=10, adaptNf=c(200,200), nChains=2 )

####################################### postprocessing #####################################################
postBeta = array(unlist(lapply(poolMcmcChains(m$postList), function(a) a$Beta)),c(nc,ns,m$samples))
plot(Beta,apply(postBeta,c(1,2),mean),main="Beta")
abline(0,1,col="red")

codaList = convertToCodaObject(m)
gelman.diag(codaList$Beta)

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

# postSigma = array(unlist(lapply(m$postList, function(a) a$sigma)),c(ns,m$samples))
# plot(sigma,apply(postSigma,1,mean),main="Sigma")
# abline(0,1,col="red")
#
LFix = m$X%*%m$postList[[m$samples]]$Beta
LRan = vector("list", m$nr)
for(r in 1:m$nr){
   LRan[[r]] = m$postList[[100]]$Eta[[r]][m$Pi[,r],]%*%m$postList[[100]]$Lambda[[r]]
}
LAll = LFix + Reduce("+", LRan)

# plot(L,LAll)
# plot(X%*%Beta,LFix)


indProbit = rep(TRUE,ns)
lB = matrix(-Inf,ny,ns)
uB = matrix(Inf,ny,ns)
lB[as.logical(Y[,indProbit])] = 0
uB[!Y[,indProbit]] = 0

Z = matrix(rtruncnorm(ny*ns,a=lB,b=uB,mean=LAll),ny,ns)
# plot(LAll,Z,col=Y+1)

# j = 8
# plot(LAll[,j],Z[,j],col=Y[,j]+1)

# plot(L[,8],Y[,8])

# library(msm)
# Z = rtnorm(ny*ns, mean=LAll, sd=1, lower=lB,upper=uB)
#
# plot(LAll,Z)
