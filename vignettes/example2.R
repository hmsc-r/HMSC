# Defining a HMSC with 2 latent factors levels (non-spatial at the level of ovservations and 2D at the level of sites),
# normal data, traits and phylogeny
rm(list=ls())
set.seed(1)

# download, install the package from GitHub and load it to the session
# library(devtools)
# install_github("gtikhonov/HMSC")
library(Hmsc)

ny = 1001L
ns = 101L
nc = 13L
nt = 1L
rho = 0.2


distr = "normal"

samples = 100
thin = 10

X = matrix(rnorm(ny*nc),ny,nc)
X[,1] = 1
Tr = matrix(rnorm(ns*nt),ns,nt)
Tr[,1] = 1
C = 1-as.matrix(dist((1:ns)/ns))
Q = rho*C + (1-rho)*diag(ns)

f0 = 1
V0 = diag(1,nc)
# V = riwish(nc+1,V0)
V = diag(nc)
mGamma0 = rep(0,nt*nc)
UGamma0 = diag(1,nt*nc)
Gamma = matrix(mvrnorm(1,mGamma0,UGamma0),nc,nt)
Mu = tcrossprod(Gamma, Tr)
Beta = matrix(NA,nc,ns)
Beta = Mu + matrix(mvrnorm(1, rep(0,nc*ns), kronecker(Q,V)), nc,ns)
# image(Q)
# image(cor(Beta))
# aaa

np = as.integer(c(ny, round(ny/10)))
nr = 1
dfPi = matrix(NA,ny,nr)
for(r in 1:nr){
   dfPi[,r] = (as.integer(((1:ny) - 1) %% np[r] + 1))
}

nf = as.integer(c(2,2))
np = apply(dfPi,2,function(a) length(unique(a)))
Eta = vector("list", nr)
Lambda = vector("list", nr)
LRan = vector("list", nr)
for(r in 1:nr){
   Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
   Lambda[[r]] = matrix(rnorm(nf[r]*ns),nf[r],ns)
   LRan[[r]] = Eta[[r]][dfPi[,r],]%*%Lambda[[r]]
}
LRanSum = Reduce("+", LRan)

# d = rep(1,ns)
d = rgamma(ns,2,1)
LFix = X%*%Beta
Y = LFix + LRanSum + matrix(rnorm(ny*ns),ny,ns)*matrix(rep(sqrt(d),each=ny),ny,ns)
if(distr == "probit")
   Y = 1*(Y>0)

for(r in 1:nr){
   dfPi[,r] = as.character(dfPi[,r])
}
BetaT = Beta
GammaT = Gamma
VT = V
sigmaT = d
LambdaT = Lambda
EtaT = Eta

# create the main model and specify data, priors, parameters
m = Hmsc$new(Y=Y, X=X, dist=distr, Pi=dfPi, Tr=Tr, C=C)

initPar = list(Beta=BetaT, Gamma=GammaT, sigma=sigmaT, Eta=EtaT, Lambda=LambdaT, V=VT)
start = proc.time()
m$sampleMcmc(samples, thin=thin, adaptNf=0*c(200,200), initPar=initPar)
stop = proc.time()

# postprocessing....

postBeta = array(unlist(lapply(m$postList, function(a) a$Beta)),c(nc,ns,m$samples))
plot(Beta,apply(postBeta,c(1,2),mean),main="Beta")
abline(0,1,col="red")

postGamma = array(unlist(lapply(m$postList, function(a) a$Gamma)),c(nc,nt,m$samples))
plot(GammaT,apply(postGamma,c(1,2),mean),main="Gamma")
abline(0,1,col="red")


# mcmcBeta = as.mcmc(matrix(as.vector(postBeta),m$samples,nc*ns,byrow=TRUE))
# plot(mcmcBeta)
# acfplot(mcmcBeta)

getOmega = function(a,r=1)
   return(crossprod(a$Lambda[[r]]))
for(r in 1:m$nr){
   postOmega = array(unlist(lapply(m$postList,getOmega,r)),c(ns,ns,m$samples))
   postOmegaMean = apply(postOmega,c(1,2),mean)
   plot(crossprod(Lambda[[r]]),postOmegaMean,main=sprintf("Omega%d",r))
   abline(0,1,col="red")
}
print(stop-start)
aaa


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
