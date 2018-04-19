# Defining a HMSC with 2 latent factors levels (non-spatial at the level of ovservations and 2D at the level of sites),
# probit data, and traits
rm(list=ls())
set.seed(3)

# download, install the package from GitHub and load it to the session
# library(devtools)
# install_github("gtikhonov/HMSC")
library(Hmsc)

ny = 201L
ns = 41L
nc = 3L
nt = 1L
np = c(ny, round(ny/10))
# np = np[1]

distr = "normal"

samples = 1000
thin = 10

X = matrix(rnorm(ny*nc),ny,nc)
X[,1] = 1
Tr = matrix(rnorm(ns*nt),ns,nt)
Tr[,1] = 1

f0 = 1
V0 = diag(1,nc)
V = riwish(nc+1,V0)
mGamma0 = rep(0,nt*nc)
UGamma0 = diag(1,nt*nc)
Gamma = matrix(mvrnorm(1,mGamma0,UGamma0),nc,nt)
Mu = tcrossprod(Gamma, Tr)
Beta = matrix(NA,nc,ns)
for(j in 1:ns)
   Beta[,j] = mvrnorm(1,Mu[,j],V)

nr = length(np)
dfPi = matrix(NA,ny,nr)
for(r in 1:nr){
   dfPi[,r] = (as.integer(((1:ny) - 1) %% np[r] + 1))
}
np = apply(dfPi,2,function(a) length(unique(a)))

sDim = c(0,2)
nf = as.integer(c(2,2))

rL = vector("list", nr)
sRL = vector("list", nr)
Alpha = vector("list", nr)
for(r in 1:nr){
   if(sDim[r]>0){
      sRL[[r]] = matrix(runif(np[r]*sDim[r]),np[r],sDim[r])
      rownames(sRL[[r]]) = as.character(1:np[r]) #sprintf('%.3d',1:np[r])
      rL[[r]] = HmscRandomLevel$new(data=sRL[[r]], priors="default")
      Alpha[[r]] = c(80, 20)
   } else{
      rL[[r]] = HmscRandomLevel$new(N=np[r])
   }
}


for(r in 1:nr){
   dfPi[,r] = as.character(dfPi[,r]) #sprintf('%.3d',dfPi[,r])
}
dfPi = as.data.frame(dfPi)

Eta = vector("list", nr)
Lambda = vector("list", nr)
LRan = vector("list", nr)
for(r in 1:nr){
   if(rL[[r]]$sDim > 0){
      D = as.matrix(dist(rL[[r]]$s[levels(dfPi[,r]),]))
      Eta[[r]] = matrix(NA,np[r],nf[r])
      for(h in 1:nf[r]){
         if(rL[[r]]$alphapw[Alpha[[r]][h],1] > 0){
            K = exp(-D/rL[[r]]$alphapw[Alpha[[r]][h],1])
         } else
            K = diag(np[r])
         Eta[[r]][,h] = crossprod(chol(K), rnorm(np[r]))
      }
   } else{
      Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
   }
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

BetaT = Beta
GammaT = Gamma
VT = V
sigmaT = d
LambdaT = Lambda
EtaT = Eta
AlphaT = Alpha

# create the main model and specify data, priors, parameters
m = Hmsc$new(Y=Y, X=X, dist=distr, dfPi=dfPi, Tr=Tr, rL=rL)

start = proc.time()
m$sampleMcmc(samples, thin=thin, adaptNf=1*rep(2000,nr) )
# , initPar=list(Alpha=AlphaT, Eta=EtaT, Lambda=LambdaT,
#    Beta=BetaT, Gamma=GammaT, sigma=sigmaT)
stop = proc.time()

postList = m$postList

postBeta = array(unlist(lapply(postList, function(a) a$Beta)),c(nc,ns,m$samples))
plot(Beta,apply(postBeta,c(1,2),mean),main="Beta")
abline(0,1,col="red")

# postGamma = array(unlist(lapply(postList, function(a) a$Gamma)),c(nc,nt,m$samples))
# plot(GammaT,apply(postGamma,c(1,2),mean),main="Gamma")
# abline(0,1,col="red")


normalized = function(x){
   (x-min(x))/(max(x)-min(x))
}

getOmega = function(a,r=1)
   return(crossprod(a$Lambda[[r]]))
getEta = function(a,r=1)
   return(a$Eta[[r]])
getAlpha = function(a,r=1)
   return(a$Alpha[[r]])
for(r in 1:m$nr){
   postOmega = array(unlist(lapply(postList,getOmega,r)),c(ns,ns,m$samples))
   postEta = array(unlist(lapply(postList,getEta,r)),c(np[r],nf[r],m$samples))
   plot(apply(postEta^2, 3, mean))
   postOmegaMean = apply(postOmega,c(1,2),mean)
   plot(crossprod(LambdaT[[r]]),postOmegaMean,main=sprintf("Omega%d",r))
   abline(0,1,col="red")
   if(m$rL[[r]]$sDim > 0){
      postAlpha = array(unlist(lapply(postList,getAlpha,r)),c(nf[r],m$samples))
      boxplot(t(postAlpha))
      points(1:nf[r], AlphaT[[r]], cex=3, col="red", pch=19)
      for(h in 1:nf[r]){
         par(mfrow=c(1,2))
         cols = colorRamp(c("blue","white","red"))(normalized(Eta[[r]][,h]))
         cols = apply(cols, 1, function(c) rgb(c[1]/255, c[2]/255, c[3]/255))
         plot(rL[[r]]$s[levels(dfPi[,r]),1], rL[[r]]$s[levels(dfPi[,r]),2], col=cols, pch = 19, cex=2, main=sprintf("True, alpha=%.3f",rL[[r]]$alphapw[AlphaT[[r]][h],1]),
            xlab="x", ylab="y")
         cols = colorRamp(c("blue","white","red"))(normalized(apply(postEta,c(1,2),mean)[,h]))
         cols = apply(cols, 1, function(c) rgb(c[1]/255, c[2]/255, c[3]/255))
         plot(rL[[r]]$s[levels(dfPi[,r]),1], rL[[r]]$s[levels(dfPi[,r]),2], col=cols, pch = 19, cex=2, xlab="x", ylab="y",
            main=sprintf("Estimated, alpha=%.3f",mean(rL[[r]]$alphapw[postAlpha[h,],1])) )
         par(mfrow=c(1,1))
      }
   }
   plot(apply(postEta, 3, function(c) sum((c-EtaT[[r]])^2)))
}




