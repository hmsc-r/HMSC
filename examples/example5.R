# Defining a HMSC with 1 latent factors levels (non-spatial at the level of ovservations),
# and traits. Defining X data through data.frame and formula.
rm(list=ls())
set.seed(1)

# download, install the package from GitHub and load it to the session
# library(devtools)
# install_github("gtikhonov/HMSC")
library(Hmsc)

ny = 201L
ns = 71L
nc1 = 2L
nc2 = 2L
nt0 = 1L
np = c(ny)

distr = "normal"

samples = 1000
thin = 1

fo = ~x1+poly(x2,degree=2,raw=TRUE)+c1+c2+x1*c2
TrFo = ~.
XData = data.frame(rnorm(ny))
for(i in 2:nc1)
   XData[,i] = rnorm(ny)
colnames(XData) = sprintf("x%d", 1:nc1)

for(i in 1:nc2){
   let = letters[1:sample(2:5,1)]
   XData[,sprintf("c%d", i)] = as.factor(let[(sample.int(ny)-1) %% length(let)+1])
}

X = model.matrix(fo,XData)
nc = ncol(X)
TrData = as.data.frame(matrix(rnorm(ns*nt0),ns,nt0))
TrMult = 100
TrData = TrMult*TrData
Tr = model.matrix(TrFo,TrData)

nt = ncol(Tr)


f0 = 1
V0 = diag(1,nc)
V = riwish(nc+1,V0)
mGamma0 = rep(0,nt*nc)
UGamma0 = diag(1,nt*nc)
Gamma = matrix(mvrnorm(1,mGamma0,UGamma0),nc,nt)
Gamma[,-1] = Gamma[,-1]/TrMult
Mu = tcrossprod(Gamma, Tr)
Beta = matrix(NA,nc,ns)
for(j in 1:ns)
   Beta[,j] = mvrnorm(1,Mu[,j],V)
# Beta[1,] = Beta[1,] + 1

nr = length(np)
dfPi = matrix(NA,ny,nr)
for(r in 1:nr){
   dfPi[,r] = (as.integer(((1:ny) - 1) %% np[r] + 1))
}
np = apply(dfPi,2,function(a) length(unique(a)))

sDim = c(2,0)
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

switch (distr,
   probit = {d = rep(1,ns)},
   poisson = {d = rep(0.01,ns)},
   normal = {d = rgamma(ns,2,1)}
)
LFix = X%*%Beta
Z = LFix + LRanSum + matrix(rnorm(ny*ns),ny,ns)*matrix(rep(sqrt(d),each=ny),ny,ns)

switch (distr,
   probit = {Y = 1*(Z>0)},
   poisson = {Y = matrix(rpois(ny*ns,exp(Z)), ny, ns)},
   normal = {Y = Z}
)

# Y[sample.int(length(Y),round(length(Y)/5))] = NA

BetaT = Beta
GammaT = Gamma
VT = V
sigmaT = d
LambdaT = Lambda
EtaT = Eta
AlphaT = Alpha

# create the main model and specify data, priors, parameters
m = Hmsc$new(Y=Y, XData=XData, XFormula=fo, XScale=TRUE, dist=distr, dfPi=dfPi, TrFormula=~., TrData=TrData, TrScale=TRUE, rL=rL)

start = proc.time()
m$sampleMcmc(samples, transient=round(samples/5), thin=thin, adaptNf=0*rep(2000,nr), updater=list(Gamma2=FALSE))
# , initPar=list(Alpha=AlphaT, Eta=EtaT, Lambda=LambdaT,
# Beta=BetaT, Gamma=GammaT, sigma=sigmaT))
stop = proc.time()

postList = m$postList[[1]]


XDataNew = XData[sample(1:10,ny,replace=TRUE),]
rLNew = vector("list", nr)
for(r in 1:nr){
   if(sDim[r]>0){
      sRLNew = rbind(rL[[r]]$s,rL[[r]]$s+0.1)
      rownames(sRLNew) = c(rownames(rL[[r]]$s), sprintf('new_%d',1:np[r]))
      rLNew[[r]] = HmscRandomLevel$new(data=sRLNew, priors="default")
   } else{
      rLNew[[r]] = HmscRandomLevel$new(N=np[r])
   }
}

dfPiNew = data.frame(sprintf("%s",dfPi[,1]))
p = m$predict(XData=XDataNew, dfPiNew=dfPiNew, rL=rLNew, predictEtaMean=TRUE)


postBeta = array(unlist(lapply(postList, function(a) a$Beta)),c(nc,ns,m$samples))
plot(Beta,apply(postBeta,c(1,2),mean),main="Beta")
abline(0,1,col="red")

postGamma = array(unlist(lapply(postList, function(a) a$Gamma)),c(nc,nt,m$samples))
plot(GammaT,apply(postGamma,c(1,2),mean),main="Gamma")
abline(0,1,col="red")

normalized = function(x){
   (x-min(x))/(max(x)-min(x))
}

Gradient = m$constructGradient(focalVariable = "x1")
predY = m$predict(XData=Gradient$XDataNew, dfPiNew=Gradient$dfPiNew, rL=Gradient$rLNew, expected=TRUE)
m$plotGradient(Gradient, pred=predY, measure="Y", index=1)

# VP = m$computeVariancePartitioning(c(1,1:(m$nc-1)),sprintf("%d",1:(m$nc-1)))
# m$plotVariancePartitioning(VP,ylim=c(0,10))
