### Hmsc-tests: unit tests for Hmsc functions

### BEGIN TESTS
suppressPackageStartupMessages(require(Hmsc))

### INITIALIZE DATA
set.seed(654)
n = 10
nPlot = 5
ns = 4
y = matrix(NA,ncol=ns,nrow=n)
y[,1] = rnorm(n,sd=1)
y[,2] = 1*((rnorm(n,sd=1))>0)
y[,3] = rpois(n=n,lambda=exp(rnorm(n,sd=1)))
y[,4] = rpois(n=n,lambda=exp(rnorm(n,sd=1)))

XData = data.frame(x1=rnorm(n),
                   x2=c('o','c','c','o','c','o','o','c','o','c'))
C = matrix(c(1,0.9718885,0.9313564,0,0.9718885,1,0.9313564,0,0.9313564,0.9313564,1,0,0,0,0,1),ncol=4,nrow=4)
Traits = data.frame(trait1 = rnorm(ns),
                    trait2 = c('a','b','a','b'))
xycoords = matrix(runif(2*nPlot),ncol=2)
rownames(xycoords) = 1:nPlot
studyDesign = data.frame(sample = as.factor(1:n),
                         plot = as.factor(sample(1:nPlot,n,replace=TRUE)))
XSelect = list()
for(k in 1:2){
   covGroup = k
   spGroup = c(1,1,2,2)
   q = rep(0.1,max(spGroup))
   XSelect[[k]] = list(covGroup=covGroup,spGroup=spGroup,qq=q)
}
rLxData = matrix(c(rep(1,n),XData$x1),ncol=2)
rownames(rLxData) = 1:n

load("~/Documents/Research/HMSC/BigSpatial/Scripts/HMSC/tests/refs/refs.Rdata")

### TEST HmscRandomLevel()
rL1 = HmscRandomLevel(xData = rLxData)
rL2 = HmscRandomLevel(sData = xycoords)
all.equal(rL1,rL1ref)
all.equal(rL2,rL2ref)

### TEST Hmsc()
m = Hmsc(Y=y,XData=XData,XFormula=~x1+x2, TrData=Traits,
         TrFormula = ~trait1+trait2, TrScale=TRUE, C=C, XSelect = XSelect,
         ranLevels=list("sample"=rL1ref,"plot"=rL2ref),
         studyDesign = studyDesign,
         distr=c('normal','probit','poisson','lognormal poisson'))
all.equal(m,mref)

### TEST sampleMCMC()
msampled = sampleMcmc(m,thin=1,samples=10,transient=2,nChains=1)
all.equal(msampled,msampledref)

### TEST computeWAIC()
WAIC = computeWAIC(msampledref)

