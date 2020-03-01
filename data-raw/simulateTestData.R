#' Script used to generate test data 'TD.Rdata'
#' Last run: 2020-02-29

### This script works only in RStudio. If you cannot or do not work in
### RStudio alternatives: use the plain R alternative that is now
### commented out on the last line. To use the plain R alternative,
### comment out the usethis::use_data() line and remove comment dashes
### from the last line starting save(). Plain R alternative should be
### run in 'data-raw' working directory.

library(ape)
library(MASS)
library(Hmsc)
## you need the next only if you want to use RStudio. 
library(usethis)

#Set random seed for reproducibility
set.seed(66)

TD = list()
TD$ns = 4
TD$units = 50
TD$plots = 10

#Simulating one continuous environmental covariate
TD$X = cbind(rep(1,TD$units),rnorm(TD$units))

#Simulating species niches for 4 species, using one phylogenetically structured trait.  one continuous which is fully structured by phylogeny
TD$phy = rcoal(n=TD$ns, tip.label = sprintf('sp_%.3d',1:TD$ns), br = "coalescent")
TD$C = vcv(TD$phy, model = "Brownian", corr = TRUE)
TD$Tr = cbind(rep(1,TD$ns),mvrnorm(n = 1, mu = rep(0,TD$ns), Sigma=TD$C))

gamma = cbind(c(-2,2),cbind(c(-1,1)))
mu = gamma%*%t(TD$Tr)
beta = matrix(mvrnorm(n = 1, mu = as.vector(mu),
                      Sigma=kronecker(TD$C,diag(2))),ncol=TD$ns)
Lf = TD$X%*%beta

#Simulating linear predictor for spatial random effect
TD$studyDesign = data.frame(sample = as.factor(1:TD$units),
                            plot = as.factor(sample(1:TD$plots,TD$units,replace=TRUE)))
TD$xycoords = matrix(runif(2*TD$plots),ncol=2)
Sigma = 2^2*exp(-as.matrix(dist(TD$xycoords))/0.35)
etaPlot = mvrnorm(mu=rep(0,TD$plots), Sigma=Sigma)
etaUnit = etaPlot[TD$studyDesign$plot]
lambda = c(-2,2,1.5,0)
Lr = etaUnit%*%t(lambda)
rownames(TD$xycoords) = 1:TD$plots
colnames(TD$xycoords) = c("x-plots","y-coordinate")

#Simulating species occurences
TD$Y = as.matrix(Lf+Lr+matrix(rnorm(TD$units*TD$ns),ncol=TD$ns))
colnames(TD$Y) = colnames(TD$C)
TD$Y = 1*(TD$Y>0)

#Testmodel
rL1 = HmscRandomLevel(sData = TD$xycoords)
rL2 = HmscRandomLevel(units = TD$studyDesign$sample)
TD$rL1 = setPriors(rL1, nfMax=2, nfMin=2)
TD$rL2 = setPriors(rL2, nfMax=2, nfMin=2)
TD$Tr = data.frame(TD$Tr)
TD$Tr = cbind(TD$Tr, as.factor(c('A','B','B','A')))
colnames(TD$Tr) = c('Intercept','T1','T2')
TD$X = data.frame(TD$X)
TD$X = cbind(TD$X, as.factor(c(rep('o',TD$units/2),rep('c',TD$units/2))))
colnames(TD$X) = c('Intercept','x1','x2')

TD$m = Hmsc(Y=TD$Y,
         XData=TD$X,
         XFormula=~x1+x2,
         TrData=TD$Tr,
         TrFormula = ~T1 + T2,
         phyloTree=TD$phy,
         ranLevels=list("sample"=TD$rL2,"plot"=TD$rL1),
         studyDesign = TD$studyDesign,
         distr=c('probit'))

TD$m = sampleMcmc(TD$m,thin=1,samples=100,transient=50,nChains=2)
## The next works only within RStudio
usethis::use_data(TD, overwrite=TRUE)
## If you are not in RStudio use instead the following in the data-raw
## directory

## save(TD, file = file.path("..", "data", "TD.rda"), version = 2)
