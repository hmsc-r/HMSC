# Example with Ricardo's bats data
rm(list=ls())

# download, install the package from GitHub and load it to the session
# library(devtools)
# install_github("gtikhonov/MyFirstRPackage")
# library(MyFirstRPackage)


# read data from csv files
dataDir = "C:/Google Drive/HMSC/testing/Bats/data"
Xo = read.csv(file.path(dataDir,"Xo.csv"))
Y = read.csv(file.path(dataDir,"Y.csv"))
Tro = read.csv(file.path(dataDir,"To.csv"))
C = read.csv(file.path(dataDir,"C.csv"), header=FALSE)
S = read.csv(file.path(dataDir,"site.csv"), header=FALSE)

# transform X and Tr matrices to HMSC input format using built-in R function
X = model.matrix(formula("~habitat+time.period+effort+habitat*time.period"), data=Xo)
Tr = model.matrix(formula("~habitat.specialization+vertical.space.use+aspect_ratio+relative_wing_loading"), data=Tro)

# define the level-correspondence matrix
Pi = data.frame(L1=as.character(1:nrow(Y)),L2=as.character(S$V1))
colnames(Pi) = c("Sampling unit", "site")

# create 1 random levels and specify both data and priors for i
rL1 = HmscRandomLevel$new(pi=unique(Pi[,2]), priors="default")
rL1$setPriors(nfMax=10, mu=3, a1=5, b1=1, a2=3, b2=1)

# create the main model and specify data, priors, parameters
m = Hmsc$new(Y=Y, X=X, Tr=Tr, C=C, dist="probit", rL=list(rL1), Pi=Pi)
m$setPriors(priors="default")
m$setMcmcParameters()

m$sample()
m$getPosterior() # returns posterior

# postprocessing....
