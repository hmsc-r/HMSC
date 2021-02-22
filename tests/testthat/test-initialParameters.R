context('Test computeInitialParameters with standard settings')

test_that("computeInitialParameters gives right number of output parameters",{
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList),12)
})

test_that("Standard initial beta is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Beta),12)
   expect_equal(round(parList$Beta),t(matrix(c(0,0,1,0,0,-1,1,-1,-2,-1,0,-1),nrow=4,ncol=3)))
})

test_that("Standard initial gamma is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Gamma),9)
   expect_equal(round(parList$Gamma),t(matrix(c(0,0,0,0,0,0,-1,1,0),nrow=3,ncol=3)))
})

test_that("Standard initial V is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$V),9)
   expect_equal(round(parList$V),t(matrix(c(0,0,0,0,1,0,0,0,0),nrow=3,ncol=3)))
})

test_that("Standard initial Sigma is correct",{
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(parList$sigma,rep(1,4))
})

test_that("Standard initial Eta is distributed as N(0,1)",{
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Eta),2)
   expect_equal(length(parList$Eta[[1]]),100)
   expect_equal(length(parList$Eta[[2]]),20)
   expect_equal(round(sd(parList$Eta[[1]])),1)
   expect_equal(round(sd(parList$Eta[[2]])),1)
   expect_equal(round(mean(parList$Eta[[1]])),0)
   expect_equal(round(mean(parList$Eta[[2]])),0)
})

test_that("Standard initial Alpha is correct",{
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Alpha),2)
   expect_equal(length(parList$Alpha[[1]]),2)
   expect_equal(length(parList$Alpha[[2]]),2)
   expect_equal(parList$Alpha[[1]],c(1,1))
   expect_equal(parList$Alpha[[2]],c(1,1))
})

test_that("Standard initial Rho is correct",{
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(parList$rho,1)
})

test_that("Standard initial Lambda is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Lambda),2)
   expect_equal(length(parList$Lambda[[1]]),8)
   expect_equal(length(parList$Lambda[[2]]),8)
   expect_equal(round(parList$Lambda[[1]]),matrix(rep(0,8),nrow=2,ncol=4))
   expect_equal(round(parList$Lambda[[2]]),matrix(rep(0,8),nrow=2,ncol=4))
})

test_that("Standard initial Psi is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Psi),2)
   expect_equal(length(parList$Psi[[1]]),8)
   expect_equal(length(parList$Psi[[2]]),8)
   expect_equal(round(parList$Psi[[1]]),t(matrix(c(0,1,0,1,1,0,0,1),nrow=4,ncol=2)))
   expect_equal(round(parList$Psi[[2]]),t(matrix(c(1,2,0,3,1,2,1,0),nrow=4,ncol=2)))
})

test_that("Standard initial Delta is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(length(parList$Delta),2)
   expect_equal(length(parList$Delta[[1]]),2)
   expect_equal(length(parList$Delta[[2]]),2)
   expect_equivalent(round(parList$Delta[[1]]),c(42,51))
   expect_equivalent(round(parList$Delta[[2]]),c(46,53))
})

test_that("Standard initial Z is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   expect_equal(nrow(parList$Z),TD$units)
   expect_equal(ncol(parList$Z),TD$ns)
   expect_equal(round(colMeans(parList$Z)),c(0,-1,0,-1))
})


## -----------------------------------------------------------------------------
context('Test computeInitialParameters based on intial parameters from fixed effects models')

test_that("Fixed effects initial beta is correct",{
   set.seed(200)
   Y=matrix(1:20,nrow=5,ncol=4)
   Y[,2]= sample(c(0,1),5,replace = TRUE)
   m = Hmsc(Y=Y,XData=data.frame(x1 = 1:5),distr = c('normal','probit','poisson','lognormal poisson'))
   parList = computeInitialParameters(m,initPar='fixed effects')
   fit1 = lm.fit(m$XScaled,m$Y[,1])
   fit2 = glm.fit(m$XScaled,m$Y[,2],family=binomial(link="probit"))
   fit3 = glm.fit(m$XScaled,m$Y[,3],family=poisson())
   fit4 = glm.fit(m$XScaled,m$Y[,4],family=poisson())
   expect_equivalent(parList$Beta[,1],fit1$coefficients)
   expect_equivalent(parList$Beta[,2],fit2$coefficients)
   expect_equivalent(parList$Beta[,3],fit3$coefficients)
   expect_equivalent(parList$Beta[,4],fit4$coefficients)
})

test_that("Fixed effects initial gamma",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar='fixed effects')
   fit1 = lm.fit(TD$m$Tr,parList$Beta[1,])
   fit2 = lm.fit(TD$m$Tr,parList$Beta[2,])
   fit3 = lm.fit(TD$m$Tr,parList$Beta[3,])
   expect_equivalent(parList$Gamma[1,],fit1$coefficients)
   expect_equivalent(parList$Gamma[2,],fit2$coefficients)
   expect_equivalent(parList$Gamma[3,],fit3$coefficients)
})

test_that("Fixed effects initial V",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar='fixed effects')
   expect_equal(round(parList$Gamma), matrix(c(-8,7,2,-17,14,5,1,-1,-12),nrow=3,ncol=3))
})

## -----------------------------------------------------------------------------
context('Test computeDataParameters')

test_that("computeDataParameters gives right number of parameters",{
   parList = computeDataParameters(TD$m)
   expect_equal(length(parList),5)
   expect_equal(length(parList$rLPar[[1]]),0)
   expect_equal(length(parList$rLPar[[2]]),4)
})

test_that("phylogenetic parameters when phylogeny is provided",{
   set.seed(200)
   parList = computeDataParameters(TD$m)
   expect_equal(length(parList$detQg),nrow(TD$m$rhopw))
   expect_equal(round(sum(parList$detQg)),-74)
   expect_equal(dim(parList$Qg),c(TD$m$ns,TD$m$ns,nrow(TD$m$rhopw)))
   expect_equal(dim(parList$iQg),c(TD$m$ns,TD$m$ns,nrow(TD$m$rhopw)))
   expect_equal(dim(parList$RQg),c(TD$m$ns,TD$m$ns,nrow(TD$m$rhopw)))
   expect_equal(round(sum(parList$Qg)),655)
   expect_equal(round(sum(parList$iQg)),280)
   expect_equal(round(sum(parList$RQg)),478)
})

test_that("phylogenetic parameters when phylogeny is not provided",{
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:20,nrow=10,ncol=2))
   parList = computeDataParameters(m)
   expect_equal(parList$detQg,0)
   expect_equal(dim(parList$Qg),c(m$ns,m$ns,1))
   expect_equal(dim(parList$iQg),c(m$ns,m$ns,1))
   expect_equal(dim(parList$RQg),c(m$ns,m$ns,1))
   expect_equal(parList$Qg[,,1],matrix(c(1,0,0,1),nrow=2,ncol=2))
   expect_equal(parList$iQg[,,1],matrix(c(1,0,0,1),nrow=2,ncol=2))
   expect_equal(parList$RQg[,,1],matrix(c(1,0,0,1),nrow=2,ncol=2))
})

test_that("spatial parameters when spatial data is not provided",{
   m = Hmsc(Y=matrix(1:20,nrow=10,ncol=2),X=matrix(1:20,nrow=10,ncol=2))
   parList = computeDataParameters(m)
   expect_equal(length(parList$rLPar),0)
})

test_that("spatial parameters when spatial data is provided",{
   set.seed(200)
   parList = computeDataParameters(TD$m)
   expect_equal(length(parList$detQg),nrow(TD$m$rhopw))
   expect_equal(round(sum(parList$rLPar[[2]]$detWg)),-659)
   expect_equivalent(dim(parList$rLPar[[2]]$Wg),c(TD$plots,TD$plots,nrow(TD$m$rL[[2]]$alphapw)))
   expect_equivalent(dim(parList$rLPar[[2]]$iWg),c(TD$plots,TD$plots,nrow(TD$m$rL[[2]]$alphapw)))
   expect_equivalent(dim(parList$rLPar[[2]]$RiWg),c(TD$plots,TD$plots,nrow(TD$m$rL[[2]]$alphapw)))
   expect_equal(round(sum(parList$rLPar[[2]]$Wg)),4875)
   expect_equal(round(sum(parList$rLPar[[2]]$iWg)),314)
   expect_equal(round(sum(parList$rLPar[[2]]$RiWg)),444)
})
