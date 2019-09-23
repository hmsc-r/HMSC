context('Test updaters')

test_that("updateGamma2 is correct", {
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   dataParList = computeDataParameters(TD$m)
   Gamma = updateGamma2(Z=parList$Z,Gamma=parList$Gamma,iV=chol2inv(chol(parList$V)),iSigma=sqrt(parList$sigma),
                        Eta=parList$Eta,Lambda=parList$Lambda, X=TD$m$X,Pi=TD$m$Pi,dfPi=TD$m$dfPi,Tr=TD$m$Tr,
                        C=TD$m$C,rL=TD$m$rL, iQg=dataParList$iQg, mGamma=TD$m$mGamma,iUGamma=chol2inv(chol(TD$m$UGamma)))
   expect_equal(ncol(Gamma),3)
   expect_equal(nrow(Gamma),3)
   expect_equal(round(sum(Gamma)),0)
})

test_that("UpdateGammaEta is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   dataParList = computeDataParameters(TD$m)
   GammaEtaList = updateGammaEta(Z=parList$Z,Gamma=parList$Gamma,V=parList$V,iV=chol2inv(chol(parList$V)),id=sqrt(parList$sigma),
                                 Eta=parList$Eta,Lambda=parList$Lambda,Alpha=parList$Alpha, X=TD$m$X,Pi=TD$m$Pi,
                                 dfPi=TD$m$dfPi,Tr=TD$m$Tr,rL=TD$m$rL, rLPar=dataParList$rLPar,Q=dataParList$Qg[,,parList$rho],
                                 iQ=dataParList$iQg[,,parList$rho],RQ=dataParList$RQg[,,parList$rho],
                                 U=TD$m$UGamma,iU=chol2inv(chol(TD$m$UGamma)))
   gamma = GammaEtaList$Gamma
   eta = GammaEtaList$Eta
   expect_equal(length(eta),2)
   expect_equal(length(eta[[1]]),100)
   expect_equal(length(eta[[2]]),20)
   expect_equal(round(sum(eta[[1]])),-29)
   expect_equal(round(sum(eta[[2]])),2)
   expect_equal(length(gamma),9)
   expect_equal(round(sum(gamma)),0)
})

test_that("updateBetaLambda is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   dataParList = computeDataParameters(TD$m)
   BetaLambdaList = updateBetaLambda(Y=TD$Y,Z=parList$Z,Gamma=parList$Gamma,iV=chol2inv(chol(parList$V)),
                                     iSigma=sqrt(parList$sigma),Eta=parList$Eta,Psi=parList$Psi,Delta=parList$Delta,
                                     iQ=dataParList$iQg[,,parList$rho],X=TD$m$X,Tr=TD$m$Tr,Pi=TD$m$Pi,dfPi=TD$m$dfPi,C=TD$m$C,rL=TD$m$rL)
   Beta = BetaLambdaList$Beta
   Lambda = BetaLambdaList$Lambda
   expect_equal(length(Lambda),2)
   expect_equal(length(Lambda[[1]]),8)
   expect_equal(length(Lambda[[2]]),8)
   expect_equal(round(sum(Lambda[[1]]),digits=1),0.3)
   expect_equal(round(sum(Lambda[[2]]),digits=1),0.3)
   expect_equal(length(Beta),12)
   expect_equal(round(sum(Beta)),0)
})

test_that("updateGammaV is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   dataParList = computeDataParameters(TD$m)
   GammaVList = updateGammaV(Beta=parList$Beta,Gamma=parList$Gamma,iV=chol2inv(chol(parList$V)),rho=TD$m$rho,
                             iQg=dataParList$iQg,RQg=dataParList$RQg, Tr=TD$m$Tr,C=TD$m$C, mGamma=TD$m$mGamma,
                             iUGamma=chol2inv(chol(TD$m$UGamma)),V0=TD$m$V0,f0=TD$m$f0)
   Gamma = GammaVList$Gamma
   iV = GammaVList$iV
   expect_equal(nrow(iV),3)
   expect_equal(ncol(iV),3)
   expect_equal(round(sum(iV)),18)
   expect_equal(nrow(Gamma),3)
   expect_equal(ncol(Gamma),3)
   expect_equal(round(sum(Gamma)),1)
})

test_that("updateRho is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   dataParList = computeDataParameters(TD$m)
   rho = updateRho(Beta=parList$Beta,Gamma=parList$Gamma,iV=chol2inv(chol(parList$V)), RQg=dataParList$RQg,
                   detQg=dataParList$detQg, Tr=TD$m$Tr, rhopw=TD$m$rhopw)
   expect_equal(rho,1)
})

test_that("updateLambdaPriors is correct",{
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   PsiDeltaList = updateLambdaPriors(Lambda=parList$Lambda,Delta=parList$Delta, rL=TD$m$rL)
   Psi = PsiDeltaList$Psi
   Delta = PsiDeltaList$Delta
   expect_equal(length(Delta),2)
   expect_equal(length(Delta[[1]]),2)
   expect_equal(length(Delta[[2]]),2)
   expect_equivalent(round(Delta[[1]]),c(47,50))
   expect_equivalent(round(Delta[[2]]),c(59,43))
   expect_equal(length(Psi),2)
   expect_equal(round(sum(Psi[[1]])),11)
   expect_equal(round(sum(Psi[[2]])),12)
})


test_that("updateEta is correct", {
   set.seed(200)
   parList = computeInitialParameters(TD$m,initPar=NULL)
   dataParList = computeDataParameters(TD$m)
   eta = updateEta(Y = TD$m$Y, Z=parList$Z,Beta=parList$Beta,iSigma=sqrt(parList$sigma),Eta=parList$Eta,
             Lambda=parList$Lambda, Alpha = parList$Alpha, rLPar = dataParList$rLPar, X = TD$m$X, Pi = TD$m$Pi,
             dfPi = TD$m$dfPi,rL=TD$m$rL)
   expect_equal(length(eta),2)
   expect_equal(length(eta[[1]]),100)
   expect_equal(length(eta[[2]]),20)
   expect_equal(round(sum(eta[[1]])),-10)
   expect_equal(round(sum(eta[[2]])),-10)
})

test_that("updateAlpha is correct",{
   set.seed(200)
   dataParList = computeDataParameters(TD$m)
   Alpha = updateAlpha(Eta=TD$m$postList[[1]][[1]]$Eta, rLPar=dataParList$rLPar, rL=TD$m$rL)
   expect_equal(length(Alpha),2)
   expect_equal(length(Alpha[[1]]),2)
   expect_equal(length(Alpha[[2]]),2)
   expect_equal(Alpha[[1]],c(1,1))
   expect_equal(Alpha[[2]],c(19,1))
})

test_that("updateInvSigma is correct",{
   set.seed(200)
   Y=matrix(1:20,nrow=5,ncol=4)
   Y[,2]= sample(c(0,1),5,replace = TRUE)
   m = Hmsc(Y=Y,XData=data.frame(x1 = 1:5),distr = c('normal','probit','poisson','lognormal poisson'))
   parList = computeInitialParameters(m,initPar='fixed effects')
   dataParList = computeDataParameters(m)
   iSigma = updateInvSigma(Y=m$Y,Z=parList$Z,Beta=parList$Beta,iSigma=sqrt(parList$sigma),
                           Eta=parList$Eta,Lambda=parList$Lambda, distr=m$distr,X=m$X,Pi=m$Pi,
                           dfPi=m$dfPi,rL=m$rL, aSigma=m$aSigma,bSigma=m$bSigma)
   expect_equal(round(iSigma),c(0,1,0,0))
})

test_that("updateNf is correct",{
   set.seed(200)
   rL2 = HmscRandomLevel(units = TD$studyDesign$sample)
   m = Hmsc(Y=TD$Y,
               XData=TD$X,
               XFormula=~x1+x2,
               TrData=TD$Tr,
               TrFormula = ~T1 + T2,
               phyloTree=TD$phy,
               ranLevels=list("sample"=rL2),
               studyDesign = TD$studyDesign,
               distr=c('probit'))
   parList = computeInitialParameters(m,initPar=NULL)
   dataParList = computeDataParameters(m)
   listPar = updateNf(eta=parList$Eta[[1]],lambda=parList$Lambda[[1]],alpha=parList$Alpha[[1]],psi=parList$Psi[[1]],
                      delta=parList$Delta[[1]],rL=m$rL[[1]], iter=100)
   Lambda = listPar$lambda
   Eta = listPar$eta
   Alpha = listPar$alpha
   Psi = listPar$psi
   Delta = listPar$delta
   expect_equal(nrow(Lambda),3)
   expect_equal(ncol(Eta),3)
   expect_equal(length(Alpha),3)
   expect_equal(nrow(Psi),3)
   expect_equal(length(Delta),3)
})
## -----------------------------------------------------------------------------
context('Test sampleMCMC')

test_that("sampleMCMC returns m object of right size",{
   set.seed(200)
   m = sampleMcmc(TD$m,samples=1)
   expect_equal(length(m),72)
   expect_equal(length(m$postList[[1]][[1]]),13)
})

