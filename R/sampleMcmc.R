#' @title sampleMCMC
#'
#' @description Samples the posterior with block-conditional Gibbs MCMC sampler
#'
#' @param samples number of MCMC steps to be recorded
#' @param transient
#' @param thin thinning between recorded MCMC samples
#' @param initPar initial parameters value
#' @param saveToDisk whether to save replicates to the disk once they are ready and discard them from RAM memory
#' @param verbose at which steps should model progress be displayed (default = samples*thin / 100)
#' @param adaptNf
#' @param nChains number of MCMC chains to run
#' @param dataParList
#' @param updater
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export

sampleMcmc = function(hM, samples, transient=0, thin=1, initPar=NULL, verbose=samples*thin/100, adaptNf=rep(transient,hM$nr), nChains=1, nCores=1, dataParList=NULL, updater=list()){
   force(adaptNf)
   if(nCores > nChains){
      warning('number of cores cannot be more than number of chains')
      nCores <- nChains
   }

   X = hM$XScaled
   Tr = hM$TrScaled
   Y = hM$Y
   distr = hM$distr
   Pi = hM$Pi
   C = hM$C
   nr = hM$nr

   mGamma = hM$mGamma
   iUGamma = chol2inv(solve(hM$UGamma))
   V0 = hM$V0
   f0 = hM$f0
   aSigma = hM$aSigma
   bSigma = hM$bSigma
   rhopw = hM$rhopw

   if(is.null(dataParList))
      dataParList = computeDataParameters(hM)
   Qg = dataParList$Qg
   iQg = dataParList$iQg
   RiQg = dataParList$RiQg
   detQg = dataParList$detQg
   rLPar = dataParList$rLPar

   hM$postList = vector("list", nChains)
   hM$repList = vector("list", nChains)
   initSeed = sample.int(.Machine$integer.max, nChains)

   sampleChain = function(chain){
      if(nChains>1)
         print(sprintf("Computing chain %d", chain))
      set.seed(initSeed[chain])
      parList = computeInitialParameters(hM,initPar)

      Gamma = parList$Gamma
      V = parList$V
      iV = chol2inv(chol(V))
      Beta = parList$Beta
      sigma = parList$sigma
      iSigma = 1 / sigma
      Lambda = parList$Lambda
      Eta = parList$Eta
      Alpha = parList$Alpha
      Psi = parList$Psi
      Delta = parList$Delta
      rho = parList$rho
      Z = parList$Z

      postList = vector("list", samples)
      for(iter in 1:(transient+samples*thin)){ #  the parallel version fails on this line
         if(!identical(updater$Gamma2, FALSE))
            Gamma = updateGamma2(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,
               Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,Tr=Tr,C=C, iQg=iQg,
               mGamma=mGamma,iUGamma=iUGamma)

         if(!identical(updater$GammaEta, FALSE) && hM$nr>0){
            if(all(unlist(lapply(hM$rL, function(a) a$sDim)))){ # temporary patch for models with non-spatial levels
               GammaEtaList = updateGammaEta(Z=Z,V=chol2inv(chol(iV)),iV=iV,id=iSigma,
                  Eta=Eta,Lambda=Lambda,Alpha=Alpha, X=X,Pi=Pi,Tr=Tr,rL=hM$rL, rLPar=rLPar,Q=Qg[,,rho],iQ=iQg[,,rho],U=hM$UGamma,iU=iUGamma)
               Gamma = GammaEtaList$Gamma
               Eta = GammaEtaList$Eta
            }
         }

         if(!identical(updater$BetaLambda, FALSE)){
            BetaLambdaList = updateBetaLambda(Y=Y,Z=Z,Gamma=Gamma,iV=iV,
               iSigma=iSigma,Eta=Eta,Psi=Psi,Delta=Delta, iQ=iQg[,,rho],
               X=X,Tr=Tr,Pi=Pi,C=C)
            Beta = BetaLambdaList$Beta
            Lambda = BetaLambdaList$Lambda
         }

         if(!identical(updater$GammaV, FALSE)){
            GammaVList = updateGammaV(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho,
               iQg=iQg, Tr=Tr,C=C, mGamma=mGamma,iUGamma=iUGamma,V0=V0,f0=f0)
            Gamma = GammaVList$Gamma
            iV = GammaVList$iV
         }

         if(!is.null(hM$C) && !identical(updater$Rho, FALSE)){
            rho = updateRho(Beta=Beta,Gamma=Gamma,iV=iV, RiQg=RiQg,
               detQg=detQg, Tr=Tr, rhopw=rhopw)
            # print(rho)
         }

         if(!identical(updater$LambdaPriors, FALSE)){
            PsiDeltaList = updateLambdaPriors(Lambda=Lambda,Delta=Delta, rL=hM$rL)
            Psi = PsiDeltaList$Psi
            Delta = PsiDeltaList$Delta
         }

         if(!identical(updater$Eta, FALSE))
            Eta = updateEta(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,
               Lambda=Lambda,Alpha=Alpha, rLPar=rLPar, X=X,Pi=Pi,rL=hM$rL)

         if(!identical(updater$Alpha, FALSE))
            Alpha = updateAlpha(Eta=Eta, rLPar=rLPar, rL=hM$rL)

         if(!identical(updater$InvSigma, FALSE))
            iSigma = updateInvSigma(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,
               Eta=Eta,Lambda=Lambda, distr=distr,X=X,Pi=Pi, aSigma=aSigma,bSigma=bSigma)

         if(!identical(updater$Z, FALSE))
            Z = updateZ(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,distr=distr)

         for(r in seq_len(nr)){
            if(iter <= adaptNf[r]){
               listPar = updateNf(eta=Eta[[r]],lambda=Lambda[[r]],alpha=Alpha[[r]],psi=Psi[[r]],delta=Delta[[r]],
                  rL=hM$rL[[r]], iter=iter)
               Lambda[[r]] = listPar$lambda
               Eta[[r]] = listPar$eta
               Alpha[[r]] = listPar$alpha
               Psi[[r]] = listPar$psi
               Delta[[r]] = listPar$delta
            }
         }

         if((iter>transient) && ((iter-transient) %% thin == 0)){
            postList[[(iter-transient)/thin]] = combineParameters(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho,iSigma=iSigma,
               Eta=Eta,Lambda=Lambda,Alpha=Alpha,Psi=Psi,Delta=Delta,
               nc=hM$nc, XScalePar=hM$XScalePar, XInterceptInd=hM$XInterceptInd, nt=hM$nt, TrScalePar=hM$TrScalePar, TrInterceptInd=hM$TrInterceptInd, rhopw=rhopw)
         }
         if((verbose > 0) && (iter%%verbose == 0)){
            if(iter > transient){
               samplingStatusString = "sampling"
            } else{
               samplingStatusString = "transient"
            }
            print(sprintf("Chain %d, iteration %d of %d, (%s)", chain, iter, transient+samples*thin, samplingStatusString) )
         }
      }
      return(postList)
   }

   if(nCores > 1){
      cl = makeCluster(nCores, type="SOCK")
      clusterExport(cl, c("hM","nChains","transient","samples","thin","adaptNf","initSeed","initPar",
         "X", "Tr", "Y", "distr", "Pi", "C", "nr",
         "mGamma", "iUGamma", "V0", "f0", "aSigma", "bSigma", "rhopw",
         "Qg", "iQg", "RiQg", "detQg", "rLPar"), envir=environment())

      clusterEvalQ(cl, {
         library(mvtnorm);
         library(BayesLogit);
         library(MCMCpack);
         library(truncnorm);
         library(Matrix) })
      hM$postList = clusterApplyLB(cl, 1:nChains, fun=sampleChain)
      stopCluster(cl)
   } else{
      for(chain in 1:nChains){
         hM$postList[[chain]] = sampleChain(chain)
      }
   }


   hM$samples = samples
   hM$transient = transient
   hM$thin = thin
   hM$adaptNf = adaptNf

   return(hM)
}
