#' @title sampleMCMC
#'
#' @description Samples the posterior with block-conditional Gibbs MCMC sampler
#'
#' @param samples number of MCMC steps to be recorded
#' @param thin thinning between recorded MCMC samples
#' @param initPar initial parameters value
#' @param replicates number of replicates the samples should be splitted
#' @param saveToDisk whether to save replicates to the disk once they are ready and discard them from RAM memory
#' @param discardDataPar whether to discard the precomputed data-based quantities once the sampling is finished
#'
#' @examples
#'

sampleMcmc = function(samples, transient=0, thin=1, initPar=NULL, repN=1, saveToDisk=FALSE, verbose=samples*thin/100, adaptNf=NULL, nChains=1, dataParList=NULL, updater=list()){
   X = self$XScaled
   Tr = self$TrScaled
   Y = self$Y
   distr = self$distr
   Pi = self$Pi
   C = self$C
   nr = self$nr

   mGamma = self$mGamma
   iUGamma = chol2inv(solve(self$UGamma))
   V0 = self$V0
   f0 = self$f0
   aSigma = self$aSigma
   bSigma = self$bSigma
   nu = self$nu
   a1 = self$a1
   b1 = self$b1
   a2 = self$a2
   b2 = self$b2
   rhopw = self$rhopw

   if(is.null(dataParList))
      dataParList = private$computeDataParameters()
   iQg = dataParList$iQg
   RiQg = dataParList$RiQg
   detQg = dataParList$detQg
   rLPar = dataParList$rLPar

   self$postList = vector("list", nChains)
   self$repList = vector("list", nChains)
   initSeed = sample.int(.Machine$integer.max, nChains)
   for(chain in 1:nChains){
      if(nChains>1)
         print(sprintf("Computing chain %d", chain))
      set.seed(initSeed[chain])
      parList = private$computeInitialParameters(initPar)

      Gamma = parList$Gamma
      V = parList$V
      iV = solve(V)
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

      repList = vector("list", repN)
      for(repN in 1:repN){
         postList = vector("list", samples)
         for(iter in 1:(transient+samples*thin)){
            if(!identical(updater$Gamma2, FALSE))
               Gamma = updateGamma2(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,Tr=Tr,C=C, iQg=iQg, mGamma=mGamma,iUGamma=iUGamma)

            if(!identical(updater$BetaLambda, FALSE)){
               BetaLambdaList = updateBetaLambda(Y=Y,Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,Eta=Eta,Psi=Psi,Delta=Delta,rho=rho, iQg=iQg, X=X,Tr=Tr,Pi=Pi,C=C)
               Beta = BetaLambdaList$Beta
               Lambda = BetaLambdaList$Lambda
            }

            if(!identical(updater$BetaLambda, FALSE)){
               GammaVList = updateGammaV(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho, iQg=iQg, Tr=Tr,C=C, mGamma=mGamma,iUGamma=iUGamma,V0=V0,f0=f0)
               Gamma = GammaVList$Gamma
               iV = GammaVList$iV
            }

            if(!is.null(self$C) && !identical(updater$Rho, FALSE)){
               rho = updateRho(Beta=Beta,Gamma=Gamma,iV=iV, RiQg=RiQg,detQg=detQg, Tr=Tr, rhopw=rhopw)
               # print(rho)
            }

            if(!identical(updater$LambdaPriors, FALSE)){
               PsiDeltaList = updateLambdaPriors(Lambda=Lambda,Delta=Delta, rL=self$rL) #nu=nu,a1=a1,b1=b1,a2=a2,b2=b2)
               Psi = PsiDeltaList$Psi
               Delta = PsiDeltaList$Delta
            }

            if(!identical(updater$Eta, FALSE))
               Eta = updateEta(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda,Alpha=Alpha, rLPar=rLPar, X=X,Pi=Pi,rL=self$rL)

            if(!identical(updater$Alpha, FALSE))
               Alpha = updateAlpha(Eta=Eta, rLPar=rLPar, rL=self$rL)

            if(!identical(updater$InvSigma, FALSE))
               iSigma = updateInvSigma(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda, distr=distr,X=X,Pi=Pi, aSigma=aSigma,bSigma=bSigma)

            if(!identical(updater$Z, FALSE))
               Z = updateZ(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,distr=distr)

            for(r in seq_len(nr)){
               if( (is.list(adaptNf) && adaptNf[[r]][1] >= repN && adaptNf[[r]][2] >= iter) || (is.vector(adaptNf) && adaptNf[r] >= iter)){
                  listPar = updateNf(eta=Eta[[r]],lambda=Lambda[[r]],alpha=Alpha[[r]],psi=Psi[[r]],delta=Delta[[r]],
                     rL=self$rL[[r]], iter=iter)
                  Lambda[[r]] = listPar$lambda
                  Eta[[r]] = listPar$eta
                  Alpha[[r]] = listPar$alpha
                  Psi[[r]] = listPar$psi
                  Delta[[r]] = listPar$delta
               }
            }

            if((iter>transient) && ((iter-transient) %% thin == 0)){
               postList[[(iter-transient)/thin]] = private$combineParameters(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho,iSigma=iSigma,
                  Eta=Eta,Lambda=Lambda,Alpha=Alpha,Psi=Psi,Delta=Delta,rhopw=rhopw)
            }
            if((verbose > 0) && (iter%%verbose == 0)){
               if(iter > transient){
                  samplingStatusString = "sampling"
               } else{
                  samplingStatusString = "transient"
               }
               print(sprintf("Chain %d, replicate %d, iteration %d of %d, (%s)", chain, repN, iter, transient+samples*thin, samplingStatusString) )
            }
         }
         repList[[repN]] = postList
         self$postList[[chain]] = postList
      }
      self$repList[[chain]] = repList
   }

   self$samples = samples
   self$transient = transient
   self$thin = thin
   self$repN = repN
   self$saveToDisk = saveToDisk
   self$adaptNf = adaptNf
}

Hmsc$set("public", "sampleMcmc", sampleMcmc, overwrite=TRUE)

