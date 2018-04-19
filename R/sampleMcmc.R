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
#' @export

sampleMcmc = function(samples, thin=1, initPar=NULL, repN=1, saveToDisk=FALSE, verbose=samples*thin/100, adaptNf=NULL){
   self$samples = samples
   self$thin = thin
   self$repN = repN
   self$saveToDisk = saveToDisk
   self$adaptNf = adaptNf

   X = self$X
   Tr = self$Tr
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

   dataParList = private$computeDataParameters()
   iQg = dataParList$iQg
   RiQg = dataParList$RiQg
   detQg = dataParList$detQg
   rLPar = dataParList$rLPar

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
      for(iter in 1:(samples*thin)){
         # Beta = updateBeta(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Tr=Tr,Pi=Pi)
         # Lambda = updateLambda(Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda,Psi=Psi,Delta=Delta, X=X,Pi=Pi,rL=self$rL)

         # if(nr>0){

         BetaLambdaList = updateBetaLambda(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,Eta=Eta,Psi=Psi,Delta=Delta,rho=rho, iQg=iQg, X=X,Tr=Tr,Pi=Pi,C=C)
         Beta = BetaLambdaList$Beta
         Lambda = BetaLambdaList$Lambda

         # } else{
            # Beta = updateBeta(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Tr=Tr,Pi=Pi)
         # }

         GammaVList = updateGammaV(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho, iQg=iQg, Tr=Tr,C=C, mGamma=mGamma,iUGamma=iUGamma,V0=V0,f0=f0)
         Gamma = GammaVList$Gamma
         iV = GammaVList$iV
         if(!is.null(self$C)){
            rho = updateRho(Beta=Beta,Gamma=Gamma,iV=iV, RiQg=RiQg,detQg=detQg, Tr=Tr, rhopw=rhopw)
            print(rho)
         }
         PsiDeltaList = updateLambdaPriors(Lambda=Lambda,Delta=Delta, rL=self$rL) #nu=nu,a1=a1,b1=b1,a2=a2,b2=b2)
         Psi = PsiDeltaList$Psi
         Delta = PsiDeltaList$Delta

         Eta = updateEta(Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda,Alpha=Alpha, rLPar=rLPar, X=X,Pi=Pi,rL=self$rL)
         Alpha = updateAlpha(Eta=Eta, rLPar=rLPar, rL=self$rL)
         iSigma = updateInvSigma(Z=Z,Beta=Beta,Eta=Eta,Lambda=Lambda, distr=distr,X=X,Pi=Pi, aSigma=aSigma,bSigma=bSigma)
         Z = updateZ(Y=Y,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,distr=distr)

         for(r in 1:nr){
            if( (is.list(adaptNf) && adaptNf[[r]][1] >= repN && adaptNf[[r]][2] >= iter) || (is.vector(adaptNf) && adaptNf[r] >= iter)){
               listPar = updateNf(eta=Eta[[r]],lambda=Lambda[[r]],psi=Psi[[r]],delta=Delta[[r]],
                  rL=self$rL[[r]], iter=iter)
               Lambda[[r]] = listPar$lambda
               Eta[[r]] = listPar$eta
               Psi[[r]] = listPar$psi
               Delta[[r]] = listPar$delta
            }
         }

         if(iter %% thin == 0){
            postList[[iter/thin]] = private$combineParameters(Beta=Beta,Gamma=Gamma,iV=iV,iSigma=iSigma,Eta=Eta,Lambda=Lambda,Alpha=Alpha,Psi=Psi,Delta=Delta)
         }
         if(iter %% verbose == 0){
            print(sprintf("Replicate %d, iteration %d of %d", repN, iter, samples*thin))
         }
      }
      repList[[repN]] = postList
      self$postList = postList
   }
   self$repList = repList
}

Hmsc$set("public", "sampleMcmc", sampleMcmc, overwrite=TRUE)

