#' @title sampleMCMC
#'
#' @description Samples the posterior with block-conditional Gibbs MCMC sampler
#'
#' @param samples the number of MCMC samples to be obtained in each chain
#' @param transient the number of MCMC steps that are executed before starting recording posterior samples
#' @param thin the number of MCMC steps between recording samples to the posterior
#' @param initPar a named list with parameter values that is used for initialization of MCMC states
#' @param verbose the interval between MCMC steps, when the MCMC progress is printed to the console
#' @param adaptNf a vector of length \eqn{n_r} with number of MCMC steps at which the adaptation of the number of latent
#'   factors is conducted
#' @param nChains number of independent MCMC chains to run
#' @param nParallel number of parallel processes that the chains are executed among
#' @param dataParList a named list with pre-computed \code{Qg}, \code{iQg}, \code{RQg}, \code{detQg}, \code{rLPar}
#'   parameters
#' @param updater a named list, specifying which conditional updaters should be ommitted
#' @param fromPrior whether prior (TRUE) or posterior (FALSE) is to be sampled
#'
#' @return an \class{Hmsc}-class object with chains of posterior samples added to the \code{postList} field
#'
#' @details The exact number of samples to be recorded in order to get a proper estimate of the full posterior with
#'   Gibbs MCMC algorithms, as well as the required thinning and cut-off of transient is very problem-specific and
#'   depends both on the model structure and the data itself. Therefore, in general it is very cahllenging to a priori
#'   provide an informed recommendation on what values should be used for a particular problem. A common recommended
#'   strategy involves executing the posterior sampling with MCMC with some guess of the values for these arguments,
#'   checking the properties of the obtained samples (primarily potential scale reduction factor and effective sample
#'   size), and adjusting the guess accordingly.
#'
#'   The value of 1 for \code{thin} argument means that at each MCMC step after the transient a sample is recorded.
#'
#'   Typically, the vlaue of \code{nParallel} equal to \code{nChains} leads to most efficient usage of available
#'   parallelization capacities. However, this may be not the case if the R is configured with multi-tread linear
#'   algebra libraries. For debug and test purposes, the \code{nParallel} shall be set to 1, since only in this case a
#'   detailization of the potentially encountered errors would be available.
#'
#'   The \code{dataParList} argument may be handy for large problems that needed to be refitted multiple times, e.g.
#'   with different prior values. In that case, the data parameters that are precomputed for Hmsc sampling scheme may
#'   require undesirably lot of storage space if they are saved for each of the model. Instead, they could be computed
#'   only once and then directly reused, therefore reducing the storing redundancy.
#'
#'   Some of available conditional updaters partially duplicate each other. In certain cases, the usage of all of them
#'   may lead to subotpimal performance, compared some subset of those. Then, it is possible to manually disable some of
#'   them, by adding a \code{$UPDATER_NAME=FALSE} pair to the {updater} argument. Another usage of this argument
#'   involves cases when some of the model parameters are known and have to be fixed. However, such tweaks of the
#'   sampling scheme should be done with caution, as if compromized they would lead to erroneuos results.
#'
#'
#' @seealso \code{\link{Hmsc}}
#'

#' @examples
#' Y = matrix(rnorm(100*10),100,10,dimnames=list(NULL,letters[1:10]))
#' X = matrix(rnorm(100*3),100,3,dimnames=list(NULL,sprintf("cov%d",1:3)))
#' m = Hmsc(Y=Y, XData=as.data.frame(X), XFormula=~cov1+cov3)
#' m = sampleMcmc(m, 100)
#'
#' @export

sampleMcmc = function(hM, samples, transient=0, thin=1, initPar=NULL,
                      verbose=samples*thin/100, adaptNf=rep(transient,hM$nr),
                      nChains=1, nParallel=1, dataParList=NULL, updater=list(),
                      fromPrior = FALSE){
   if(fromPrior)
      nParallel = 1
   force(adaptNf)
   if(nParallel > nChains){
      warning('number of cores cannot be more than number of chains')
      nParallel <- nChains
   }

   X = hM$XScaled
   Tr = hM$TrScaled
   Y = hM$YScaled
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
   RQg = dataParList$RQg
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
      for(iter in seq_len(transient+samples*thin)){

         if(!identical(updater$Gamma2, FALSE) && is.matrix(X))
            Gamma = updateGamma2(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,
               Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,Tr=Tr,C=C,rL=hM$rL, iQg=iQg,
               mGamma=mGamma,iUGamma=iUGamma)

         if(!identical(updater$GammaEta, FALSE) && hM$nr>0 && identical(mGamma,rep(0,hM$nc*hM$nt)) && is.matrix(X)){ # assumes mGamma = 0
            GammaEtaList = updateGammaEta(Z=Z,Gamma=Gamma,V=chol2inv(chol(iV)),iV=iV,id=iSigma,
               Eta=Eta,Lambda=Lambda,Alpha=Alpha, X=X,Pi=Pi,Tr=Tr,rL=hM$rL, rLPar=rLPar,Q=Qg[,,rho],iQ=iQg[,,rho],RQ=RQg[,,rho],U=hM$UGamma,iU=iUGamma)
            Gamma = GammaEtaList$Gamma
            Eta = GammaEtaList$Eta
         }

         if(!identical(updater$BetaLambda, FALSE)){
            BetaLambdaList = updateBetaLambda(Y=Y,Z=Z,Gamma=Gamma,iV=iV,
               iSigma=iSigma,Eta=Eta,Psi=Psi,Delta=Delta, iQ=iQg[,,rho],
               X=X,Tr=Tr,Pi=Pi,C=C,rL=hM$rL)
            Beta = BetaLambdaList$Beta
            Lambda = BetaLambdaList$Lambda
         }

         if(!identical(updater$GammaV, FALSE)){
            GammaVList = updateGammaV(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho,
               iQg=iQg,RQg=RQg, Tr=Tr,C=C, mGamma=mGamma,iUGamma=iUGamma,V0=V0,f0=f0)
            Gamma = GammaVList$Gamma
            iV = GammaVList$iV
         }

         if(!is.null(hM$C) && !identical(updater$Rho, FALSE)){
            rho = updateRho(Beta=Beta,Gamma=Gamma,iV=iV, RQg=RQg,
               detQg=detQg, Tr=Tr, rhopw=rhopw)
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
               Eta=Eta,Lambda=Lambda, distr=distr,X=X,Pi=Pi,rL=hM$rL, aSigma=aSigma,bSigma=bSigma)

         if(!identical(updater$Z, FALSE))
            Z = updateZ(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,distr=distr,rL=hM$rL)

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

   if(nParallel > 1){
      cl = makeCluster(nParallel, type="SOCK")
      clusterExport(cl, c("hM","nChains","transient","samples","thin","adaptNf","initSeed","initPar","updater",
         "X", "Tr", "Y", "distr", "Pi", "C", "nr",
         "mGamma", "iUGamma", "V0", "f0", "aSigma", "bSigma", "rhopw",
         "Qg", "iQg", "RQg", "detQg", "rLPar"), envir=environment())

      clusterEvalQ(cl, {
         library(mvtnorm);
         library(BayesLogit);
         library(MCMCpack);
         library(truncnorm);
         library(Matrix);
         library(abind);
         library(Hmsc)})
      hM$postList = clusterApplyLB(cl, 1:nChains, fun=sampleChain)
      stopCluster(cl)
   } else {
      for(chain in 1:nChains){
         if (fromPrior){
            postList = vector("list", samples)
            for (iter in 1:samples){
               postList[[iter]] = samplePrior(hM,dataParList = dataParList)
            }
            hM$postList[[chain]] = postList
         } else {
            hM$postList[[chain]] = sampleChain(chain)
         }
      }
   }

   hM$samples = samples
   hM$transient = transient
   hM$thin = thin
   hM$adaptNf = adaptNf

   return(hM)
}
