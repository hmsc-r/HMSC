#' @title sampleMCMC
#'
#' @description Samples the posterior with block-conditional Gibbs MCMC sampler
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param samples the number of MCMC samples to be obtained in each chain
#' @param transient the number of MCMC steps that are executed before starting recording posterior samples
#' @param thin the number of MCMC steps between each recording of samples from the posterior
#' @param initPar a named list of parameter values used for
#'     initialization of MCMC states, or alternatively text
#'     \code{"fixed effects"} to use linear Maximum Likelihood model
#'     instead of randomizing from prior; the \code{"fixed effects"}
#'     can shorten the transient phase of sampling, but will
#'     initialize all chains to the same starting values
#' @param verbose the interval between MCMC steps printed to the console (default is an interval that prints ca. 50 reports)
#' @param adaptNf a vector of length \eqn{n_r} with number of MCMC steps at which the adaptation of the
#' number of latent factors is conducted
#' @param nChains number of independent MCMC chains to be run
#'
#' @param nParallel number of parallel processes by which the chains
#'     are executed.
#' @param useSocket (logical) use socket clusters in parallel
#'     processing; in Windows this is the only option, but in other
#'     operating systems fork clusters are a better alternative, and
#'     this should be set \code{FALSE}.
#'
#' @param dataParList a named list with pre-computed \code{Qg}, \code{iQg}, \code{RQg}, \code{detQg}, \code{rLPar}
#'   parameters
#' @param updater a named list, specifying which conditional updaters should be ommitted
#' @param fromPrior whether prior (TRUE) or posterior (FALSE) is to be sampled
#' @param alignPost boolean flag indicating whether the posterior of each chains should be aligned
#'
#' @return An \code{Hmsc}-class object with chains of posterior samples added to the \code{postList} field
#'
#' @details The exact number of samples to be recorded in order to get a proper estimate of the full posterior with
#'   Gibbs MCMC algorithms, as well as the required thinning and cut-off of transient is very problem-specific and
#'   depends both on the model structure and the data itself. Therefore, in general it is very challenging to a priori
#'   provide an informed recommendation on what values should be used for a particular problem. A common recommended
#'   strategy involves executing the posterior sampling with MCMC with some guess of the values for these arguments,
#'   checking the properties of the obtained samples (primarily potential scale reduction factor and effective sample
#'   size), and adjusting the guess accordingly.
#'
#'   The value of 1 for \code{thin} argument means that at each MCMC step after the transient a sample is recorded.
#'
#'   Typically, the value of \code{nParallel} equal to \code{nChains} leads to most efficient usage of available
#'   parallelization capacities. However, this may be not the case if R is configured with multi-tread linear
#'   algebra libraries. For debug and test purposes, the \code{nParallel} should be set to 1, since only in this case a
#'   details of the potentially encountered errors would be available.
#'
#'   The \code{dataParList} argument may be handy for large problems that needs to be refitted multiple times, e.g.
#'   with different prior values. In that case, the data parameters that are precomputed for the Hmsc sampling
#'   scheme may require an undesirably lot of storage space if they are saved for each of the model.
#'   Instead, they could be computed only once and then directly reused, therefore reducing the storing redundancy.
#'
#'   Some of the available conditional updaters partially duplicate each other. In certain cases, the usage of all
#'   of them may lead to suboptimal performance, compared to some subset of those. Then, it is possible to manually
#'   disable some of them, by adding a \code{$UPDATER_NAME=FALSE} pair to the {updater} argument. Another usage of
#'   this argument involves cases when some of the model parameters are known and have to be fixed. However, such
#'   tweaks of the sampling scheme should be done with caution, as if compromized they would lead to erroneuos
#'   results.
#'
#'
#' @seealso \code{\link{Hmsc}}
#'
#' @examples
#' ## you need 1000 or more samples, but that will take too long
#' ## in an example
#' m = sampleMcmc(TD$m, samples=10)
#'
#' \dontrun{
#' ## Record 1000 posterior samples while skipping 1 MCMC step between samples
#' ## from 2 chains after discarding the first 500 MCMC steps
#' m = sampleMcmc(TD$m, samples=1000, transient=500, thin=2, nChains=2, nParallel=1)
#' }
#'
#' @importFrom parallel makeCluster clusterExport clusterEvalQ clusterApplyLB
#'     stopCluster mclapply
#' @export

sampleMcmc =
    function(hM, samples, transient=0, thin=1, initPar=NULL,
             verbose, adaptNf=rep(transient,hM$nr),
             nChains=1, nParallel=1,
             useSocket = TRUE,
             dataParList=NULL, updater=list(),
             fromPrior = FALSE, alignPost = TRUE)
{
   ## use socket cluster if requested or in Windows
   if (nParallel > 1 && .Platform$OS.type == "windows" && !useSocket) {
       useSocket <- TRUE
       message("only socket clusters can be used in Windows: setting useSocket = TRUE")
   }
   if (missing(verbose)) {
       if (samples*thin <= 50) # report every sampling
           verbose <- 1
       else                    # report ~50 steps of sampling
           verbose <- samples*thin/50
   }
   verbose <- as.integer(verbose) # truncate to integer
   if(fromPrior)
      nParallel = 1
   force(adaptNf)

   if(nParallel > nChains) {
      nParallel <- nChains
      message('using ', nParallel, ' cores for ', nChains, ' chains')
   }
   if(any(adaptNf > transient))
      stop('transient parameter should be no less than any element of adaptNf parameter')

   #X1 is the original X matrix (scaled version).
   #X used in computations is modified from X1 by variable selection and dimension reduction.
   X1 = hM$XScaled
   if (hM$ncsel>0){
      if(is.matrix(X1)){
         X2=X1
         X1=list()
         for (j in 1:hM$ns){
            X1[[j]] = X2
         }
      }
   }

   Tr = hM$TrScaled
   Y = hM$YScaled
   distr = hM$distr
   Pi = hM$Pi
   dfPi = hM$dfPi
   C = hM$C
   nr = hM$nr

   mGamma = hM$mGamma
   iUGamma = chol2inv(chol(hM$UGamma))
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
   ## save random seed that is used to generate initSeed[s]
   if (!exists(".Random.seed")) # may not exist, so generate
      runif(1)
   hM$randSeed <- .Random.seed
   initSeed = sample.int(.Machine$integer.max, nChains)

   ######## switching of the augmented updaters if the required conditions are not met
   EPS = 1e-6
   updaterWarningFlag = TRUE
   # updater$Gamma2
   if(!identical(updater$Gamma2, FALSE) && any(abs(mGamma) > EPS)){
      updater$Gamma2 = FALSE
      if(updaterWarningFlag)
         message("setting updater$Gamma2=FALSE due to non-zero mGamma")
   }
   if(!identical(updater$Gamma2, FALSE) && any(abs(iUGamma - kronecker(iUGamma[1:hM$nc,1:hM$nc], diag(hM$nt))) > EPS)){
      updater$Gamma2 = FALSE
      if(updaterWarningFlag)
         message("setting updater$Gamma2=FALSE due to non-kronecker structure of UGamma matrix")
   }
   if(!identical(updater$Gamma2, FALSE) && (!is.null(C))){
      updater$Gamma2 = FALSE
      if(updaterWarningFlag)
         message("setting updater$Gamma2=FALSE due to specified phylogeny matrix")
   }
   # updater$GammaEta
   if(!identical(updater$GammaEta, FALSE) && any(abs(mGamma) > EPS)){
      updater$GammaEta = FALSE
      if(updaterWarningFlag)
         message("setting updater$GammaEta=FALSE due to non-zero mGamma")
   }
   if(!identical(updater$GammaEta, FALSE) && hM$nr==0){
      updater$GammaEta = FALSE
      if(updaterWarningFlag)
         message("setting updater$GammaEta=FALSE due to absence of random effects included to the model")
   }
   # NNGP & GPP models will give an error in updateGammaEta()
   if (!identical(updater$GammaEta, FALSE) &&
       any(sapply(hM$rL,
                  function(s) !is.null(s$spatialMethod) &&
                  s$spatialMethod %in% c("GPP", "NNGP")))) {
       updater$GammaEta = FALSE
       if (updaterWarningFlag)
           message("setting updater$GammaEta=FALSE: not implemented for spatial methods 'GPP' and 'NNGP'")
   }

   ## latentLoadingOrderSwap: as of version 3.0.10 this is still an
   ## experimental feature and we do not advertise it by broadcasting
   ## messages that it is disabled
   if(identical(updater$latentLoadingOrderSwap, NULL)){
      updater$latentLoadingOrderSwap = 0
      if(FALSE && updaterWarningFlag) # do not advertise yet
         message("setting updater$latentLoadingOrderSwap=0 disabling full-conditional swapping of consecutive latent loadings")
   }


   sampleChain = function(chain){
      if(nChains>1)
         cat(sprintf("Computing chain %d\n", chain))
      set.seed(initSeed[chain])
      parList = computeInitialParameters(hM,initPar)

      Gamma = parList$Gamma
      V = parList$V
      iV = chol2inv(chol(V))
      Beta = parList$Beta
      BetaSel = parList$BetaSel
      PsiRRR = parList$PsiRRR
      DeltaRRR = parList$DeltaRRR
      wRRR = parList$wRRR
      sigma = parList$sigma
      iSigma = 1 / sigma
      Lambda = parList$Lambda
      Eta = parList$Eta
      Alpha = parList$Alpha
      Psi = parList$Psi
      Delta = parList$Delta
      rho = parList$rho
      Z = parList$Z

      X1A = X1

      if(hM$ncsel>0){
         for (i in 1:hM$ncsel){
            XSel = hM$XSelect[[i]]
            for (spg in 1:length(XSel$q)){
               if(!BetaSel[[i]][spg]){
                  fsp = which(XSel$spGroup==spg)
                  for (j in fsp){
                     X1A[[j]][,XSel$covGroup]=0
                  }
               }
            }
         }
      }

      X = X1A
      if(hM$ncRRR>0){
         XB=hM$XRRRScaled%*%t(wRRR)
         if(is.matrix(X)){
            X=cbind(X,XB)
         } else {
            for (j in 1:hM$ns){
               X[[j]] = cbind(X[[j]],XB)
            }
         }
      }

      if(!identical(updater$Gamma2, FALSE) && (!is.matrix(X))){
         updater$Gamma2 = FALSE
         if(updaterWarningFlag)
            message("setting updater$Gamma2=FALSE due to X is not a matrix")
      }
      if(!identical(updater$GammaEta, FALSE) && (!is.matrix(X))){
         updater$GammaEta = FALSE
         if(updaterWarningFlag)
            message("setting updater$GammaEta=FALSE due to X is not a matrix")
      }

      postList = vector("list", samples)
      failed <- numeric(14) # counts of failed try(update*())s
      names(failed) <- c("Gamma2", "GammaEta", "BetaLambda", "wRRR",
                         "BetaSel", "GammaV", "Rho", "LambdaPriors",
                         "wRRRPriors", "Eta", "Alpha",
                         "invSigma", "Nf", "LatentLoadingOrder")
      for(iter in seq_len(transient+samples*thin)){

         if(!identical(updater$Gamma2, FALSE)) {
            out = try(updateGamma2(Z=Z,Gamma=Gamma,iV=iV,iSigma=iSigma,
               Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,dfPi=dfPi,Tr=Tr,C=C,rL=hM$rL, iQg=iQg,
               mGamma=mGamma,iUGamma=iUGamma), silent = TRUE)
            if (!inherits(out, "try-error")) {
                Gamma <- out
            } else {
                failed["Gamma2"] <- failed["Gamma2"] + 1
            }
         }
         if(!identical(updater$GammaEta, FALSE)){
            GammaEtaList = try(updateGammaEta(Z=Z,Gamma=Gamma,V=chol2inv(chol(iV)),iV=iV,id=iSigma,
               Eta=Eta,Lambda=Lambda,Alpha=Alpha, X=X,Pi=Pi,dfPi=dfPi,Tr=Tr,rL=hM$rL, rLPar=rLPar,Q=Qg[,,rho],iQ=iQg[,,rho],RQ=RQg[,,rho],
               mGamma=mGamma,U=hM$UGamma,iU=iUGamma), silent = TRUE)
            if (!inherits(GammaEtaList, "try-error")) {
                Gamma = GammaEtaList$Gamma
                Eta = GammaEtaList$Eta
            } else {
                failed["GammaEta"] <- failed["GammaEta"] + 1
            }
         }

         if(!identical(updater$BetaLambda, FALSE)){
            BetaLambdaList = try(updateBetaLambda(Y=Y,Z=Z,Gamma=Gamma,iV=iV,
               iSigma=iSigma,Eta=Eta,Psi=Psi,Delta=Delta, iQ=iQg[,,rho],
               X=X,Tr=Tr,Pi=Pi,dfPi=dfPi,C=C,rL=hM$rL), silent = TRUE)
            if (!inherits(BetaLambdaList, "try-error")) {
                Beta = BetaLambdaList$Beta
                Lambda = BetaLambdaList$Lambda
            } else {
                failed["BetaLambda"] <- failed["BetaLambda"] + 1
            }
         }

         if(!identical(updater$wRRR, FALSE) &&  hM$ncRRR>0){
            wRRRXList = try(updatewRRR(Z=Z, Beta=Beta, iSigma=iSigma,
                                 Eta=Eta, Lambda=Lambda, X1A=X1A, XRRR=hM$XRRRScaled,
                                 Pi=Pi, dfPi=dfPi,rL = hM$rL, PsiRRR=PsiRRR, DeltaRRR=DeltaRRR), silent = TRUE)
            if (!inherits(wRRRXList, "try-error")) {
                wRRR = wRRRXList$wRRR
                X = wRRRXList$X
            } else {
                failed["wRRR"] <- failed["wRRR"] + 1
            }
         }

         if(!identical(updater$BetaSel, FALSE) &&  hM$ncsel>0){
            BetaSelXList = try(updateBetaSel(Z=Z,XSelect = hM$XSelect, BetaSel=BetaSel,Beta=Beta, iSigma=iSigma,
                                             Lambda=Lambda, Eta=Eta, X1=X1,Pi=Pi,dfPi=dfPi,rL=hM$rL), silent = TRUE)
            if (!inherits(BetaSelXList, "try-error")) {
                BetaSel = BetaSelXList$BetaSel
                X = BetaSelXList$X
            } else {
                failed["BetaSel"] <- failed["BetaSel"] + 1
            }
         }

         if(!identical(updater$GammaV, FALSE)){
            GammaVList = try(updateGammaV(Beta=Beta,Gamma=Gamma,iV=iV,rho=rho,
                                          iQg=iQg,RQg=RQg, Tr=Tr,C=C, mGamma=mGamma,iUGamma=iUGamma,V0=V0,f0=f0), silent = TRUE)
            if (!inherits(GammaVList, "try-error")) {
                Gamma = GammaVList$Gamma
                iV = GammaVList$iV
            } else {
                failed["GammaV"] <- failed["GammaV"] + 1
            }
         }

         if(!is.null(hM$C) && !identical(updater$Rho, FALSE)){
            out = try(updateRho(Beta=Beta,Gamma=Gamma,iV=iV, RQg=RQg,
                                detQg=detQg, Tr=Tr, rhopw=rhopw), silent = TRUE)
            if (!inherits(out, "try-error"))
                rho <- out
            else
                failed["Rho"] <- failed["Rho"] + 1
         }

         if(!identical(updater$LambdaPriors, FALSE)){
            PsiDeltaList = try(updateLambdaPriors(Lambda=Lambda,Delta=Delta,
                                                  rL=hM$rL), silent = TRUE)
            if (!inherits(PsiDeltaList, "try-error")) {
                Psi = PsiDeltaList$Psi
                Delta = PsiDeltaList$Delta
            } else {
                failed["LambdaPriors"] <- failed["LambdaPriors"] + 1
            }
         }
         if(!identical(updater$wRRRPriors, FALSE) &&  hM$ncRRR>0){
            PsiDeltaList = try(updatewRRRPriors(wRRR=wRRR,Delta=DeltaRRR,
                                           nu=hM$nuRRR,a1=hM$a1RRR,
                                           b1=hM$b1RRR,a2=hM$a2RRR,b2=hM$b2RRR),
                               silent = TRUE)
            if (!inherits(PsiDeltaList, "try-error")) {
                PsiRRR = PsiDeltaList$Psi
                DeltaRRR = PsiDeltaList$Delta
            } else {
                failed["wRRRPriors"] <- failed["wRRRPriors"] + 1
            }
         }

         if(!identical(updater$Eta, FALSE))
            out = try(updateEta(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,
                                Lambda=Lambda,Alpha=Alpha, rLPar=rLPar, X=X,
                                Pi=Pi,dfPi=dfPi,rL=hM$rL), silent = TRUE)
            if (!inherits(out, "try-error"))
                Eta <- out
            else
                failed["Eta"] <- failed["Eta"] + 1

         if(!identical(updater$Alpha, FALSE))
            out = try(updateAlpha(Eta=Eta, rLPar=rLPar, rL=hM$rL), silent = TRUE)
            if (!inherits(out, "try-error"))
                Alpha <- out
            else
                failed["Alpha"] <- failed["Alpha"] + 1

         if(!identical(updater$InvSigma, FALSE))
            out = try(updateInvSigma(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,
                                    Eta=Eta,Lambda=Lambda, distr=distr,X=X,
                                    Pi=Pi,dfPi=dfPi,rL=hM$rL, aSigma=aSigma,
                                    bSigma=bSigma), silent = TRUE)
            if (!inherits(out, "try-error"))
                iSigma <- out
            else
                failed["invSigma"] <- failed["invSigma"] + 1

         if(!identical(updater$Z, FALSE)){
            Z = updateZ(Y=Y,Z=Z,Beta=Beta,iSigma=iSigma,Eta=Eta,Lambda=Lambda, X=X,Pi=Pi,dfPi=dfPi,distr=distr,rL=hM$rL)
         }

         for(r in seq_len(nr)){
            if(iter <= adaptNf[r]){
               listPar = try(updateNf(eta=Eta[[r]],lambda=Lambda[[r]],alpha=Alpha[[r]],psi=Psi[[r]],delta=Delta[[r]],
                                      rL=hM$rL[[r]], iter=iter), silent = TRUE)
               if (!inherits(listPar, "try-error")) {
                   Lambda[[r]] = listPar$lambda
                   Eta[[r]] = listPar$eta
                   Alpha[[r]] = listPar$alpha
                   Psi[[r]] = listPar$psi
                   Delta[[r]] = listPar$delta
               } else {
                   failed["Nf"] <- failed["Nf"] + 1
               }
            }
         }

         if(updater$latentLoadingOrderSwap>0 && (iter %% updater$latentLoadingOrderSwap == 0)){
            for(r in seq_len(nr)){
               listPar = try(updateLatentLoadingOrder(eta=Eta[[r]],lambda=Lambda[[r]],alpha=Alpha[[r]],delta=Delta[[r]],rL=hM$rL[[r]]), silent = TRUE)
               if (!inherits(listPar, "try-error")) {
                   Lambda[[r]] = listPar$lambda
                   Eta[[r]] = listPar$eta
                   Alpha[[r]] = listPar$alpha
                   Delta[[r]] = listPar$delta
               } else {
                   failed["LatentLoadingOrder"] + failed["LatentLoadingOrder"] + 1
               }
            }
            PsiDeltaList = try(updateLambdaPriors(Lambda=Lambda,Delta=Delta, rL=hM$rL))
            if (!inherits(PsiDeltaList, "try-error")) {
                Psi = PsiDeltaList$Psi
                Delta = PsiDeltaList$Delta
            } else {
                failed["PsiDelta"] <- failed["PsiDelta"] + 1
            }
         }

         if((iter>transient) && ((iter-transient) %% thin == 0)){
            postList[[(iter-transient)/thin]] = combineParameters(Beta=Beta,BetaSel=BetaSel,wRRR = wRRR, Gamma=Gamma,iV=iV,rho=rho,iSigma=iSigma,
               Eta=Eta,Lambda=Lambda,Alpha=Alpha,Psi=Psi,Delta=Delta,
               PsiRRR=PsiRRR,DeltaRRR=DeltaRRR,
               ncNRRR=hM$ncNRRR, ncRRR=hM$ncRRR, ncsel = hM$ncsel, XSelect = hM$XSelect,
               XScalePar=hM$XScalePar, XInterceptInd=hM$XInterceptInd, XRRRScalePar=hM$XRRRScalePar,
               nt=hM$nt, TrScalePar=hM$TrScalePar, TrInterceptInd=hM$TrInterceptInd, rhopw=rhopw)
         }
         postList$failedUpdates <- failed
         if((verbose > 0) && (iter%%verbose == 0)){
            if(iter > transient){
               samplingStatusString = "sampling"
            } else{
               samplingStatusString = "transient"
            }
            cat(sprintf("Chain %d, iteration %d of %d (%s)\n", chain, iter, transient+samples*thin, samplingStatusString) )
         }
      }
      return(postList)
   }

    if (nParallel > 1) {
        if (useSocket) {
            cl = makeCluster(nParallel)
            clusterExport(cl, c("hM","nChains","transient","samples","thin","verbose",
                                "adaptNf","initSeed","initPar","updater",
                                "X1", "Tr", "Y", "distr", "Pi", "C", "nr",
                                "mGamma", "iUGamma", "V0", "f0", "aSigma", "bSigma",
                                "rhopw","Qg", "iQg", "RQg", "detQg", "rLPar"),
                          envir=environment())
            clusterEvalQ(cl, {
                library(BayesLogit);
                library(MCMCpack);
                library(truncnorm);
                library(Matrix);
                library(abind);
                library(Hmsc)})
            hM$postList = clusterApplyLB(cl, 1:nChains, fun=sampleChain)
            stopCluster(cl)
        } else { # fork cluster
            verbose <- 0
            hM$postList <- mclapply(seq_len(nChains), function(i) sampleChain(i),
                                    mc.cores = nParallel)
        }
    } else {
       for(chain in 1:nChains ){
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
   ## warn on failed updaters
   for(chain in seq_len(nChains)) {
       ntries <- transient + samples * thin
       if (any(hM$postList[[chain]]$failedUpdates > 0)) {
           cat("Failed updaters and their counts in chain", chain,
               " (", ntries, " attempts):\n")
           failures <- hM$postList[[chain]]$failedUpdates
           failures <- failures[failures > 0]
           print(failures)
       }
       attr(hM$postList[[chain]], "failedUpdates") <-
           hM$postList[[chain]]$failedUpdates     # save as an attribute
       hM$postList[[chain]]$failedUpdates <- NULL # remove from postList
   }
   hM$samples = samples
   hM$transient = transient
   hM$thin = thin
   hM$verbose = verbose
   hM$adaptNf = adaptNf
   if (alignPost){
      for (i in 1:5){
         hM = alignPosterior(hM)
      }
   }

   return(hM)
}
