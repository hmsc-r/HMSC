#' @title samplePrior
#'
#' @description Samples the parameter vector from prior
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param dataParList list of data parameters (see \code{\link{computeDataParameters}})
#'
#' @return A named list containing the Hmsc model parameters
#'
#' @importFrom stats rgamma rnorm
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack riwish

samplePrior = function(hM, dataParList=NULL){

   if(is.null(dataParList))
      dataParList = computeDataParameters(hM)
   Qg = dataParList$Qg
   iQg = dataParList$iQg
   RQg = dataParList$RQg
   detQg = dataParList$detQg
   rLPar = dataParList$rLPar


   #PRIMARY PARAMETERS NOT RELATED TO RANDOM EFFECTS
   Gamma = matrix(mvrnorm(1, hM$mGamma, hM$UGamma), hM$nc, hM$nt)

   V = riwish(hM$f0, hM$V0)

   sigma = rep(NA, hM$ns)
   for(j in 1:hM$ns){
      if(hM$distr[j,2] == 1){
         sigma[j] = 1 / rgamma(1, shape=hM$aSigma[j], rate=hM$bSigma[j])
      } else{
         switch(hM$distr[j,1],
                sigma[j] <- 1,
                sigma[j] <- 1,
                sigma[j] <- 1e-2
         )
      }
   }

   if(is.null(hM$C)){
      rho = 1
   } else {
      rho = sample(x = 1:dim(hM$rhopw)[1], size = 1, prob = hM$rhopw[,2])
   }

   #PRIMARY PARAMETERS FOR RANDOM EFFECTS

   nf = rep(NA, hM$nr)
   ncr = rep(NA, hM$nr)
   Delta = vector("list", hM$nr)
   Psi = vector("list", hM$nr)
   Lambda = vector("list", hM$nr)
   Eta = vector("list", hM$nr)
   np = hM$np
   Alpha = vector("list", hM$nr)

   for(r in seq_len(hM$nr)){
      if(hM$rL[[r]]$nfMax==Inf){
         nf[r] = 10
      } else {
         nf[r] = hM$rL[[r]]$nfMax
      }

      ncr[r] = max(hM$rL[[r]]$xDim, 1)
      if(hM$rL[[r]]$xDim == 0){
         Delta[[r]] = matrix(c(rgamma(1,hM$rL[[r]]$a1,hM$rL[[r]]$b1), rgamma(nf[r]-1,hM$rL[[r]]$a2,hM$rL[[r]]$b2)))
      } else
         Delta[[r]] = matrix(c(rgamma(ncr[r],hM$rL[[r]]$a1,hM$rL[[r]]$b1), rgamma(ncr[r]*(nf[r]-1),hM$rL[[r]]$a2,hM$rL[[r]]$b2)), nf[r],ncr[r],byrow=TRUE)

      if(hM$rL[[r]]$xDim == 0){
         Psi[[r]] = matrix(rgamma(nf[r]*hM$ns, hM$rL[[r]]$nu/2, hM$rL[[r]]$nu/2), nf[r], hM$ns)
      } else
         Psi[[r]] = array(rgamma(nf[r]*hM$ns*ncr[r], hM$rL[[r]]$nu/2, hM$rL[[r]]$nu/2), dim=c(nf[r],hM$ns,ncr[r]))

      tau = matrix(apply(Delta[[r]], 2, cumprod), nf[r], ncr[r])
      if(hM$rL[[r]]$xDim == 0){
         tauMat = matrix(tau,nf[r],hM$ns)
         mult = sqrt(Psi[[r]]*tauMat)^-1
         Lambda[[r]] = matrix(rnorm(nf[r]*hM$ns)*mult, nf[r], hM$ns)
      } else{
         tauArray = array(tau,dim=c(nf[r],1,ncr[r]))[,rep(1,hM$ns),,drop=FALSE]
         mult = sqrt(Psi[[r]]*tauArray)^-1
         Lambda[[r]] = array(rnorm(nf[r]*hM$ns*ncr[r])*mult, dim=c(nf[r],hM$ns,ncr[r]))
      }
      if(hM$rL[[r]]$sDim == 0){
         Alpha[[r]] = rep(1,nf[r])
         Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
      } else {
         Alpha[[r]] = sample(x = 1:dim(hM$rL[[r]]$alphapw)[1], size = nf[r], prob = hM$rL[[r]]$alphapw[,2], replace = TRUE)
         Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
         Wg = rLPar[[r]]$Wg
         alpha = Alpha[[r]]
         for(i in 1:nf[r]){
            Eta[[r]][,i]=mvrnorm(mu=rep(0,np[r]),Sigma=Wg[,,alpha[i]])
         }
      }
   }

   #DERIVED PARAMETERS

   Mu = tcrossprod(Gamma,hM$TrScaled)
   if(is.null(hM$C)){
      Beta = matrix(NA, hM$nc, hM$ns)
      for(j in 1:hM$ns)
         Beta[,j] = mvrnorm(1, Mu[,j], V)
   }
   else {
      Beta = t(matrix(mvrnorm(mu=as.vector(t(Mu)),Sigma = kronecker(V,Qg[,,rho])), hM$ns, hM$nc))
   }

   switch(class(hM$XScaled)[1L],
          matrix = {
             LFix = hM$XScaled %*% Beta
          },
          list = {
             LFix = matrix(NA,hM$ny,hM$ns)
             for(j in 1:hM$ns)
                LFix[,j] = hM$XScaled[[j]] %*% Beta[,j]
          }
   )
   LRan = vector("list", hM$nr)
   for(r in seq_len(hM$nr)){
      if(hM$rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][hM$Pi[,r],]%*%Lambda[[r]]
      } else{
         LRan[[r]] = matrix(0,hM$ny,hM$ns)
         for(k in 1:hM$rL[[r]]$xDim)
            LRan[[r]] = LRan[[r]] + (Eta[[r]][hM$Pi[,r],,drop=FALSE]*hM$rL[[r]]$x[as.character(hM$dfPi[,r]),r]) %*% Lambda[[r]][,,k]
      }
   }

   iSigma = 1 / sigma
   iV = chol2inv(chol(V))

   sample = combineParameters(Beta=Beta,BetaSel=NULL,wRRR=NULL,Gamma=Gamma,iV=iV,rho=rho,iSigma=iSigma,
                     Eta=Eta,Lambda=Lambda,Alpha=Alpha,Psi=Psi,Delta=Delta,
                     PsiRRR=NULL, DeltaRRR=NULL,ncNRRR=hM$ncNRRR, ncRRR=hM$ncRRR, ncsel=hM$ncsel, XSelect=NULL,
                     XScalePar=hM$XScalePar, XInterceptInd=hM$XInterceptInd, nt=hM$nt, TrScalePar=hM$TrScalePar,
                     TrInterceptInd=hM$TrInterceptInd, rhopw=hM$rhopw)
   return(sample)
}

