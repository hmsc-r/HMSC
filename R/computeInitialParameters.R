#' @title computeInitialParameters
#'
#' @description Computes initial parameter values before the sampling starts
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param initPar a list of initial parameter values
#'
#' @return a list of Hmsc model parameters
#'
#' @importFrom stats glm.fit lm.fit coef rgamma rnorm binomial poisson
#'    cov runif
#' @importFrom abind abind
#' @importFrom MASS mvrnorm
#' @importFrom MCMCpack riwish


computeInitialParameters = function(hM, initPar){
   parList = list()

   if(hM$ncRRR>0){
      DeltaRRR = matrix(c(rgamma(1,hM$a1RRR,hM$b1RRR), rgamma(hM$ncRRR-1,hM$a2RRR,hM$b2RRR)))
      PsiRRR = matrix(rgamma(hM$ncRRR*hM$ncORRR, hM$nuRRR/2, hM$nuRRR/2), hM$ncRRR, hM$ncORRR)
      tauRRR = matrix(apply(DeltaRRR, 2, cumprod), hM$ncRRR, 1)
      tauMatRRR = matrix(tauRRR,hM$ncRRR,hM$ncORRR)
      multRRR = sqrt(PsiRRR*tauMatRRR)^-1
      wRRR = matrix(rnorm(hM$ncRRR*hM$ncORRR)*multRRR, hM$ncRRR, hM$ncORRR)
      XB=hM$XRRRScaled%*%t(wRRR)
   } else {
      wRRR = NULL
      PsiRRR = NULL
      DeltaRRR = NULL
   }

   switch(class(hM$X)[1L],
          matrix = {
             XScaled = hM$XScaled
             if(hM$ncRRR>0){
                XScaled=cbind(XScaled,XB)
             }
          },
          list = {
             XScaled=list()
             for(j in 1:hM$ns){
                XScaled[[j]] = hM$XScaled[[j]]
                if(hM$ncRRR>0){
                   XScaled[[j]]=cbind(XScaled[[j]],XB)
                }
             }
          }
   )

   if(identical(initPar, "fixed effects")){
      cat("Hmsc::computeInitialParameter - initializing fixed effects with SSDM estimates\n")
      Beta = matrix(NA,hM$nc,hM$ns)
      for(j in 1:hM$ns){
         switch(class(hM$X)[1L],
            matrix = {
              XEff = XScaled
            },
            list = {
               XEff = XScaled[[j]]
            }
            )
         ## Y can contain NA values
         kk <- !is.na(hM$Y[,j])
         if(hM$distr[j,1] == 1)
            fm = lm.fit(XEff[kk,, drop=FALSE], hM$Y[kk,j])
         if(hM$distr[j,1] == 2)
            fm = glm.fit(XEff[kk,, drop=FALSE], hM$Y[kk,j], family=binomial(link="probit"))
         if(hM$distr[j,1] == 3)
            fm = glm.fit(XEff[kk,, drop=FALSE], hM$Y[kk,j], family=poisson())
         Beta[,j] = coef(fm)
      }
      Gamma = matrix(NA,hM$nc,hM$nt)
      for(k in 1:hM$nc){
         fm = lm.fit(hM$Tr, Beta[k,])
         Gamma[k,] = coef(fm)
      }
      V = rowSums(abind(cov(t(Beta-tcrossprod(Gamma,hM$Tr))), diag(hM$nc), along=3), na.rm=TRUE, dims=2)

      initPar = NULL
   } else{

      if(!is.null(initPar$Gamma)){
         Gamma = initPar$Gamma
      } else{
         Gamma = matrix(mvrnorm(1, hM$mGamma, hM$UGamma), hM$nc, hM$nt)
      }

      if(!is.null(initPar$V)){
         V = initPar$V
      } else{
         V = riwish(hM$f0, hM$V0)
      }

      if(!is.null(initPar$Beta)){
         Beta = initPar$Beta
      } else{
         Beta = matrix(NA, hM$nc, hM$ns)
         Mu = tcrossprod(Gamma,hM$Tr)
         for(j in 1:hM$ns){
            Beta[,j] = mvrnorm(1, Mu[,j], V)
         }
      }
   }

   BetaSel=list()
   for (i in seq_len(hM$ncsel)){
      XSel = hM$XSelect[[i]]
      BetaSel[[i]] = runif(length(XSel$q))<XSel$q
   }

   if(!is.null(initPar$sigma)){
      sigma = initPar$sigma
   } else{
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
   }

   nf = rep(NA, hM$nr)
   ncr = rep(NA, hM$nr)
   Delta = vector("list", hM$nr)
   Psi = vector("list", hM$nr)
   Lambda = vector("list", hM$nr)
   Eta = vector("list", hM$nr)
   np = hM$np
   if(!is.null(initPar$Delta)){
      Delta = initPar$Delta
   } else{
      Delta = vector("list", hM$nr)
   }
   if(!is.null(initPar$Psi)){
      Psi = initPar$Psi
   } else{
      Psi = vector("list", hM$nr)
   }
   if(!is.null(initPar$Lambda)){
      Lambda = initPar$Lambda
   } else{
      Lambda = vector("list", hM$nr)
   }
   if(!is.null(initPar$Eta)){
      Eta = initPar$Eta
   } else{
      Eta = vector("list", hM$nr)
   }

   for(r in seq_len(hM$nr)){
      if(!is.null(initPar$Delta[[r]])){
         nf[r] = nrow(Delta[[r]])
      }
      if(!is.null(initPar$Psi[[r]])){
         nf[r] = nrow(Psi[[r]])
      }
      if(!is.null(initPar$Lambda[[r]])){
         nf[r] = nrow(Lambda[[r]])
      }
      if(!is.null(initPar$Eta[[r]])){
         nf[r] = ncol(Eta[[r]])
      }
      if(is.na(nf[r]))
         nf[r] = hM$rL[[r]]$nfMin

      ncr[r] = max(hM$rL[[r]]$xDim, 1)
      if(is.null(initPar$Delta[[r]])){
         if(hM$rL[[r]]$xDim == 0){
            Delta[[r]] = matrix(c(rgamma(1,hM$rL[[r]]$a1,hM$rL[[r]]$b1), rgamma(nf[r]-1,hM$rL[[r]]$a2,hM$rL[[r]]$b2)))
         } else{
            Delta[[r]] = matrix(c(rgamma(ncr[r],hM$rL[[r]]$a1,hM$rL[[r]]$b1), rgamma(ncr[r]*(nf[r]-1),hM$rL[[r]]$a2,hM$rL[[r]]$b2)), nf[r],ncr[r],byrow=TRUE)
         }

      }
      if(is.null(initPar$Psi[[r]])){
         if(hM$rL[[r]]$xDim == 0){
            Psi[[r]] = matrix(rgamma(nf[r]*hM$ns, hM$rL[[r]]$nu/2, hM$rL[[r]]$nu/2), nf[r], hM$ns)
         } else {
            Psi[[r]] = array(rgamma(nf[r]*hM$ns*ncr[r], hM$rL[[r]]$nu/2, hM$rL[[r]]$nu/2), dim=c(nf[r],hM$ns,ncr[r]))
         }
      }
      if(is.null(initPar$Lambda[[r]])){
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
      }
   }


   for(r in seq_len(hM$nr)){
      if(!is.null(initPar$Eta[[r]])){
         Eta[[r]] = initPar$Eta[[r]]
      } else{
         Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
      }
   }
   if(!is.null(initPar$Eta)){
      Alpha = initPar$Alpha
   } else{
      Alpha = vector("list", hM$nr)
   }
   for(r in seq_len(hM$nr)){
      if(!is.null(initPar$Alpha[[r]])){
         Alpha[[r]] = initPar$Alpha[[r]]
      } else{
         Alpha[[r]] = rep(1,nf[r])
      }
   }

   if(!is.null(initPar$rho)){
      rho = which.min(abs(initPar$rho-hM$rhopw[,1]))
   } else{
      rho = 1
   }

   switch(class(XScaled)[1L],
      matrix = {
         LFix = XScaled %*% Beta
      },
      list = {
         LFix = matrix(NA,hM$ny,hM$ns)
         for(j in 1:hM$ns)
            LFix[,j] = XScaled[[j]] %*% Beta[,j]
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
   if(hM$nr > 0){
      Z = LFix + Reduce("+", LRan)
   } else
      Z = LFix

   Z = updateZ(Y=hM$Y,Z=Z,Beta=Beta,iSigma=sigma^-1,Eta=Eta,Lambda=Lambda, X=XScaled,Pi=hM$Pi,dfPi=hM$dfPi,distr=hM$distr,rL=hM$rL)

   parList$Gamma = Gamma
   parList$V = V
   parList$Beta = Beta
   parList$BetaSel = BetaSel
   parList$PsiRRR = PsiRRR
   parList$DeltaRRR = DeltaRRR
   parList$wRRR = wRRR
   parList$sigma = sigma
   parList$Eta = Eta
   parList$Lambda = Lambda
   parList$Psi = Psi
   parList$Delta = Delta
   parList$Alpha = Alpha
   parList$rho = rho
   parList$Z = Z

   return(parList)
}

