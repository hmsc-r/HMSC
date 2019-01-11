#' @title computeInitialParameters
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'


computeInitialParameters = function(hM, initPar){
   parList = list()

   if(identical(initPar, "fixed effects")){
      cat("Hmsc::computeInitialParameter - initializing fixed effects with SSDM estimates\n")
      Beta = matrix(NA,hM$nc,hM$ns)
      for(j in 1:hM$ns){
         switch(class(hM$X),
            matrix = {
              XEff = hM$X
            },
            list = {
               XEff = hM$X[[j]]
            }
         )
         if(hM$distr[j,1] == 1)
            fm = lm.fit(XEff, hM$Y[,j])
         if(hM$distr[j,1] == 2)
            fm = glm.fit(XEff, hM$Y[,j], family=binomial(link="probit"))
         if(hM$distr[j,1] == 3)
            fm = glm.fit(XEff, hM$Y[,j], family=poisson())
         Beta[,j] = coef(fm)
      }
      Gamma = matrix(NA,hM$nc,hM$nt)
      for(k in 1:hM$nc){
         fm = lm.fit(hM$Tr, Beta[k,])
         Gamma[k,] = coef(fm)
      }
      V = cov(t(Beta - tcrossprod(Gamma, hM$Tr)))
      initPar = NULL
   } else{

      if(!is.null(initPar$Gamma)){
         Gamma = initPar$Gamma
      } else{
         Gamma = matrix(rmvnorm(1, hM$mGamma, hM$UGamma), hM$nc, hM$nt)
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
         for(j in 1:hM$ns)
            Beta[,j] = rmvnorm(1, Mu[,j], V)
      }
   }

   if(!is.null(initPar$sigma)){
      sigma = initPar$sigma
   } else{
      sigma = rep(NA, hM$ns)
      for(j in 1:hM$ns){
         if(hM$distr[j,2] == 1){
            sigma[j] = rgamma(1, shape=hM$aSigma[j], rate=hM$bSigma[j])
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
   Delta = vector("list", hM$nr)
   Psi = vector("list", hM$nr)
   Lambda = vector("list", hM$nr)
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
      if(is.na(nf[r]))
         nf[r] = 2

      if(is.null(initPar$Delta[[r]])){
         Delta[[r]] = matrix(c(rgamma(1,hM$rL[[r]]$a1,hM$rL[[r]]$b1), rgamma(nf[r]-1,hM$rL[[r]]$a2,hM$rL[[r]]$b2)))
      }
      if(is.null(initPar$Psi[[r]])){
         Psi[[r]] = matrix(rgamma(nf[r]*hM$ns, hM$rL[[r]]$nu/2, hM$rL[[r]]$nu/2), nf[r], hM$ns)
      }
      if(is.null(initPar$Lambda[[r]])){
         tau = apply(Delta[[r]], 2, cumprod)
         tauMat = matrix(tau,nf[r],hM$ns)
         mult = sqrt(Psi[[r]]*tauMat)^-1
         Lambda[[r]] = matrix(rnorm(nf[r]*hM$ns)*mult, nf[r], hM$ns)
      }
   }

   Eta = vector("list", hM$nr)
   np = hM$np

   if(!is.null(initPar$Eta)){
      Eta = initPar$Eta
   } else{
      Eta = vector("list", hM$nr)
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

   switch(class(hM$X),
      matrix = {
         LFix = hM$X %*% Beta
      },
      list = {
         LFix = matrix(NA,hM$ny,hM$ns)
         for(j in 1:hM$ns)
            LFix[,j] = hM$X[[j]] %*% Beta[,j]
      }
   )
   LRan = vector("list", hM$nr)
   for(r in seq_len(hM$nr)){
      LRan[[r]] = Eta[[r]][hM$Pi[,r],]%*%Lambda[[r]]
   }
   if(hM$nr > 0){
      Z = LFix + Reduce("+", LRan)
   } else
      Z = LFix

   Z = updateZ(Y=hM$Y,Z=Z,Beta=Beta,iSigma=sigma^-1,Eta=Eta,Lambda=Lambda, X=hM$X,Pi=hM$Pi,distr=hM$distr)

   parList$Gamma = Gamma
   parList$V = V
   parList$Beta = Beta
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

