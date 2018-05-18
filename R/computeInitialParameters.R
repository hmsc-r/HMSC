#' @title computeInitialParameters
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'



computeInitialParameters = function(initPar){
   parList = list()

   if(!is.null(initPar$Gamma)){
      Gamma = initPar$Gamma
   } else{
      Gamma = matrix(rmvnorm(1, self$mGamma, self$UGamma), self$nc, self$nt)
   }

   if(!is.null(initPar$V)){
      V = initPar$V
   } else{
      V = riwish(self$f0, self$V0)
   }

   if(!is.null(initPar$Beta)){
      Beta = initPar$Beta
   } else{
      Beta = matrix(NA, self$nc, self$ns)
      Mu = tcrossprod(Gamma,self$Tr)
      for(j in 1:self$ns)
         Beta[,j] = rmvnorm(1, Mu[,j], V)
   }

   if(!is.null(initPar$sigma)){
      sigma = initPar$sigma
   } else{
      sigma = rep(NA, self$ns)
      for(j in 1:self$ns){
         if(self$distr[j,2] == 1){
            sigma[j] = rgamma(1, shape=self$aSigma[j], rate=self$bSigma[j])
         } else{
            switch(self$distr[j,1],
               sigma[j] <- 1,
               sigma[j] <- 1,
               sigma[j] <- 1e-3
            )
         }
      }
   }

   nf = rep(NA, self$nr)
   Delta = vector("list", self$nr)
   Psi = vector("list", self$nr)
   Lambda = vector("list", self$nr)
   if(!is.null(initPar$Delta)){
      Delta = initPar$Delta
   } else{
      Delta = vector("list", self$nr)
   }
   if(!is.null(initPar$Psi)){
      Psi = initPar$Psi
   } else{
      Psi = vector("list", self$nr)
   }
   if(!is.null(initPar$Lambda)){
      Lambda = initPar$Lambda
   } else{
      Lambda = vector("list", self$nr)
   }

   for(r in seq_len(self$nr)){
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
         Delta[[r]] = matrix(c(rgamma(1,self$rL[[r]]$a1,self$rL[[r]]$b1), rgamma(nf[r]-1,self$rL[[r]]$a2,self$rL[[r]]$b2)))
      }
      if(is.null(initPar$Psi[[r]])){
         Psi[[r]] = matrix(rgamma(nf[r]*self$ns, self$rL[[r]]$nu/2, self$rL[[r]]$nu/2), nf[r], self$ns)
      }
      if(is.null(initPar$Lambda[[r]])){
         tau = apply(Delta[[r]], 2, cumprod)
         tauMat = matrix(tau,nf[r],self$ns)
         mult = sqrt(Psi[[r]]*tauMat)^-1
         Lambda[[r]] = matrix(rnorm(nf[r]*self$ns)*mult, nf[r], self$ns)
      }
   }

   Eta = vector("list", self$nr)
   np = self$np

   if(!is.null(initPar$Eta)){
      Eta = initPar$Eta
   } else{
      Eta = vector("list", self$nr)
   }
   for(r in seq_len(self$nr)){
      if(!is.null(initPar$Eta[[r]])){
         Eta[[r]] = initPar$Eta[[r]]
      } else{
         Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
      }
   }
   if(!is.null(initPar$Eta)){
      Alpha = initPar$Alpha
   } else{
      Alpha = vector("list", self$nr)
   }
   for(r in seq_len(self$nr)){
      if(!is.null(initPar$Alpha[[r]])){
         Alpha[[r]] = initPar$Alpha[[r]]
      } else{
         Alpha[[r]] = rep(1,nf[r])
      }
   }

   if(!is.null(initPar$rho)){
      rho = which.min(abs(initPar$rho-self$rhopw[,1]))
   } else{
      rho = 1
   }

   LFix = self$X%*%Beta
   LRan = vector("list", self$nr)
   for(r in 1:self$nr){
      LRan[[r]] = Eta[[r]][self$Pi[,r],]%*%Lambda[[r]]
   }
   Z = LFix + Reduce("+", LRan)
   # Z = matrix(0,self$ny,self$ns)
   Z = updateZ(Y=self$Y,Z=Z,Beta=Beta,iSigma=sigma^-1,Eta=Eta,Lambda=Lambda, X=self$X,Pi=self$Pi,distr=self$distr)

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

Hmsc$set("private", "computeInitialParameters", computeInitialParameters, overwrite=TRUE)

