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
      Gamma = t(rmvnorm(self$nt, self$mGamma, self$UGamma))
   }

   if(!is.null(initPar$Gamma)){
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
         Beta[,j] = rmvnorm(1, Mu[,j], self$UGamma)
   }

   if(!is.null(initPar$sigma)){
      sigma = initPar$sigma
   } else{
      sigma = rep(NA, ns)
      indVarFix = (self$distr[,2] == 0)
      sigma[indVarFix] = 1
      sigma[!indVarFix] = rgamma(sum(!indVarFix), shape=self$aSigma[!indVarFix], rate=self$bSigma[!indVarFix])
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

   for(r in 1:self$nr){
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
         Delta[[r]] = matrix(c(rgamma(1,self$a1[r],self$b1[r]), rgamma(nf[r]-1,self$a2[r],self$b2[r])))
      }
      if(is.null(initPar$Psi[[r]])){
         Psi[[r]] = matrix(rgamma(nf[r]*ns, self$nu[r]/2, self$nu[r]/2), nf[r], ns)
      }
      if(is.null(initPar$Lambda[[r]])){
         tau = apply(Delta[[r]], 2, cumprod)
         tauMat = matrix(tau,nf[r],ns)
         mult = sqrt(Psi[[r]]*tauMat)^-1
         Lambda[[r]] = matrix(rnorm(nf[r]*ns)*mult, nf[r], ns)
      }
   }

   Eta = vector("list", self$nr)
   np = self$np

   if(!is.null(initPar$Eta)){
      Eta = initPar$Eta
   } else{
      Eta = vector("list", self$nr)
   }
   for(r in 1:self$nr){
      if(!is.null(initPar$Eta[[r]])){
         Eta[[r]] = initPar$Eta[[r]]
      } else{
         Eta[[r]] = matrix(rnorm(np[r]*nf[r]),np[r],nf[r])
      }
   }


   Z = self$Y

   parList$Gamma = Gamma
   parList$V = V
   parList$Beta = Beta
   parList$sigma = sigma
   parList$Eta = Eta
   parList$Lambda = Lambda
   parList$Psi = Psi
   parList$Delta = Delta
   parList$Z = Z

   return(parList)
}

Hmsc$set("private", "computeInitialParameters", computeInitialParameters, overwrite=TRUE)

