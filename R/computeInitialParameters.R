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
      sigma = rgamma(self$ns, shape=self$aSigma, rate=self$bSigma)
   }

   Z = self$Y

   parList$Gamma = Gamma
   parList$V = V
   parList$Beta = Beta
   parList$sigma = sigma
   parList$Z = Z

   return(parList)
}

Hmsc$set("private", "computeInitialParameters", computeInitialParameters, overwrite=TRUE)

