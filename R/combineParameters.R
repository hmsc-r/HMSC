combineParameters = function(Beta,Gamma,iV,rho,iSigma,Eta,Lambda,Alpha,Psi,Delta, rhopw){

   for(k in 1:self$nc){
      if(any(self$XScalePar[,k] != 0)){
         m = self$XScalePar[1,k]
         s = self$XScalePar[2,k]
         Beta[k,] = Beta[k,]/s
         Gamma[k,] = Gamma[k,]/s
         if(!is.null(self$XInterceptInd)){
            Beta[self$XInterceptInd,] = Beta[self$XInterceptInd,] - m*Beta[k,]
            Gamma[self$XInterceptInd,] = Gamma[self$XInterceptInd,] - m*Gamma[k,]
         }
         iV[k,] = iV[k,]*s
         iV[,k] = iV[,k]*s
      }
   }


   V = chol2inv(chol(iV))
   sigma = 1/iSigma
   par = list(Beta=Beta, Gamma=Gamma, V=V, rho=rhopw[rho,1], sigma=sigma, Eta=Eta, Lambda=Lambda, Alpha=Alpha, Psi=Psi, Delta=Delta)
}

Hmsc$set("private", "combineParameters", combineParameters, overwrite=TRUE)


