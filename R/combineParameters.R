combineParameters = function(Beta,Gamma,iV,rho,iSigma,Eta,Lambda,Alpha,Psi,Delta, rhopw){

   for(p in 1:self$nt){
      m = self$TrScalePar[1,p]
      s = self$TrScalePar[2,p]
      if(m!=0 || s!=1){
         Gamma[,p] = Gamma[,p]/s
         if(!is.null(self$TrInterceptInd)){
            Gamma[,self$TrInterceptInd] = Gamma[,self$TrInterceptInd] - m*Gamma[,p]
         }
      }
   }

   for(k in 1:self$nc){
      m = self$XScalePar[1,k]
      s = self$XScalePar[2,k]
      if(m!=0 || s!=1){
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


