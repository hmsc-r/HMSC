combineParameters = function(Beta, Gamma, iV, rho, iSigma, Eta,
   Lambda,Alpha,Psi,Delta, nc, XScalePar, XInterceptInd, nt, TrScalePar, TrInterceptInd, rhopw){
   for(p in 1:nt){
      m = TrScalePar[1,p]
      s = TrScalePar[2,p]
      if(m!=0 || s!=1){
         Gamma[,p] = Gamma[,p]/s
         if(!is.null(TrInterceptInd)){
            Gamma[,TrInterceptInd] = Gamma[,TrInterceptInd] - m*Gamma[,p]
         }
      }
   }

   for(k in 1:nc){
      m = XScalePar[1,k]
      s = XScalePar[2,k]
      if(m!=0 || s!=1){
         Beta[k,] = Beta[k,]/s
         Gamma[k,] = Gamma[k,]/s
         if(!is.null(XInterceptInd)){
            Beta[XInterceptInd,] = Beta[XInterceptInd,] - m*Beta[k,]
            Gamma[XInterceptInd,] = Gamma[XInterceptInd,] - m*Gamma[k,]
         }
         iV[k,] = iV[k,]*s
         iV[,k] = iV[,k]*s
      }
   }

   V = chol2inv(chol(iV))
   sigma = 1/iSigma
   par = list(Beta=Beta, Gamma=Gamma, V=V, rho=rhopw[rho,1], sigma=sigma, Eta=Eta, Lambda=Lambda, Alpha=Alpha, Psi=Psi, Delta=Delta)
}


