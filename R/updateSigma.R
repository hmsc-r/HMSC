updateInvSigma = function(Y,Z,Beta,iSigma,Eta,Lambda, distr,X,Pi, aSigma,bSigma){
   indVarSigma = (distr[,2]==1)
   if(any(indVarSigma)){
      ny = nrow(Z)
      ns = ncol(Z)
      nr = ncol(Pi)

      switch(X,
         matrix = {
            LFix = X%*%Beta
         },
         list = {
            LFix = matrix(NA,ny,ns)
            for(j in 1:ns)
               LFix[,j] = X[[j]]%*%Beta[,j]
         }
      )
      LRan = vector("list", nr)
      for(r in seq_len(nr)){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
      }
      if(nr > 0){
         Eps = Z - (LFix + Reduce("+", LRan))
      } else
         Eps = Z - LFix

      Yx = !is.na(Y)
      nyx = colSums(Yx)
      shape = aSigma + nyx/2
      rate = bSigma + apply((Eps*Yx)^2, 2, sum)/2

      iSigma[indVarSigma] = rgamma(sum(indVarSigma), shape[indVarSigma], rate[indVarSigma])
   }
   return(iSigma)
}

