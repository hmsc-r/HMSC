updateInvSigma = function(Y,Z,Beta,iSigma,Eta,Lambda, distr,X,Pi,rL, aSigma,bSigma){
   indVarSigma = (distr[,2]==1)
   if(any(indVarSigma)){
      ny = nrow(Z)
      ns = ncol(Z)
      nr = ncol(Pi)

      switch(class(X),
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
         if(rL[[r]]$xDim == 0){
            LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
         } else{
            LRan[[r]] = matrix(0,ny,ns)
            for(k in 1:rL[[r]]$xDim)
               LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[Pi[,r],k]) %*% Lambda[[r]][,,k]
         }
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
