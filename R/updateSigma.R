updateInvSigma = function(Z,Beta,Eta,Lambda, distr,X,Pi, aSigma,bSigma){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)

   LFix = X%*%Beta
   LRan = vector("list", nr)
   for(r in 1:nr){
      LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
   }
   Eps = Z - (LFix + Reduce("+", LRan))

   shape = aSigma + ny/2
   rate = bSigma + apply(Eps^2, 2, sum)/2

   iSigma = rep(1,ns)
   indVarSigma = (distr[,2]==1)
   iSigma[indVarSigma] = rgamma(sum(indVarSigma), shape[indVarSigma], rate[indVarSigma])

   return(iSigma)
}

