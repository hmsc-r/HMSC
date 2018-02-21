updateZ = function(Y,Z,Beta,iSigma,Eta,Lambda, X,Pi,distr){
   if(any(distr[,1]!=1)){
      ny = nrow(Z)
      ns = ncol(Z)
      nr = ncol(Pi)
      np = apply(Pi, 2, function(a) length(unique(a)))

      LFix = X%*%Beta
      LRan = vector("list", nr)
      for(r in 1:nr){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
      }
      E = LFix + Reduce("+", LRan)

      Z = matrix(NA,ny,ns)

      indNormal = (distr[,1]==1)
      Z[,indNormal] = Y[,indNormal]

      indProbit = (distr[,1]==2)
      pN = sum(indProbit)
      lB = matrix(-Inf,ny,pN)
      uB = matrix(Inf,ny,pN)
      lB[as.logical(Y[,indProbit])] = 0
      uB[!Y[,indProbit]] = 0
      sigma = 1/iSigma
      std = matrix(sqrt(sigma[indProbit]),ny,pN,byrow=TRUE)
      Z[,indProbit] = rtruncnorm(ny*pN, a=lB, b=uB, mean=E[,indProbit], sd=std)
      return(Z)
   } else{
      return(Y)
   }
}


