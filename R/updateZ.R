updateZ = function(Y,Beta,iSigma,Eta,Lambda, X,Pi,distr){
   if(any(distr[,1]!=1)){
      ny = nrow(Y)
      ns = ncol(Y)
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
      # Z[,indProbit] = rtnorm(ny*pN, mean=E[,indProbit], sd=1, lower=lB,upper=uB)
      # Z[Z>30] = 30
      # Z[Z<-30] = -30
      # plot(E,Z)
      # a=1
      return(Z)
   } else{
      return(Y)
   }
}


