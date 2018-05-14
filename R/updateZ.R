updateZ = function(Y,Beta,iSigma,Eta,Lambda, X,Pi,distr, ind){
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
   indNA = is.na(Y)
   std = matrix(iSigma^-0.5,ny,ns,byrow=TRUE)

   indColNormal = (distr[,1]==1)
   Z[,indColNormal] = Y[,indColNormal]

   indColProbit = (distr[,1]==2)
   pN = sum(indColProbit)
   if(pN > 0){
      ZProbit = matrix(NA,ny,pN)
      YProbit = Y[,indColProbit]
      EProbit = E[,indColProbit]
      stdProbit = std[,indColProbit]
      indCellProbit = !indNA[,indColProbit]
      YProbit = as.logical(YProbit[indCellProbit])
      e = EProbit[indCellProbit]
      s = stdProbit[indCellProbit]
      lB = rep(-Inf, length(YProbit))
      uB = rep(Inf, length(YProbit))
      lB[YProbit] = 0
      uB[!YProbit] = 0
      z = rtruncnorm(length(YProbit), a=lB, b=uB, mean=e, sd=s) # this is often the bottleneck for performance
      ZProbit[indCellProbit] = z
      Z[,indColProbit] = ZProbit
   }

   Z[indNA] = rnorm(E[indNA], std[indNA])
   return(Z)
}


