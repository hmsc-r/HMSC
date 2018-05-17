updateZ = function(Y,Z,Beta,iSigma,Eta,Lambda, X,Pi,distr, ind){
   ZPrev = Z
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

   indColPoisson = (distr[,1]==3)
   pN = sum(indColPoisson)
   if(pN > 0){
      r = 1000 # acquiring Poisson as limit of negative-binomial
      ZPoisson = matrix(NA,ny,pN)
      YPoisson = Y[,indColPoisson]
      EPoisson = E[,indColPoisson]
      stdPoisson = std[,indColPoisson]
      indCellPoisson = !indNA[,indColPoisson]
      y = YPoisson[indCellPoisson]
      e = EPoisson[indCellPoisson]
      s = stdPoisson[indCellPoisson]
      zPrev = ZPrev[,indColPoisson][indCellPoisson]
      w = rpg(length(y), r+y, zPrev)
      prec = s^-2
      sigmaZ = (prec + w)^-1
      muZ = sigmaZ*((y-r)/2 + prec*e)
      z = rnorm(length(y), muZ, sqrt(sigmaZ))
      if(any(is.na(z) | is.nan(z))){
         print("Fail in Poisson Z update")
      }
      ZPoisson[indCellPoisson] = z
      Z[,indColPoisson] = ZPoisson
   }

   Z[indNA] = rnorm(sum(indNA), E[indNA], std[indNA])
   return(Z)
}


