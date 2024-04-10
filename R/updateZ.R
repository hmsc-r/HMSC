#' @importFrom stats dnorm pnorm rnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom BayesLogit rpg
updateZ = function(Y,Z,Beta,iSigma,Eta,Lambda, Loff,X,Pi,dfPi,distr,rL, ind){
   ZPrev = Z
   ny = nrow(Y)
   ns = ncol(Y)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))

   switch(class(X)[1L],
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
            LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),k]) %*% Lambda[[r]][,,k]
      }
   }
   E = Reduce("+", c(list(LFix), LRan))
   if(!is.null(Loff)) E = E + Loff

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
      if(any(indCellProbit)){
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
   }

   indColPoisson = (distr[,1]==3)
   pN = sum(indColPoisson)
   if(pN > 0){
      r = 1e3# acquiring Poisson as limit of negative-binomial
      ZPoisson = matrix(NA,ny,pN)
      YPoisson = Y[,indColPoisson]
      EPoisson = E[,indColPoisson]
      stdPoisson = std[,indColPoisson]
      indCellPoisson = !indNA[,indColPoisson]
      if(any(indCellPoisson)){
         y = YPoisson[indCellPoisson]
         e = EPoisson[indCellPoisson]
         s = stdPoisson[indCellPoisson]
         zPrev = ZPrev[,indColPoisson][indCellPoisson]
         w = rpg(num=length(y), h=y+r, z=zPrev-1*log(r))
         prec = s^-2
         sigmaZ = (prec + w)^-1
         muZ = sigmaZ*((y-r)/(2) + prec*(e-log(r))) + 1*log(r)
         z = rnorm(length(y), muZ, sqrt(sigmaZ))
         if(any(is.na(z) | is.nan(z))){
            warning("failure in Poisson Z update")
         }
         ZPoisson[indCellPoisson] = z
         Z[,indColPoisson] = ZPoisson
      }
   }

   Z[indNA] = rnorm(sum(indNA), E[indNA], std[indNA])
   return(Z)
}
