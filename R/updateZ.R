#' @importFrom stats dnorm pnorm rnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom BayesLogit rpg
#' @importFrom tensorflow tf_probability tf_function
updateZ = function(Y,Z,Beta,iSigma,Eta,Lambda, X,Pi,dfPi,distr,rL, ind,TensorFlowAccelerationFlag=FALSE){
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
      LRan[[r]] = computePredictor.HmscRandomLevel(rL[[r]], Eta[[r]], Lambda[[r]], Pi[,r], dfPi[,r])
   }
   if(nr > 0){
      E = LFix + Reduce("+", LRan)
   } else
      E = LFix

   if(!exists("updateZ_NA", envir = parent.frame())) {
      indNA = is.na(Y)
      lenNA = sum(indNA)
      tmpList = list(indNA=indNA, lenNA=lenNA)
      assign("updateZ_NA", tmpList, envir=parent.frame())
   } else{
      tmpList = get("updateZ_NA", envir=parent.frame())
      indNA = tmpList$indNA
      lenNA = tmpList$lenNA
   }
   Z = matrix(NA,ny,ns)
   std = matrix(iSigma^-0.5,ny,ns,byrow=TRUE)

   if(!exists("updateZ_normal", envir = parent.frame())) {
      indColNormal = (distr[,1]==1)
      nN = sum(indColNormal)
      tmpList = list(indColNormal=indColNormal, nN=nN)
      assign("updateZ_normal", tmpList, envir=parent.frame())
   } else {
      tmpList = get("updateZ_normal", envir=parent.frame())
      indColNormal = tmpList$indColNormal
      nN = tmpList$nN
   }
   if(nN > 0){
      if(nN < ns){
         Z[,indColNormal] = Y[,indColNormal]
      } else
         Z = Y
   }

   if(!exists("updateZ_probit_outer", envir = parent.frame())) {
      indColProbit = (distr[,1]==2)
      pN = sum(indColProbit)
      tmpList = list(indColProbit=indColProbit, pN=pN)
      assign("updateZ_probit_outer", tmpList, envir=parent.frame())
   } else {
      tmpList = get("updateZ_probit_outer", envir=parent.frame())
      indColProbit = tmpList$indColProbit
      pN = tmpList$pN
   }
   if(pN > 0){
      ZProbit = matrix(NA,ny,pN)
      if(!exists("updateZ_probit_inner", envir = parent.frame())) {
         YProbit = Y[,indColProbit]
         indCellProbit = !indNA[,indColProbit]
         YProbit = as.logical(YProbit[indCellProbit])
         lB = rep(-Inf, length(YProbit))
         uB = rep(Inf, length(YProbit))
         lB[YProbit] = 0
         uB[!YProbit] = 0
         anyIndCellProbit = any(indCellProbit)
         lengthYProbit = length(YProbit)
         tmpList = list(indCellProbit=indCellProbit, lB=lB, uB=uB, anyIndCellProbit=anyIndCellProbit, lengthYProbit=lengthYProbit)
         assign("updateZ_probit_inner", tmpList, envir=parent.frame())
      } else {
         tmpList = get("updateZ_probit_inner", envir=parent.frame())
         indCellProbit = tmpList$indCellProbit
         lB = tmpList$lB
         uB = tmpList$uB
         anyIndCellProbit = tmpList$anyIndCellProbit
         lengthYProbit = tmpList$lengthYProbit
      }
      if(anyIndCellProbit){
         if(pN < ns){
            EProbit = E[,indColProbit]
            stdProbit = std[,indColProbit]
         } else{
            EProbit = E
            stdProbit = std
         }
         if(lengthYProbit < pN*ny){
            e = EProbit[indCellProbit]
            s = stdProbit[indCellProbit]
         } else{
            e = as.vector(EProbit)
            s = as.vector(stdProbit)
         }
         if(TensorFlowAccelerationFlag==FALSE){
            z = rtruncnorm(lengthYProbit, a=lB, b=uB, mean=e, sd=s) # this is often the bottleneck for performance
         } else{
            if(!exists("updateZ_probit_tf", envir = parent.frame())){
               f = function(e,s, lB,uB){
                  tfp = tf_probability()
                  ZNew = tfp$distributions$TruncatedNormal(e,s, lB,uB, name='TruncatedNormal')$sample()
                  return(ZNew)
               }
               assign("updateZ_probit_tf", tf_function(f), envir = parent.frame())
            } else{
               updateZ_probit_tf = get("updateZ_probit_tf", envir=parent.frame())
            }
            ZNew = updateZ_probit_tf(tf$constant(e),tf$constant(s), tf$constant(lB),tf$constant(uB))
            z = ZNew$numpy()
         }
         if(lengthYProbit < pN*ny){
            ZProbit[indCellProbit] = z
         } else
            ZProbit = matrix(z,ny,pN)
      }
      if(pN < ns){
         Z[,indColProbit] = ZProbit
      } else
         Z = ZProbit
   }


   if(!exists("updateZ_poisson_outer", envir = parent.frame())) {
      indColPoisson = (distr[,1]==3)
      pN = sum(indColPoisson)
      tmpList = list(indColPoisson=indColPoisson, pN=pN)
      assign("updateZ_poisson_outer", tmpList, envir=parent.frame())
   } else {
      tmpList = get("updateZ_poisson_outer", envir=parent.frame())
      indColPoisson = tmpList$indColPoisson
      pN = tmpList$pN
   }
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
         w = rpg(num=length(y), h=y+r, z=zPrev-log(r))
         prec = s^-2
         sigmaZ = (prec + w)^-1
         muZ = sigmaZ*((y-r)/(2) + prec*(e-log(r))) + log(r)
         z = rnorm(length(y), muZ, sqrt(sigmaZ))
         if(any(is.na(z) | is.nan(z))){
            warning("Fail in Poisson Z update")
         }
         ZPoisson[indCellPoisson] = z
         Z[,indColPoisson] = ZPoisson
      }
   }

   if(lenNA > 0)
      Z[indNA] = rnorm(sum(indNA), E[indNA], std[indNA])
   return(Z)
}
