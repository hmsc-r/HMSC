updateRho = function(rho,Beta,Gamma,iV, phyloPar, Tr, rhopw,rhoAlphaDP){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)
   rhoN = nrow(rhopw)

   Mu = tcrossprod(Gamma, Tr)
   BT = t(Beta - Mu)
   RiV = chol(iV)
   if(rhoAlphaDP==0){
      E = tcrossprod(BT,RiV);
      qF = rep(NA, rhoN)
      VTE = crossprod(phyloPar$VC, E)
      for(rN in 1:rhoN){
         # qF[rN] = sum(backsolve(phyloPar$RQg[,,rN],E,transpose=TRUE)^2)
         qF[rN] = sum(((rhopw[rN,1]*phyloPar$eC+(1-rhopw[rN,1]))^-0.5 * VTE)^2)
      }
      logdetg = nc*phyloPar$detQg # -ns*log(det(iV)) part dropped as constant
      logLike = log(rhopw[,2]) - 0.5*logdetg - 0.5*qF
      logLike = logLike - max(logLike)
      rho = sample(rhoN, 1, prob=exp(logLike))
   } else{
      VCB = crossprod(VC,BT)
      VCBArray = array(VCB,c(ns,nc,rhoN))
      countTable = table(c(rho,1:rhoN)) - 1
      eQ = matrix(eC,ns,rhoN) * matrix(rhopw[,1],ns,rhoN,byrow=TRUE) + (1-matrix(rhopw[,1],ns,rhoN,byrow=TRUE))
      ieQ05 = eQ^(-0.5)
      logdetg = colSums(log(eQ))
      for(k in 1:nc){
         countTable[rho[k]] = countTable[rho[k]] - 1
         if(!is.infinite(rhoAlphaDP)){
            condPriorProb = (rhoAlphaDP*rhopw[,2] + countTable) / (nc-1+rhoAlphaDP)
         } else{
            condPriorProb = rhopw[,2]
         }
         # qF = rep(NA, rhoN)
         # rhoTemp = rho
         # for(rN in 1:rhoN){
         #    rhoTemp[k] = rN
         #    qF[rN] = sum(tcrossprod(RiV,VCB*ieQ05[,rhoTemp])^2)
         # }
         ieQ05Array = array(ieQ05[,rho],c(ns,nc,rhoN))
         ieQ05Array[,k,] = ieQ05
         tmp1 = matrix(aperm(VCBArray*ieQ05Array,c(2,1,3)), nc,ns*rhoN)
         tmp2 = colSums((RiV %*% tmp1)^2)
         qF = colSums(matrix(tmp2,ns,rhoN))
         logLike = log(condPriorProb) - 0.5*logdetg - 0.5*qF
         logLike = logLike - max(logLike)
         rho[k] = sample(rhoN, 1, prob=exp(logLike))
         countTable[rho[k]] = countTable[rho[k]] + 1
      }
   }
   return(rho)
}


