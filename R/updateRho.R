updateRho = function(Beta,Gamma,iV,rhoInd, Tr, phyloFast,phyloTreeList,phyloTreeRoot, RQg,detQg, rhopw){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)
   rhoGridN = nrow(rhopw)
   rhoLen = length(rhoInd)

   Mu = tcrossprod(Gamma, Tr)
   E = Beta - Mu
   if(rhoLen == 1){
      RiV = chol(iV)
      RiV_E = RiV %*% E
      vg = rep(NA, rhoGridN)
      if(phyloFast == FALSE){
         RiV_E_trans = t(RiV_E)
         for(rN in 1:rhoGridN){
            vg[rN] = sum(backsolve(RQg[,,rN], RiV_E_trans, transpose=TRUE)^2)
         }
         logDetg = nc*detQg
      } else{
         logDetg = rep(NA, rhoGridN)
         RiV_E_arr = array(RiV_E, c(nc,ns,1))
         for(rN in 1:rhoGridN){
            res = fastPhyloBilinearDet(phyloTreeList, RiV_E_arr, RiV_E_arr, phyloTreeRoot, 1, rhopw[rN,1])
            vg[rN] = sum(diag(res$XiSY))
            logDetg[rN] = nc*res$logDet
         }
      }
      logLike = log(rhopw[,2]) - 0.5*logDetg - 0.5*vg
      logLike = logLike - max(logLike)
      like = exp(logLike)
      rhoInd = sample.int(rhoGridN, 1, prob=like)
   } else{
      stop("TODO vector rho case to be finished later")
      # rhoNew = rho
      # for(rhoInd in 1:rhoLen){
      #    rhoVec = rhoNew
      #    like = rep(NA, rhoN)
      #    for(rN in 1:rhoN){
      #       rhoVec[rhoInd] = rN
      #       rho = rhoVec[covRhoInd]
      #       V1 = diag(sqrt(rho)) %*% V %*% diag(sqrt(rho))
      #       V2 = diag(sqrt(1-rho)) %*% V %*% diag(sqrt(1-rho))
      #       QHat = kronecker(C,V1) + kronecker(diag(ns),V2)
      #       RQHat = chol(QHat)
      #       v = sum(backsolve(RQHat, as.vector(E), transpose=TRUE)^2)
      #       logDetQHat = 2*sum(log(diag(RQHat)))
      #       like[rN] = log(rhopw[rN,2]) - 0.5*logDetQHat - 0.5*v
      #    }
      #    rhoNew[rhoInd] = sample.int(rhoN, 1, prob=like)
      # }
      # rho = rhoNew



      # rhoVec = rep(rhopw[rhoInd,1], nc)
      # E_arr = array(E, c(nc,ns,1))
      # res = fastPhyloBilinearDet(phyloTreeList, E_arr, E_arr, phyloTreeRoot, iV, rhoVec)
      # vg[rN] = res$XiSY
      # logDetg[rN] = res$logDet
   }

   return(rhoInd)
}


