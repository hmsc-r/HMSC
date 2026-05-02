updateRho = function(Beta,Gamma,iV, RQg,detQg, Tr, rhopw){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)
   rhoN = nrow(rhopw)

   Mu = tcrossprod(Gamma, Tr)
   E = t(Beta - Mu)

   RiV = chol(iV)
   E = tcrossprod(E,RiV);
   v = rep(NA, rhoN)
   for(rN in 1:rhoN){
      v[rN] = sum(backsolve(RQg[,,rN],E,transpose=TRUE)^2)
   }
   # logdetg = -ns*log(det(iV))+nc*detQg
   logdetg = nc*detQg
   logLike = log(rhopw[,2]) - 0.5*logdetg - 0.5*v;
   logLike = logLike - max(logLike)
   like = exp(logLike)
   rhoInd = sample.int(rhoN, 1, prob=like)
   return(rhoInd)
}


