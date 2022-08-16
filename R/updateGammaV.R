#' @importFrom stats rnorm
#' @importFrom MCMCpack rwish
#'
updateGammaV = function(Beta,Gamma,iV,rho, Tr,C, iQg,RQg, mGamma,iUGamma,V0,f0){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)

   if(is.null(C)){
      iQ = diag(ns)
      RQ = diag(ns)
   } else{
      iQ = iQg[,,rho]
      RQ = RQg[,,rho]
   }

   Mu = tcrossprod(Gamma, Tr)
   E = Beta - Mu
   A = E %*% tcrossprod(iQ, E)
   Vn = chol2inv(chol(A+V0))
   iV = rwish(f0+ns, Vn)

   # tmp = crossprod(Tr, iQ)
   # R = chol(iUGamma + kronecker(iV, tmp %*% Tr))
   # res = tmp %*% crossprod(Beta, iV)
   # mGammaS = chol2inv(R) %*% (iUGamma %*% mGamma + as.vector(res))
   # Gamma = mGammaS + backsolve(R,rnorm(nc*nt))
   # Gamma = matrix(Gamma,nc,nt,byrow=TRUE)

   R = chol(iUGamma + kronecker(crossprod(backsolve(RQ,Tr,transpose=TRUE)), iV))
   mg = chol2inv(R) %*% (iUGamma%*%mGamma + as.vector((iV%*%Beta)%*%(iQ%*%Tr)))
   Gamma = matrix(mg + backsolve(R,rnorm(nc*nt)),nc,nt)
   return(list(Gamma=Gamma, iV=iV))
}


