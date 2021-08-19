#' @importFrom stats rnorm
#' @importFrom MCMCpack rwish
#'
# @importFrom tensorflow tf

updateGammaV = function(Beta,Gamma,iV,rho, Tr,C, VC,eC,iQg,RQg, mGamma,iUGamma,V0,f0,rhopw){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)
   Mu = tcrossprod(Gamma, Tr)
   E = Beta - Mu

   if(length(rho)==1){
      if(is.null(C)){
         iQ = diag(ns)
         RQ = diag(ns)
      } else{
         iQ = iQg[,,rho]
         RQ = RQg[,,rho]
      }
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
   } else{
      eQ = matrix(eC,ns,nc) * matrix(rhopw[rho,1],ns,nc,byrow=TRUE) + (1-matrix(rhopw[rho,1],ns,nc,byrow=TRUE))
      ieQ05 = eQ^(-0.5)
      M = ieQ05 * t(E %*% VC)
      A = crossprod(M)
      Vn = chol2inv(chol(A+V0))
      iV = rwish(f0+ns, Vn)
      VCT = crossprod(VC,Tr) # possible to precompute
      # TiQHatT_perm = matrix(tf$linalg$einsum("jt,ja,ab,jb,jq->atbq",VCT,ieQ05,iV,ieQ05,VCT)$numpy(), nt*nc, nt*nc)
      ieQ05VCT = array(ieQ05,c(ns,nc,nt)) * array(VCT[rep(1:ns,nc),],c(ns,nc,nt))
      tmp = matrix(ieQ05VCT,ns,nt*nc)
      TiQHatT_perm = crossprod(tmp) * kronecker(matrix(1,nt,nt),iV)
      R = chol(iUGamma + TiQHatT_perm)
      tmp = as.vector(crossprod(ieQ05*((ieQ05*t(Beta%*%VC))%*%iV),VCT))
      mg = chol2inv(R) %*% (iUGamma%*%mGamma + tmp)
      Gamma = matrix(mg + backsolve(R,rnorm(nc*nt)),nc,nt)
   }
   return(list(Gamma=Gamma, iV=iV))
}


