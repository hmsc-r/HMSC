updateGammaV = function(Beta,Gamma,iV, Tr, mGamma,iUGamma,V0,f0){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)

   iQ = diag(ns)

   Mu = tcrossprod(Gamma, Tr)
   E = Beta - Mu
   A = E %*% tcrossprod(iQ, E)
   Vn = chol2inv(chol(A+V0))
   iV = rwish(f0+ns, Vn);

   tmp = crossprod(Tr, iQ)
   R = chol(Matrix(iUGamma + kronecker(iV, tmp %*% Tr)))
   res = tmp %*% crossprod(Beta, iV)
   mGammaS = chol2inv(R) %*% (iUGamma %*% mGamma + as.vector(res))
   Gamma = mGammaS + solve(R) %*% rnorm(nc*nt)
   return(list(Gamma=Gamma, iV=iV))
}


