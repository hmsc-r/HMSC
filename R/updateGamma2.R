
### updating Gamma after marginalizing Beta

#' @importFrom stats rnorm
#'
updateGamma2 = function(Z,Gamma=Gamma,iV,iSigma,Eta,Lambda, X,Pi,dfPi,Tr,C,rL, iQg, mGamma,iUGamma){
   ns = ncol(Z)
   ny = nrow(Z)
   nc = nrow(iV)
   nt = ncol(Tr)
   nr = ncol(Pi)

   inv = function(A){
      return(chol2inv(chol(A)))
   }
   cholL = function(A){
      return(t(chol(A)))
   }

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
   if(nr > 0){
      S = Z - Reduce("+", LRan)
   } else
      S = Z

   if(is.null(C)){
      if(all(iSigma==1)){
         iV0 = iUGamma[1:nc,1:nc]
         V0 = inv(iV0)
         XX = crossprod(X)
         TT = crossprod(Tr)
         iP = inv(iV + XX)
         LiP = cholL(iP)

         R = inv( kronecker(diag(nt),iV0) + kronecker(TT,iV-tcrossprod(iV%*%LiP)) )
         LR = cholL(R)
         XZT = crossprod(X,S%*%Tr)
         iPXZT = iP%*%XZT
         tmp = kronecker(TT, V0%*%XX%*%iP%*%iV)
         muG = as.vector(V0%*%(XZT - XX%*%iPXZT)) - tmp %*% R %*% as.vector(iV%*%iPXZT)
         SigmaG = kronecker(diag(nt),V0) - kronecker(TT,tcrossprod(V0%*%t(X))-tcrossprod(V0%*%XX%*%LiP)) + tcrossprod(tmp%*%LR)
         # SogmaG = (SigmaG+t(SigmaG))/2
         LSigmaG = cholL(SigmaG)
         Gamma = muG + LSigmaG%*%rnorm(nc*nt)
         Gamma = matrix(Gamma,nc,nt)
      }
   } else{
      # to be implemented later
   }
   return(Gamma)
}


