# updating Gamma after marginalizing Beta

updateGamma2 = function(Z,Gamma,iV,iSigma, X,Pi,Tr,C, iQg, mGamma,iUGamma){
   ns = ncol(Z)
   nc = nrow(iV)
   nt = ncol(Tr)

   inv = function(A){
      return(chol2inv(chol(A)))
   }
   cholL = function(A){
      return(t(chol(A)))
   }

   LFix = 0
   LRan = vector("list", nr)
   for(r in 1:nr){
      LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
   }
   S = Z - (LFix + Reduce("+", LRan))

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

   }
   return(Gamma)
}


