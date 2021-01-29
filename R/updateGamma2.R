
### updating Gamma after marginalizing Beta

#' @importFrom stats rnorm
#'
updateGamma2 = function(Z,Gamma=Gamma,iV,iSigma,Eta,Lambda, X,Pi,dfPi,Tr,C,rL, iQg, mGamma,iUGamma){
   ns = ncol(Z)
   ny = nrow(Z)
   nc = nrow(iV)
   nt = ncol(Tr)
   nr = ncol(Pi)

   # inv = function(A){
   #    return(chol2inv(chol(A)))
   # }
   # cholL = function(A){
   #    return(t(chol(A)))
   # }

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
      # if(all(iSigma==1)){
      #    iV0 = iUGamma[1:nc,1:nc]
      #    V0 = inv(iV0)
      #    XtX = crossprod(X)
      #    TT = crossprod(Tr)
      #    iP = inv(iV + XtX)
      #    LiP = cholL(iP)
      #
      #    R = inv( kronecker(diag(nt),iV0) + kronecker(TT,iV-tcrossprod(iV%*%LiP)) )
      #    LR = cholL(R)
      #    XZT = crossprod(X,S%*%Tr)
      #    iPXZT = iP%*%XZT
      #    tmp = kronecker(TT, V0%*%XtX%*%iP%*%iV)
      #    muG = as.vector(V0%*%(XZT - XtX%*%iPXZT)) - tmp %*% R %*% as.vector(iV%*%iPXZT)
      #    SigmaG = kronecker(diag(nt),V0) - kronecker(TT,tcrossprod(V0%*%t(X))-tcrossprod(V0%*%XtX%*%LiP)) + tcrossprod(tmp%*%LR)
      #    LSigmaG = cholL(SigmaG)
      #    Gamma = muG + LSigmaG%*%rnorm(nc*nt)
      #    Gamma = matrix(Gamma,nc,nt)
      # }

      XtX = crossprod(X)
      if(all(iSigma==1)){
         TtT = crossprod(Tr)
         B = iV + XtX
         RB = chol(B)
         iSg = iUGamma + kronecker(TtT,XtX-crossprod(backsolve(RB,XtX,transpose=TRUE)))
         XtST = crossprod(X,S%*%Tr)
         iBXtST = backsolve(RB,backsolve(RB,XtST,transpose=TRUE))
         XtXiBXtST = XtX%*%iBXtST
         mg0 = iUGamma%*%mGamma + as.vector(XtST - XtXiBXtST)
      } else{
         iDT = matrix(iSigma,ns,nt)*Tr
         iD05T = matrix(sqrt(iSigma),ns,nt)*Tr
         TtiDT = crossprod(iD05T)
         Bst = array(rep(iV,each=ns),c(ns,nc,nc)) + array(iSigma,c(ns,nc,nc))*array(rep(XtX,each=ns),c(ns,nc,nc))
         XtSiD = matrix(iSigma,nc,ns,byrow=TRUE)*crossprod(X,S)

         RBst = array(NA, c(ns,nc,nc))
         iLBstXtX = array(NA, c(ns,nc,nc))
         XtXiBstXtX = array(NA, c(ns,nc,nc))
         tmp1 = matrix(0,nt*nc,nt*nc)
         iBXtSiD = matrix(NA,nc,ns)
         for(j in 1:ns){ # TODO this cycle shall be redone as batched operations
            RBst[j,,] = chol(Bst[j,,])
            iLBstXtX[j,,] = backsolve(RBst[j,,],XtX,transpose=TRUE)
            XtXiBstXtX[j,,] = crossprod(iLBstXtX[j,,])
            tmp1 = tmp1 + iSigma[j]^2 * kronecker(crossprod(Tr[j,,drop=FALSE]), XtXiBstXtX[j,,])
            iBXtSiD[,j] = backsolve(RBst[j,,],backsolve(RBst[j,,],XtSiD[,j],transpose=TRUE))
         }
         iSg = iUGamma + kronecker(TtiDT,XtX) - tmp1
         tmp2 = XtX %*% ((matrix(iSigma,nc,ns,byrow=TRUE)*iBXtSiD) %*% Tr)
         mg0 = iUGamma%*%mGamma + as.vector(XtSiD%*%Tr - tmp2)
      }
      RiSg = chol(iSg)
      mg1 = backsolve(RiSg, mg0, transpose=TRUE)
      Gamma = matrix(backsolve(RiSg,mg1+rnorm(nc*nt)),nc,nt)
   } else{
      # to be implemented later
   }
   return(Gamma)
}


