#' @importFrom stats rnorm
#'
updateGamma2 = function(Z,Gamma,iV,rho,iD,Eta,Lambda, X,Pi,dfPi,Tr,C,rL, iQ, mGamma,iUGamma){
   ns = ncol(Z)
   ny = nrow(Z)
   nc = nrow(iV)
   nt = ncol(Tr)
   nr = ncol(Pi)

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
   iDS = iD*S
   iDS[is.na(Z)] = 0
   XtiDS = crossprod(X,iDS)
   XtiDX = array(NA, c(ns,nc,nc))
   for(j in 1:ns){
      XtiDX[j,,] = crossprod(sqrt(iD[,j])*X)
   }

   if(is.null(C)){
      Bst = RBst = iLBstXtiDX = XtiDXiBstXtiDX = array(NA, c(ns,nc,nc))
      tmp1 = matrix(0,nc*nt,nc*nt)
      iBXtiDS = XtiDXiBXtiDS = matrix(NA,nc,ns)
      for(j in 1:ns){
         Bst[j,,] = iV + XtiDX[j,,]
         RBst[j,,] = chol(Bst[j,,])
         iLBstXtiDX[j,,] = backsolve(RBst[j,,],XtiDX[j,,],transpose=TRUE)
         XtiDXiBstXtiDX[j,,] = crossprod(iLBstXtiDX[j,,])
         tmp1 = tmp1 + kronecker(crossprod(Tr[j,,drop=FALSE]), XtiDX[j,,]-XtiDXiBstXtiDX[j,,])
         iBXtiDS[,j] = backsolve(RBst[j,,],backsolve(RBst[j,,],XtiDS[,j],transpose=TRUE))
         XtiDXiBXtiDS[,j] = XtiDX[j,,] %*% iBXtiDS[,j]
      }
      iSigmaG = iUGamma + tmp1
      mG0 = as.vector(iUGamma%*%mGamma) + as.vector(XtiDS%*%Tr - XtiDXiBXtiDS%*%Tr)
      RiSigmaG = chol(iSigmaG)
      mG1 = backsolve(RiSigmaG, mG0, transpose=TRUE)
      Gamma = matrix(backsolve(RiSigmaG,mG1+rnorm(nc*nt)),nc,nt)
   } else{
      XtiDXbd = bdiag(lapply(split(XtiDX, rep(1:ns,nc^2)), matrix, nrow=nc))
      W = kronecker(iQ,iV) + XtiDXbd
      RW = chol(W)
      RiW_XtiDX_TkI = solve(t(RW), XtiDXbd) %*% kronecker(Tr,Diagonal(nc))
      TtkXt_iD_TkX = matrix(0,nc*nt,nc*nt)
      for(j in 1:ns){
         TtkXt_iD_TkX = TtkXt_iD_TkX + kronecker(crossprod(Tr[j,,drop=FALSE]), XtiDX[j,,])
      }
      iSigmaG = iUGamma + TtkXt_iD_TkX - Matrix::crossprod(RiW_XtiDX_TkI)
      RiSigmaG = chol(iSigmaG)
      XtiDST = XtiDS %*% Tr
      iWXtiDS = solve(RW, solve(t(RW), as.vector(XtiDS)))
      TtkXt_iD_IkX_iWXtiDS = kronecker(t(Tr),Diagonal(nc)) %*% (XtiDXbd %*% iWXtiDS)
      mG0 = as.vector(iUGamma%*%mGamma) + as.vector(XtiDST) - TtkXt_iD_IkX_iWXtiDS
      mG1 = backsolve(RiSigmaG, mG0, transpose=TRUE)
      Gamma = matrix(backsolve(RiSigmaG,mG1+rnorm(nc*nt)),nc,nt)
   }
   return(Gamma)
}


