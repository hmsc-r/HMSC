#' @importFrom stats rnorm
#' @importFrom tensorflow tf
#'
updateGamma2 = function(Z,Gamma,iV,rho,iD,Eta,Lambda, X,Pi,dfPi,Tr,C,rL, iQ, mGamma,iUGamma, tfCompFlag=FALSE){
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
   if(tfCompFlag==TRUE){
      tfla = tf$linalg
      ic = function(...) as.integer(c(...))
      XtiDX = tfla$einsum("ic,ij,ik->jck",X,iD,X)
      TtkXt_iD_XkT = tf$reshape(tfla$einsum("ja,jck,jb->acbk",Tr,XtiDX,Tr), ic(nt*nc,nt*nc))
      XiDS = tfla$einsum("ik,ij->jk",X,iDS)
   } else{
      XtiDX = array(NA, c(ns,nc,nc))
      for(j in 1:ns){
         XtiDX[j,,] = crossprod(sqrt(iD[,j])*X)
      }
   }
   if(is.null(C)){
      if(tfCompFlag==TRUE){
         Bst = iV + XtiDX
         LBst = tfla$cholesky(Bst)
         iLBstXtiDX = tfla$triangular_solve(LBst, XtiDX)
         XtiDXiBstXtiDX = tf$matmul(iLBstXtiDX,iLBstXtiDX,transpose_a=TRUE)
         tmp1 = tf$reshape(tfla$einsum("jt,jck,jv->tcvk",Tr,XtiDXiBstXtiDX,Tr), ic(nt*nc,nt*nc))
         iSigmaG = iUGamma + TtkXt_iD_XkT - tmp1
         tmp2 = tf$squeeze(tfla$cholesky_solve(LBst, XiDS[,,NULL]), ic(-1))
         tmp3 = tfla$einsum("jt,jck,jk->tc",Tr,XtiDX,tmp2)
         mG0 = as.vector(iUGamma%*%mGamma) + tf$reshape(tfla$einsum("ik,ij,jt->tk",X,iDS,Tr),ic(nc*nt)) - tf$reshape(tmp3, ic(nc*nt))
         iSigmaG = iSigmaG$numpy(); mG0 = mG0$numpy()
      } else{
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
      }
      RiSigmaG = chol(iSigmaG)
      mG1 = backsolve(RiSigmaG, mG0, transpose=TRUE)
      Gamma = matrix(backsolve(RiSigmaG,mG1+rnorm(nc*nt)),nc,nt)
   } else{
      if(tfCompFlag==TRUE){
         W = kronecker(iQ,iV) + tf$reshape(tf$transpose(tfla$diag(tf$transpose(XtiDX,ic(1,2,0))), ic(2,0,3,1)), ic(ns*nc,ns*nc))
         LW = tfla$cholesky(W)
         iW = tfla$cholesky_solve(LW, tf$eye(ic(ns*nc),dtype=tf$float64))
         tmp01 = tfla$einsum("jac,jcgk,gbk->jagb",XtiDX,tf$reshape(iW, ic(ns,nc,ns,nc)),XtiDX)
         tmp02 = tfla$einsum("jt,jcgk,gv->tcvk",Tr,tmp01,Tr)
         iSigmaG = iUGamma + TtkXt_iD_XkT - tf$reshape(tmp02, ic(nt*nc,nt*nc))
         tmp2 = tf$reshape(tfla$cholesky_solve(LW, tf$reshape(XiDS,ic(ns*nc,1))), ic(ns,nc))
         tmp3 = tfla$einsum("jt,jck,jk->tc",Tr,XtiDX,tmp2)
         mG0 = as.vector(iUGamma%*%mGamma) + tf$reshape(tfla$einsum("ik,ij,jt->tk",X,iDS,Tr),ic(nc*nt)) - tf$reshape(tmp3, ic(nc*nt))
         iSigmaG = iSigmaG$numpy(); mG0 = mG0$numpy()
         RiSigmaG = chol(iSigmaG)
         mG1 = backsolve(RiSigmaG, mG0, transpose=TRUE)
         Gamma = matrix(backsolve(RiSigmaG,mG1+rnorm(nc*nt)),nc,nt)
      }else{
         warning("updateGamma2() is not yet implemented for model with phylogeny and tfCompFlag=FALSE. Adjust the tfCompFlag or disable this sampler. ")
      }

   }
   return(Gamma)
}


