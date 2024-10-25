#' @importFrom stats rnorm
#' @importFrom MCMCpack rwish
#'
updateGammaV = function(Beta,Gamma,iV,rhoInd, Tr,phyloFlag,phyloFast,phyloTreeList,phyloTreeRoot, iQg,RQg, mGamma,iUGamma,V0,f0,rhopw){
   ns = ncol(Beta)
   nc = nrow(Beta)
   nt = ncol(Tr)

   Mu = tcrossprod(Gamma, Tr)
   E = Beta - Mu
   if(phyloFlag == FALSE){
      A = tcrossprod(E)
   } else{
      if(phyloFast == FALSE){
         A = E %*% tcrossprod(iQg[,,rhoInd], E)
      } else{
         E_arr = array(E, c(nc,ns,1))
         res = fastPhyloBilinearDet(phyloTreeList, E_arr, E_arr, phyloTreeRoot, 1, rhopw[rhoInd,1])
         A = res$XiSY
      }
   }
   Vn = chol2inv(chol(A+V0))
   iV = rwish(f0+ns, Vn)

   if(phyloFlag == FALSE){
      R = chol(iUGamma + kronecker(crossprod(Tr), iV))
      mg = chol2inv(R) %*% (iUGamma%*%mGamma + as.vector((iV%*%Beta)%*%Tr))
   } else{
      if(phyloFast == FALSE){
         TrT_iQ_Tr = crossprod(backsolve(RQg[,,rhoInd], Tr, transpose=TRUE))
         iVBeta_iQ_Tr = (iV%*%Beta) %*% (iQg[,,rhoInd]%*%Tr)
      } else{
         Tr_arr = array(t(Tr), c(nt,ns,1))
         TrT_iQ_Tr = fastPhyloBilinearDet(phyloTreeList, Tr_arr, Tr_arr, phyloTreeRoot, 1, rhopw[rhoInd,1])$XiSY
         BetaT_iV_arr = array(iV%*%Beta, c(nc,ns,1))
         iVBeta_iQ_Tr = fastPhyloBilinearDet(phyloTreeList, BetaT_iV_arr, Tr_arr, phyloTreeRoot, 1, rhopw[rhoInd,1])$XiSY
      }
      R = chol(iUGamma + kronecker(TrT_iQ_Tr, iV))
      mg = chol2inv(R) %*% (iUGamma%*%mGamma + as.vector(iVBeta_iQ_Tr))
   }
   Gamma = matrix(mg + backsolve(R,rnorm(nc*nt)),nc,nt)
   return(list(Gamma=Gamma, iV=iV))
}


