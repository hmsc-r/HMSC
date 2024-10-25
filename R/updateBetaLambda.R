# @title updateBetaLambda
#
# @description updates beta lambda
#
#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix
#'
updateBetaLambda = function(Y,Z,Gamma,iV,iSigma,Eta,Psi,Delta,rhoInd, Loff,X,Tr,Pi,dfPi, phyloFlag,phyloFast,phyloTreeList,phyloTreeRoot,covRhoGroup, iQg,rhopw, rL, sdMult=1){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)

   Lambda = vector("list", nr)
   if(nr > 0){
      EtaFull = vector("list", nr)
      nf = rep(NA,nr)
      ncr = rep(NA,nr)
      for(r in seq_len(nr)){
         if(rL[[r]]$xDim == 0){
            EtaFull[[r]] = Eta[[r]][Pi[,r],]
         } else{
            EtaFull[[r]] = vector("list", rL[[r]]$xDim)
            for(k in 1:rL[[r]]$xDim)
               EtaFull[[r]][[k]] = Eta[[r]][Pi[,r],] * rL[[r]]$x[as.character(dfPi[,r]),k]
         }
         nf[r] = ncol(Eta[[r]])
         ncr[r] = max(rL[[r]]$xDim, 1)
      }
      nfSum = sum(nf*ncr)
      EtaSt = matrix(unlist(EtaFull),ny,nfSum)
      if(is.matrix(X)){
         XEta = cbind(X,EtaSt)
      } else if(is.list(X)){
         XEta = lapply(X, cbind, EtaSt)
      }
      PsiT = vector("list",nr)
      for(r in 1:nr){
         if(rL[[r]]$xDim == 0){
            PsiT[[r]] = t(Psi[[r]])
         } else{
            PsiT[[r]] = aperm(Psi[[r]], c(2,1,3))
         }
      }
      psiSt = matrix(unlist(PsiT),nfSum,ns,byrow=TRUE)
      Tau = lapply(Delta, function(a) apply(a,2,cumprod))
      tauSt = matrix(unlist(Tau),nfSum,1)
      priorLambda = psiSt*matrix(tauSt,nfSum,ns)
   } else{
      nf = c()
      ncr = c()
      nfSum = 0
      XEta = X
      priorLambda = matrix(numeric(0),0,ns)
   }

   Mu = tcrossprod(Gamma,Tr)
   if(is.null(Loff)) S = Z else S = Z - Loff
   if(is.matrix(X)){
      M = X %*% Mu
      S = S - M
      XEtaTXEta = crossprod(XEta)
      isXTS = crossprod(XEta,S) * matrix(iSigma,nc+nfSum,ns,byrow=TRUE)
   } else if(is.list(X)){
      M = matrix(NA,ny,ns)
      XEtaTXEta = lapply(XEta, crossprod)
      isXTS = matrix(NA,nc+nfSum,ns)
      for(j in 1:ns){
         M[,j] = X[[j]] %*% Mu[,j]
         S[,j] = S[,j] - M[,j]
         isXTS[,j] = crossprod(XEta[[j]],S[,j]) * iSigma[j]
      }
   }

   if(phyloFlag==FALSE){
      # no phylogeny information
      Yx = !is.na(Y)
      indColFull = apply(Yx,2,all)
      indColNA = !indColFull
      diagiV = diag(iV)
      P0 = matrix(0,nc+nfSum,nc+nfSum)
      P0[1:nc,1:nc] = iV
      BetaLambda = matrix(NA, nc+nfSum, ns)
      #TODO: there is a special case for no latent factors and full data as U can be computed once
      for(j in which(indColFull)){ # test whether worthy to rewrite with tensorA?
         P = P0
         diag(P) = c(diagiV, priorLambda[,j])
         if(is.matrix(X)){
            iU = P + XEtaTXEta*iSigma[j]
         } else if(is.list(X)){
            iU = P + XEtaTXEta[[j]]*iSigma[j]
         }
         RiU = chol(iU)
         U = chol2inv(RiU)
         m = U %*% isXTS[,j]
         BetaLambda[,j] = m + backsolve(RiU, sdMult*rnorm(nc+nfSum))
      }
      for(j in which(indColNA)){
         indObs = Yx[,j]
         if(is.matrix(X)){
            XEtaTXEta = crossprod(XEta[indObs,,drop=FALSE])
            isXTS = crossprod(XEta[indObs,,drop=FALSE],S[indObs,j]) * iSigma[j]
         } else if(is.list(X)){
            XEtaTXEta = crossprod(XEta[[j]][indObs,,drop=FALSE])
            isXTS = crossprod(XEta[[j]][indObs,,drop=FALSE],S[indObs,j]) * iSigma[j]
         }
         P = P0
         diag(P) = c(diagiV, priorLambda[,j])
         iU = P + XEtaTXEta*iSigma[j]
         RiU = chol(iU)
         U = chol2inv(RiU)
         m = U %*% isXTS
         BetaLambda[,j] = m + backsolve(RiU, sdMult*rnorm(nc+nfSum))
      }
   } else{
      # available phylogeny information
      if(phyloFast==FALSE){
         P = bdiag(kronecker(iV,iQg[,,rhoInd]), Diagonal(x=t(priorLambda)))
         if(is.matrix(X)){
            RiU = chol(kronecker(XEtaTXEta,diag(iSigma)) + P)
         } else if(is.list(X)){
            tmp = vector("list",ns)
            for(j in 1:ns)
               tmp[[j]] = XEtaTXEta[[j]] * iSigma[j]
            tmpMat = Reduce(rbind, tmp)
            ind1 = rep(rep(1:ns,each=nc+nfSum)+ns*rep(0:(nc+nfSum-1),ns), nc+nfSum)
            ind2 = rep(1:((nc+nfSum)*ns), each=nc+nfSum)
            mat = sparseMatrix(ind1, ind2, x=as.vector(tmpMat))
            RiU = chol(as.matrix(mat + P))
         }
         m1 = backsolve(RiU, as.vector(t(isXTS)), transpose=TRUE)
         BetaLambda = matrix(backsolve(RiU, m1 + sdMult*rnorm(ns*(nc+nfSum))), nc+nfSum, ns, byrow=TRUE)
      } else{
         rhoVec = rhopw[rhoInd[covRhoGroup],1]
         iV_e = as.matrix(bdiag(iV, diag(rep(1,nfSum))))
         V_e = as.matrix(bdiag(chol2inv(chol(iV)), diag(rep(1,nfSum))))
         rhoVec_e = c(rhoVec, rep(0,nfSum))
         rho2Mat_e = rbind(matrix(1-rhoVec, nc, ns), priorLambda^-1)
         if(is.matrix(X)){
            XTiDX = array(XEtaTXEta, c(nc+nfSum,nc+nfSum,ns))
         } else if(is.list(X)){
            XTiDX = abind(XEtaTXEta, rev.along=0)
         }
         XTiDX = XTiDX * array(rep(iSigma, each=(nc+nfSum)^2), c(nc+nfSum,nc+nfSum,ns))
         BetaLambda = recFunSample(phyloTreeList, phyloTreeRoot, V_e, iV_e, rhoVec_e, rho2Mat_e, XTiDX, isXTS, sdMult)
      }
   }
   Beta = Mu + BetaLambda[1:nc,,drop=FALSE]
   nfCumSum = c(0,cumsum(nf*ncr)) + nc
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         Lambda[[r]] = BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE]
      } else
         Lambda[[r]] = aperm(array(BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE], c(nf[r],ncr[r],ns)),c(1,3,2))
   }
   return(list(Beta=Beta, Lambda=Lambda))
}

# if(phyloFlag==FALSE){
#    # no phylogeny information
#    diagiV = diag(iV)
#    P0 = matrix(0,nc+nfSum,nc+nfSum)
#    P0[1:nc,1:nc] = iV
#    BetaLambda = matrix(NA, nc+nfSum, ns)
#    #TODO: there is a special case for no latent factors and full data as U can be computed once
#    for(j in 1:ns){ # test whether worthy to rewrite with tensorA?
#       P = P0
#       diag(P) = c(diagiV, priorLambda[,j])
#       if(is.matrix(X)){
#          iU = P + XEtaTXEta*iSigma[j]
#       } else if(is.list(X)){
#          iU = P + XEtaTXEta[[j]]*iSigma[j]
#       }
#       RiU = chol(iU)
#       m = backsolve(RiU, backsolve(RiU, isXTS[,j], transpose=TRUE))
#       BetaLambda[,j] = M[,j] + m + backsolve(RiU, sdMult*rnorm(nc+nfSum))
#    }
# }
