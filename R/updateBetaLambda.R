# @title updateBetaLambda
#
# @description updates beta lambda
#
#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix
#'
updateBetaLambda = function(Y,Z,Gamma,iV,iSigma,Eta,Psi,Delta,rho, iQ, X,Tr,Pi,dfPi,C,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)

   S = Z
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
      switch(class(X),
         matrix = {
            XEta = cbind(X,EtaSt)
         },
         list = {
            XEta = lapply(X, cbind, EtaSt)
         }
      )
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

   Mu = rbind(tcrossprod(Gamma,Tr), matrix(0,nfSum,ns))
   switch(class(X),
      matrix = {
         XEtaTXEta = crossprod(XEta)
         isXTS = crossprod(XEta,S) * matrix(iSigma,nc+nfSum,ns,byrow=TRUE)
      },
      list = {
         XEtaTXEta = lapply(XEta, crossprod)
         isXTS = matrix(NA,nc+nfSum,ns)
         for(j in 1:ns)
            isXTS[,j] = crossprod(XEta[[j]],S[,j]) * iSigma[j]
      }
   )

   if(is.null(C)){
      # no phylogeny information
      Yx = !is.na(Y)
      indColFull = apply(Yx,2,all)
      indColNA = !indColFull

      diagiV = diag(iV)
      P0 = matrix(0,nc+nfSum,nc+nfSum)
      P0[1:nc,1:nc] = iV
      BetaLambda = matrix(NA, nc+nfSum, ns)

      for(j in which(indColFull)){ # test whether worthy to rewrite with tensorA?
         P = P0
         diag(P) = c(diagiV, priorLambda[,j])
         switch(class(X),
            matrix = {
               iU = P + XEtaTXEta*iSigma[j]
            },
            list = {
               iU = P + XEtaTXEta[[j]]*iSigma[j]
            }
         )
         RiU = chol(iU)
         U = chol2inv(RiU)
         m = U %*% (P%*%Mu[,j] + isXTS[,j]);
         BetaLambda[,j] = m + backsolve(RiU, rnorm(nc+nfSum))
      }
      for(j in which(indColNA)){
         indObs = Yx[,j]
         switch(class(X),
            matrix = {
               XEtaTXEta = crossprod(XEta[indObs,,drop=FALSE])
               isXTS = crossprod(XEta[indObs,,drop=FALSE],S[indObs,j]) * iSigma[j]
            },
            list = {
               XEtaTXEta = crossprod(XEta[[j]][indObs,,drop=FALSE])
               isXTS = crossprod(XEta[[j]][indObs,,drop=FALSE],S[indObs,j]) * iSigma[j]
            }
         )
         P = P0
         diag(P) = c(diagiV, priorLambda[,j])
         iU = P + XEtaTXEta*iSigma[j]
         RiU = chol(iU)
         U = chol2inv(RiU)
         m = U %*% (P%*%Mu[,j] + isXTS);
         BetaLambda[,j] = m + backsolve(RiU, rnorm(nc+nfSum))
      }

   } else{
      # available phylogeny information
      P = bdiag(kronecker(iV,iQ), Diagonal(x=t(priorLambda)))
      switch(class(X),
         matrix = {
            RiU = chol(kronecker(XEtaTXEta,diag(iSigma)) + P)
         },
         list = {
            tmp = vector("list",ns)
            for(j in 1:ns)
               tmp[[j]] = XEtaTXEta[[j]] * iSigma[j]
            tmpMat = Reduce(rbind, tmp)
            ind1 = rep(rep(1:ns,each=nc+nfSum)+ns*rep(0:(nc+nfSum-1),ns), nc+nfSum)
            ind2 = rep(1:((nc+nfSum)*ns), each=nc+nfSum)
            mat = sparseMatrix(ind1, ind2, x=as.vector(tmpMat))
            RiU = chol(as.matrix(mat) + P)
         }
      )
      # U = chol2inv(RiU)
      # m = U %*% (P%*%as.vector(t(Mu)) + as.vector(t(isXTS)))
      # BetaLambda = matrix(m + backsolve(RiU, rnorm(ns*(nc+nfSum))), nc+nfSum, ns, byrow=TRUE)
      m1 = backsolve(RiU, P%*%as.vector(t(Mu)) + as.vector(t(isXTS)), transpose=TRUE)
      BetaLambda = matrix(backsolve(RiU, m1 + rnorm(ns*(nc+nfSum))), nc+nfSum, ns, byrow=TRUE)
   }
   Beta = BetaLambda[1:nc,,drop=FALSE]
   nfCumSum = c(0,cumsum(nf*ncr)) + nc
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         Lambda[[r]] = BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE]
      } else
         Lambda[[r]] = aperm(array(BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE], c(nf[r],ncr[r],ns)),c(1,3,2))
   }
   return(list(Beta=Beta, Lambda=Lambda))
}

