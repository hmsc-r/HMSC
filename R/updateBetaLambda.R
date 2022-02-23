# @title updateBetaLambda
#
# @description updates beta lambda
#
#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix
#' @importFrom tensorflow tf
#'
updateBetaLambda = function(Z,Gamma,iV,iD,Eta,Psi,Delta, iQ, X,Tr,Pi,dfPi,C,rL, tfCompFlag=FALSE){
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
      switch(class(X)[1L],
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
   iDS = iD*S
   iDS[is.na(Z)] = 0
   switch(class(X)[1L],
      matrix = {
         XEtaTXEta = crossprod(XEta)
         XTiDS = crossprod(XEta,iDS)
      },
      list = {
         XEtaTXEta = lapply(XEta, crossprod)
         XTiDS = matrix(NA,nc+nfSum,ns)
         for(j in 1:ns)
            XTiDS[,j] = crossprod(XEta[[j]],iDS[,j])
      }
   )

   if(is.null(C)){ # no phylogeny information
      BetaLambda = matrix(NA, nc+nfSum, ns)
      randEps = matrix(rnorm((nc+nfSum)*ns), nc+nfSum, ns)
      if(tfCompFlag==TRUE){
         tfla = tf$linalg
         ic = function(...) as.integer(c(...))
         XEtaTiDXEta = tfla$einsum("ia,ij,ib->jab",XEta,iD,XEta)
         iV_op = tfla$LinearOperatorFullMatrix(iV)
         DiagPriorLambda_op = tfla$LinearOperatorDiag(t(priorLambda))
         P = tfla$LinearOperatorBlockDiag(list(iV_op,DiagPriorLambda_op))$to_dense()
         iU = P + XEtaTiDXEta
         LiU = tfla$cholesky(iU)
         m = tfla$cholesky_solve(LiU, tf$matmul(P,tf$expand_dims(t(Mu),ic(-1))) + tf$expand_dims(t(XTiDS),ic(-1)))
         res = m + tfla$triangular_solve(LiU, tf$expand_dims(t(randEps),ic(-1)), adjoint=TRUE)
         BetaLambda = t(tf$squeeze(res,ic(-1))$numpy())
      } else{
         diagiV = diag(iV)
         P0 = matrix(0,nc+nfSum,nc+nfSum)
         P0[1:nc,1:nc] = iV
         for(j in 1:ns){
            if(class(X)[1L]=="matrix"){
               XEtaTiDXEta = crossprod(XEta,iD[,j]*XEta)
            } else if(class(X)[1L]=="list"){
               XEtaTiDXEta = crossprod(XEta[[j]],iD[,j]*XEta[[j]])
            }
            P = P0
            diag(P) = c(diagiV, priorLambda[,j])
            iU = P + XEtaTiDXEta
            RiU = chol(iU)
            U = chol2inv(RiU)
            m = U %*% (P%*%Mu[,j] + XTiDS[,j]);
            BetaLambda[,j] = m + backsolve(RiU, randEps[,j])
         }
      }

   } else{ # available phylogeny information
      P = bdiag(kronecker(iV,iQ), Diagonal(x=t(priorLambda)))
      tmp = vector("list",ns)
      for(j in 1:ns){
         if(class(X)[1L]=="matrix"){
            XEtaTiDXEta = crossprod(XEta,iD[,j]*XEta)
         } else if(class(X)[1L]=="list"){
            XEtaTiDXEta = crossprod(XEta[[j]],iD[,j]*XEta[[j]])
         }
         tmp[[j]] = XEtaTiDXEta
      }
      tmpMat = Reduce(rbind, tmp)
      ind1 = rep(rep(1:ns,each=nc+nfSum)+ns*rep(0:(nc+nfSum-1),ns), nc+nfSum)
      ind2 = rep(1:((nc+nfSum)*ns), each=nc+nfSum)
      mat = sparseMatrix(ind1, ind2, x=as.vector(tmpMat))
      RiU = chol(as.matrix(mat) + P)
      m1 = backsolve(RiU, P%*%as.vector(t(Mu)) + as.vector(t(XTiDS)), transpose=TRUE)
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

