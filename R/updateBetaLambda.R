#' @title updateBetaLambda
#'
#' @description updates beta lambda
#'
updateBetaLambda = function(Y,Z,Gamma,iV,iSigma,Eta,Psi,Delta,rho, iQ, X,Tr,Pi,C){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)

   LFix = 0
   S = Z - LFix

   Lambda = vector("list", nr)
   if(nr > 0){
      EtaFull = vector("list", nr)
      nf = rep(NA,nr)
      for(r in seq_len(nr)){
         EtaFull[[r]] = Eta[[r]][Pi[,r],]
         nf[r] = ncol(Eta[[r]])
      }
      nfSum = sum(nf)
      EtaSt = matrix(unlist(EtaFull),ny,nfSum)
      switch(class(X),
         matrix = {
            XEta = cbind(X,EtaSt)
         },
         list = {
            XEta = lapply(X, cbind, EtaSt)
         }
      )

      psiSt = matrix(unlist(lapply(Psi, t)),nfSum,ns,byrow=TRUE)
      Tau = lapply(Delta, function(a) apply(a,2,cumprod))
      tauSt = matrix(unlist(Tau),nfSum,1)
      priorLambda = psiSt*matrix(tauSt,nfSum,ns)
   } else{
      nfSum = 0
      XEta = X
      priorLambda = matrix(numeric(0),0,ns)
   }

   Mu = rbind(tcrossprod(Gamma,Tr), matrix(0,nfSum,ns))
   switch(class(X),
      matrix = {
         XEtatXEta = crossprod(XEta)
         isXTS = crossprod(XEta,S) * matrix(iSigma,nc+nfSum,ns,byrow=TRUE)
      },
      list = {
         XEtatXEta = lapply(XEta, crossprod)
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
               iU = P + XEtatXEta*iSigma[j]
            },
            list = {
               iU = P + XEtatXEta[[j]]*iSigma[j]
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
               XEtatXEta = crossprod(XEta[indObs,,drop=FALSE])
               isXTS = crossprod(XEta[indObs,,drop=FALSE],S[indObs,j]) * iSigma[j]
            },
            list = {
               XEtatXEta = crossprod(XEta[[j]][indObs,,drop=FALSE])
               isXTS = crossprod(XEta[[j]][indObs,,drop=FALSE],S[indObs,j]) * iSigma[j]
            }
         )
         P = P0
         diag(P) = c(diagiV, priorLambda[,j])
         iU = P + XEtatXEta*iSigma[j]
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
            RiU = chol(kronecker(XEtatXEta,diag(iSigma)) + P)
         },
         list = {
            tmp = vector("list",ns)
            for(j in 1:ns)
               tmp[[j]] = XEtatXEta[[j]] * iSigma[j]
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
   Beta = BetaLambda[1:nc,,drop=F]
   for(r in seq_len(nr)){
      Lambda[[r]] = BetaLambda[nc+sum(nf[seq_len(r-1)])+(1:nf[r]),,drop=FALSE]
   }
   return(list(Beta=Beta, Lambda=Lambda))
}

