updateBetaLambda = function(Z,Gamma,iV,iSigma,Eta,Psi,Delta,rho, iQg, X,Tr,Pi,C){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)

   LFix = 0
   S = Z - LFix

   if(!is.null(nr)){
      EtaFull = vector("list", nr)
      Lambda = vector("list", nr)
      nf = rep(NA,nr)
      for(r in seq_len(nr)){
         EtaFull[[r]] = Eta[[r]][Pi[,r],]
         nf[r] = ncol(Eta[[r]])
      }
      nfSum = sum(nf)
   } else
      nfSum = 0

   EtaSt = matrix(unlist(EtaFull),ny,nfSum)
   XEta = cbind(X,EtaSt)
   Q = crossprod(XEta)

   Mu = rbind(tcrossprod(Gamma,Tr), matrix(0,nfSum,ns))
   isXTS = crossprod(XEta,S) * matrix(iSigma,nc+nfSum,ns,byrow=TRUE)
   psiSt = matrix(unlist(lapply(Psi, t)),nfSum,ns,byrow=TRUE)
   Tau = lapply(Delta, function(a) apply(a,2,cumprod))
   tauSt = matrix(unlist(Tau),nfSum,1)
   priorLambda = psiSt*matrix(tauSt,nfSum,ns)

   if(is.null(C)){
      # no phylogeny information
      diagiV = diag(iV)
      P0 = matrix(0,nc+nfSum,nc+nfSum)
      P0[1:nc,1:nc] = iV
      BetaLambda = matrix(NA, nc+nfSum, ns)
      for(j in 1:ns){
         P = P0
         diag(P) = c(diagiV, priorLambda[,j])
         iU = P + Q*iSigma[j]
         RiU = chol(iU)
         U = chol2inv(RiU)
         m = U %*% (P%*%Mu[,j] + isXTS[,j]);
         BetaLambda[,j] = m + backsolve(RiU, rnorm(nc+nfSum))
      }
   } else{
      # available phylogeny information
      P = bdiag(kronecker(iV,iQg[,,rho]), Diagonal(x=t(priorLambda)))
      RiU = chol(kronecker(Q,diag(iSigma)) + P)
      U = chol2inv(RiU)
      m = U %*% (P%*%as.vector(t(Mu)) + as.vector(t(isXTS)))
      BetaLambda = matrix(m + backsolve(RiU, rnorm(ns*(nc+nfSum))), nc+nfSum, ns, byrow=TRUE)
   }
   Beta = BetaLambda[1:nc,]
   for(r in seq_len(nr)){
      Lambda[[r]] = BetaLambda[nc+sum(nf[seq_len(r-1)])+(1:nf[r]),]
   }
   return(list(Beta=Beta, Lambda=Lambda))
}

