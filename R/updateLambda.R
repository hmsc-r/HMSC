updateLambda = function(Z,Beta,iSigma,Eta,Lambda,Psi,Delta, X,Pi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)

   LFix = X%*%Beta
   LRan = vector("list", nr)
   for(r in 1:nr){
      LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
   }

   Lambda = vector("list", nr)
   for(r in 1:nr){
      S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      delta = Delta[[r]]
      psi = Psi[[r]]
      eta = Eta[[r]][Pi[,r],]
      nf = ncol(eta)
      tau = apply(delta, 2, cumprod)
      iSigmaEtaS = crossprod(eta, S) * matrix(iSigma,nf,ns,byrow=TRUE)

      priorLambda = psi*matrix(tau,nf,ns)
      E = crossprod(eta)
      lambda = matrix(NA,nf,ns)
      for(j in 1:ns){
         Q = diag(priorLambda[,j]) + E*iSigma[j]
         R = chol(Q);
         z = rnorm(nf)
         lambda[,j] = backsolve(R, backsolve(R,iSigmaEtaS[,j],transpose=T)+z)
      }
      j = 1
      Lambda[[r]] = lambda
   }
   return(Lambda)
}
