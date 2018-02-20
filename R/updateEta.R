updateEta = function(Y,Z,Beta,iSigma,Eta,Lambda, X,Pi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))

   LFix = X%*%Beta
   LRan = vector("list", nr)
   for(r in 1:nr){
      LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
   }

   Eta = vector("list", nr)
   for(r in 1:nr){
      S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      lambda = Lambda[[r]]
      nf = nrow(lambda)
      lPi = Pi[,r]

      eta = matrix(NA,np[r],nf)
      if(np[r] == ny){
         iV = diag(nf) + tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
         RiV = chol(iV)
         V = chol2inv(RiV)
         mu = tcrossprod(S,lambda*matrix(iSigma,nf,ns,byrow=TRUE)) %*% V
         eta[lPi,] = mu + t(backsolve(RiV,matrix(rnorm(ny*nf),nf,ny)))
      } else{
         eta = 0*matrix(rnorm(np[r]*nf),np[r],nf)
      }
      Eta[[r]] = eta
   }
   return(Eta)
}
