updateEta = function(Y,Z,Beta,iSigma,Eta,Lambda,Alpha, rLPar, X,Pi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))

   LFix = X%*%Beta
   LRan = vector("list", nr)
   for(r in seq_len(nr)){
      LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
   }

   Eta = vector("list", nr)
   for(r in seq_len(nr)){
      if(nr > 1){
         S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      } else{
         S = Z - LFix
      }
      lambda = Lambda[[r]]
      nf = nrow(lambda)
      lPi = Pi[,r]

      if(rL[[r]]$sDim == 0){
         LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
         eta = matrix(NA,np[r],nf)
         if(np[r] == ny){
            Yx = !is.na(Y)
            indRowFull = apply(Yx,1,all)
            indRowNA = !indRowFull
            nyFull = sum(indRowFull)

            iV = diag(nf) + LamInvSigLam
            RiV = chol(iV)
            V = chol2inv(RiV)
            mu = tcrossprod(S[indRowFull,],lambda*matrix(iSigma,nf,ns,byrow=TRUE)) %*% V
            eta[lPi[indRowFull],] = mu + t(backsolve(RiV,matrix(rnorm(nyFull*nf),nf,nyFull)))

            for(i in which(indRowNA)){
               indSp = Yx[i,]
               lam = lambda[,indSp,drop=FALSE]
               iSig = iSigma[indSp]
               nsx = sum(indSp)
               LiSL = tcrossprod(lam*matrix(sqrt(iSig),nf,nsx,byrow=TRUE))
               iV = diag(nf) + LiSL
               RiV = chol(iV)
               V = chol2inv(RiV)
               mu = tcrossprod(S[i,indSp,drop=FALSE],lam*matrix(iSig,nf,nsx,byrow=TRUE)) %*% V
               eta[lPi[i],] = mu + t(backsolve(RiV,rnorm(nf)))
            }
         } else{
            # eta = 0*matrix(rnorm(np[r]*nf),np[r],nf)
            unLPi = unique(lPi)
            for(q in 1:np[r]){
               rows = (Pi[,r]==unLPi[q])
               iV = diag(nf) + LamInvSigLam*sum(rows)
               RiV = chol(iV)
               V = chol2inv(RiV)
               mu = tcrossprod(apply(S[rows,,drop=FALSE],2,sum),lambda*matrix(iSigma,nf,ns,byrow=TRUE)) %*% V
               eta[unLPi[q],] = mu + t(backsolve(RiV,rnorm(ny)))
            }
         }
      } else{
         eta = matrix(0,np[r],nf)
         iWg = rLPar[[r]]$iWg
         alpha = Alpha[[r]]
         iWs = bdiag(lapply(seq_len(nf), function(x) iWg[,,alpha[x]]))
         LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
         if(np[r] == ny){
            tmp1 = kronecker(LamInvSigLam, Diagonal(ny))
            Rtmp1 = chol(tmp1)
            fS = tcrossprod(S[order(lPi),,drop=FALSE],lambda*matrix(iSigma,nf,ns,byrow=TRUE))
            iUEta = iWs + tmp1
            R = chol(iUEta)
            L = t(R)
            tmp2 = forwardsolve(L, as.vector(fS)) + rnorm(np[r]*nf)
            feta = backsolve(R, tmp2);
            eta = matrix(feta,np[r],nf);
         } else{
            Pr = diag(np[r])[lPi,]
            tmp1 = kronecker(LamInvSigLam, Diagonal(x=apply(Pr,2,sum)))
            fS = tcrossprod(crossprod(Pr,S), lambda*matrix(iSigma,nf,ns,byrow=TRUE))
            iUEta = iWs + tmp1
            R = chol(iUEta)
            L = t(R)
            tmp2 = forwardsolve(L, as.vector(fS)) + rnorm(np[r]*nf)
            feta = backsolve(R, tmp2);
            eta = matrix(feta,np[r],nf);
         }
      }
      Eta[[r]] = eta
   }
   return(Eta)
}

