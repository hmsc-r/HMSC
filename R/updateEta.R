#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix
#' 
updateEta = function(Y,Z,Beta,iSigma,Eta,Lambda,Alpha, rLPar, X,Pi,dfPi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   Yx = !is.na(Y)

   switch(class(X),
          matrix = {
             LFix = X%*%Beta
          },
          list = {
             LFix = matrix(NA,ny,ns)
             for(j in 1:ns)
                LFix[,j] = X[[j]]%*%Beta[,j]
          }
   )
   LRan = vector("list", nr)
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
      } else{
         LRan[[r]] = matrix(0,ny,ns)
         for(k in 1:rL[[r]]$xDim)
            LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),r]) %*% Lambda[[r]][,,r]
      }
   }
   for(r in seq_len(nr)){
      rnames=rownames(Eta[[r]])
      if(nr > 1){
         S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      } else{
         S = Z - LFix
      }
      lambda = Lambda[[r]]
      nf = dim(lambda)[1]
      lPi = Pi[,r]
      ldfPi = dfPi[,r]
      if(rL[[r]]$sDim == 0){
         eta = matrix(NA,np[r],nf)
         if(rL[[r]]$xDim == 0){
            LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
            if(np[r] == ny){
               indRowFull = apply(Yx,1,all)
               indRowNA = !indRowFull
               nyFull = sum(indRowFull)

               if(nyFull > 0){
                  iV = diag(nf) + LamInvSigLam
                  RiV = chol(iV)
                  V = chol2inv(RiV)
                  mu = tcrossprod(S[indRowFull,],lambda*matrix(iSigma,nf,ns,byrow=TRUE)) %*% V
                  eta[lPi[indRowFull],] = mu + t(backsolve(RiV,matrix(rnorm(nyFull*nf),nf,nyFull)))
               }

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
               unLPi = unique(lPi)
               for(q in 1:np[r]){
                  rows = which(lPi==unLPi[q])
                  if(all(Yx[rows,])){
                     iV = diag(nf) + LamInvSigLam*length(rows)
                     RiV = chol(iV)
                     V = chol2inv(RiV)
                     mu = tcrossprod(apply(S[rows,,drop=FALSE],2,sum), lambda*matrix(iSigma,nf,ns,byrow=TRUE)) %*% V
                  } else{
                     LiSL = matrix(0,nf,nf)
                     for(p in 1:length(rows))
                        LiSL = LiSL + tcrossprod(lambda*matrix(sqrt(iSigma)*Yx[rows[p],],nf,ns,byrow=TRUE))
                     iV = diag(nf) + LiSL
                     RiV = chol(iV)
                     V = chol2inv(RiV)
                     mu = colSums(tcrossprod(S[rows,,drop=FALSE]*matrix(iSigma,length(rows),ns,byrow=TRUE)*Yx[rows,], lambda)) %*% V
                  }

                  eta[unLPi[q],] = mu + t(backsolve(RiV,rnorm(ny)))
               }
            }
         } else{
            ncr = rL[[r]]$xDim
            unLPi = unique(lPi)
            unLdfPi = unique(as.character(ldfPi))
            for(q in 1:np[r]){
               lambdaLocal = rowSums(lambda * array(rep(rL[[r]]$x[unLdfPi[q],],each=nf*ns), c(nf,ns,ncr)), dims=2)
               rows = which(lPi==unLPi[q])
               LiSL = matrix(0,nf,nf)
               for(p in 1:length(rows))
                  LiSL = LiSL + tcrossprod(lambdaLocal * matrix(sqrt(iSigma)*Yx[rows[p],],nf,ns,byrow=TRUE))
               iV = diag(nf) + LiSL
               RiV = chol(iV)
               V = chol2inv(RiV)
               mu = colSums(tcrossprod(S[rows,,drop=FALSE]*matrix(iSigma,length(rows),ns,byrow=TRUE)*Yx[rows,], lambdaLocal)) %*% V
               eta[unLPi[q],] = mu + t(backsolve(RiV,rnorm(ny)))
            }
         }
      } else{
         eta = matrix(0,np[r],nf)
         alpha = Alpha[[r]]
         iWg = rLPar[[r]]$iWg
         switch(rL[[r]]$spatialMethod,
                "Full" = {
                   iWs = bdiag(lapply(seq_len(nf), function(x) iWg[,,alpha[x]]))
                   LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
                   if(np[r] == ny){
                      tmp1 = kronecker(LamInvSigLam, Diagonal(ny))
                      Rtmp1 = chol(tmp1)
                      fS = tcrossprod(S[order(lPi),,drop=FALSE],lambda*matrix(iSigma,nf,ns,byrow=TRUE))
                      iUEta = iWs + tmp1
                      R = chol(iUEta)
                      tmp2 = backsolve(R, as.vector(fS), transpose=TRUE) + rnorm(np[r]*nf)
                      feta = backsolve(R, tmp2);
                      eta = matrix(feta,np[r],nf);
                   } else{
                      P = sparseMatrix(i=1:ny,j=lPi)
                      tmp1 = kronecker(LamInvSigLam, Diagonal(x=Matrix::colSums(P)))
                      fS = Matrix::tcrossprod(Matrix::crossprod(P,S), lambda*matrix(iSigma,nf,ns,byrow=TRUE))
                      iUEta = iWs + tmp1
                      R = chol(iUEta)
                      tmp2 = backsolve(R, as.vector(fS), transpose=TRUE) + rnorm(np[r]*nf)
                      feta = backsolve(R, tmp2);
                      eta = matrix(feta,np[r],nf);
                   }
                },
                "NNGP" = {
                   iWs = bdiag(lapply(seq_len(nf), function(x) iWg[[alpha[x]]]))
                   LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
                   P = sparseMatrix(i=1:ny,j=lPi)
                   tmp1 = kronecker(LamInvSigLam, Diagonal(x=Matrix::colSums(P)))
                   fS = Matrix::tcrossprod(Matrix::crossprod(P,S), lambda*matrix(iSigma,nf,ns,byrow=TRUE))
                   iUEta = iWs + tmp1
                   R = chol(iUEta)
                   tmp2 = backsolve(R, as.vector(fS), transpose=TRUE) + rnorm(np[r]*nf)
                   feta = backsolve(R, tmp2);
                   eta = matrix(feta,np[r],nf);
                }
                )
      }
      rownames(eta)=rnames
      Eta[[r]] = eta
      if(r < nr){
         if(rL[[r]]$xDim == 0){
            LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
         } else{
            LRan[[r]] = matrix(0,ny,ns)
            for(k in 1:rL[[r]]$xDim)
               LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),r]) %*% Lambda[[r]][,,r]
         }
      }
   }
   return(Eta)
}

