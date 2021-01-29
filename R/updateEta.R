#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#'
updateEta = function(Y,Z,Beta,iSigma,Eta,Lambda,Alpha, rLPar, X,Pi,dfPi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   Yx = !is.na(Y)

   switch(class(X)[1L],
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
               lambdaLocal = rowSums(lambda * array(unlist(rep(rL[[r]]$x[unLdfPi[q],],each=nf*ns)), c(nf,ns,ncr)), dims=2)
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
                   iWs = sparseMatrix(c(),c(),dims=c(np[r]*nf,np[r]*nf))
                   for(h in seq_len(nf))
                      iWs = iWs + kronecker(iWg[[alpha[h]]], Diagonal(x=c(rep(0,h-1),1,rep(0,nf-h))))
                   LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
                   # TODO it is possible here to eliminate dependence on missing data at minor cost by using row-specific iSigma
                   if(np[r] == ny){
                      tmp1 = kronecker(Diagonal(ny), LamInvSigLam)
                      fS = tcrossprod(S[order(lPi),,drop=FALSE], lambda*matrix(iSigma,nf,ns,byrow=TRUE))
                   }else{
                      P = sparseMatrix(i=1:ny,j=lPi)
                      tmp1 = kronecker(Diagonal(x=Matrix::colSums(P)), LamInvSigLam)
                      fS = Matrix::tcrossprod(Matrix::crossprod(P,S), lambda*matrix(iSigma,nf,ns,byrow=TRUE))
                   }
                   iUEta = iWs + tmp1
                   R = chol(iUEta)
                   tmp2 = backsolve(R, as.vector(t(fS)), transpose=TRUE) + rnorm(nf*np[r])
                   feta = backsolve(R, tmp2)
                   eta = matrix(feta,np[r],nf,byrow=TRUE)
                },
                "GPP" = {
                   if(np[r] == ny){
                      idDg = rLPar[[r]]$idDg
                      idDW12g = rLPar[[r]]$idDW12g
                      Fg = rLPar[[r]]$Fg
                      nK = nrow(Fg)
                      idD = idDg[,alpha]
                      Fmat = matrix(0,nrow=(nK*nf),ncol=(nK*nf))
                      idD1W12 = matrix(0,nrow=(np[r]*nf),ncol=(nK*nf))
                      for(h in 1:nf){
                         Fmat[(h-1)*nK+(1:nK), (h-1)*nK+(1:nK)] = Fg[,,alpha[h]]
                         idD1W12[(h-1)*np[r]+(1:np[r]), (h-1)*nK+(1:nK)] = idDW12g[,,alpha[h]]
                      }
                      tmp = diag(iSigma,length(iSigma))%*%t(lambda)
                      fS = S[order(lPi),,drop=FALSE]%*%tmp
                      fS = matrix(fS,ncol=1)
                      LamSigLamT = lambda%*%tmp

                      B0 = array(LamSigLamT,c(nrow(LamSigLamT),ncol(LamSigLamT),ny))
                      idDV = matrix(t(idD),nrow=1)
                      tmp = t(matrix((c(1:ny)-1),nrow=ny,ncol=nf))
                      ind = matrix(tmp,nrow=1)*nf^2 + rep(nf*((1:nf)-1)+(1:nf),ny)
                      B0[ind] = B0[ind] + idDV

                      B1 = array(NA,c(ny,nf,nf))
                      LB1 = array(NA,c(ny,nf,nf))
                      for(i in 1:ny){
                         B1[i,,] = chol2inv(chol((B0[,,i])))
                         LB1[i,,] = t(chol(B1[i,,]))
                      }
                      ind1 = rep(1:(nf*ny),nf)
                      tmp1 = t(matrix((1:nf-1),nrow=nf,ncol=(ny*nf))) * ny
                      ind2 = rep(1:ny,nf^2) + t(matrix(tmp1,nrow=1))
                      iA = Matrix(0,nrow=nf*ny, ncol=nf*ny,sparse=TRUE)
                      iA[t(matrix(rbind(ind1,as.vector(ind2)),nrow=2))] = as.vector(B1)
                      LiA = Matrix(0,nrow=nf*ny, ncol=nf*ny,sparse=TRUE)
                      LiA[t(matrix(rbind(ind1,as.vector(ind2)),nrow=2))] = as.vector(LB1)
                      iAidD1W12 = iA %*% idD1W12
                      H = Fmat - t(idD1W12)%*%iAidD1W12
                      RH = chol(as.matrix(H))
                      iRH = solve(RH)

                      mu1 = iA%*%fS
                      tmp1 = iAidD1W12 %*% iRH
                      mu2 = tmp1%*%(Matrix::t(tmp1)%*%fS)

                      etaR = LiA%*%rnorm(np[r]*nf)+tmp1%*%rnorm(nK*nf)
                      eta = matrix(mu1+mu2+etaR,ncol=nf,nrow=np[r])
                   } else {
                      idDg = rLPar[[r]]$idDg
                      idDW12g = rLPar[[r]]$idDW12g
                      Fg = rLPar[[r]]$Fg
                      nK = nrow(Fg)
                      idD = idDg[,alpha]
                      Fmat = matrix(0,nrow=(nK*nf),ncol=(nK*nf))
                      idD1W12 = matrix(0,nrow=(np[r]*nf),ncol=(nK*nf))
                      for(h in 1:nf){
                         Fmat[(h-1)*nK+(1:nK), (h-1)*nK+(1:nK)] = Fg[,,alpha[h]]
                         idD1W12[(h-1)*np[r]+(1:np[r]), (h-1)*nK+(1:nK)] = idDW12g[,,alpha[h]]
                      }
                      tmp = diag(iSigma,length(iSigma))%*%t(lambda)
                      LamSigLamT = lambda%*%tmp

                      P = sparseMatrix(i=1:ny,j=lPi)
                      f = Matrix::crossprod(P,S)
                      fS = f%*%tmp
                      fS = matrix(fS,ncol=1)

                      tmp1 = kronecker(LamSigLamT, Diagonal(x=Matrix::colSums(P)))
                      tmp2 = tmp1 + Diagonal(x=idD[])
                      iA = solve(tmp2)
                      LiA = chol(iA)

                      iAidD1W12 = iA %*% idD1W12
                      H = Fmat - t(idD1W12)%*%iAidD1W12
                      RH = chol(as.matrix(H))
                      iRH = solve(RH)

                      mu1 = iA%*%fS
                      tmp1 = iAidD1W12 %*% iRH
                      mu2 = tmp1%*%(Matrix::t(tmp1)%*%fS)

                      etaR = LiA%*%rnorm(np[r]*nf)+tmp1%*%rnorm(nK*nf)
                      eta = matrix(mu1+mu2+etaR,ncol=nf,nrow=np[r])
                   }
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

