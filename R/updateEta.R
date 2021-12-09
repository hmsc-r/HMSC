#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix chol
#' @importFrom SparseM chol backsolve
#' @importFrom plyr mdply
#'
updateEta = function(Z,Beta,iSigma,Eta,Lambda,Alpha, rLPar, Y,X,Pi,dfPi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   Yx = !is.na(Y)
   ic = function(...){
      return(as.integer(c(...)))
   }
   tfla = tf$linalg
   tfm = tf$math
   EtaFullList = vector("list",nr)

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
      LRan[[r]] = computePredictor.HmscRandomLevel(rL[[r]], Eta[[r]], Lambda[[r]], Pi[,r], dfPi[,r])
   }
   for(r in seq_len(nr)){
      EtaFull = NULL
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
      if(inherits(rL[[r]],"HmscRandomLevel",TRUE)==1){
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
               for(i in which(indRowNA)){ #TODO write a batched version to avoid iterating in the loop
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
               if(!exists(sprintf("updateEta_%.3d",r), envir=parent.frame())){
                  fun = function(Lambda,iSigma,S,np,nf){
                     iSigmaObs = tf$einsum("j,ij->ij",iSigma,Yx)
                     iVSt = tf$scatter_nd(matrix(ic(lPi-1)),tf$einsum("hj,ij,fj->ihf",Lambda,iSigmaObs,Lambda),tf$stack(c(np,nf,nf))) + tf$eye(nf,dtype=tf$float64)
                     LiVSt = tfla$cholesky(iVSt)
                     lambdaiSigmaS = tf$scatter_nd(matrix(ic(lPi-1)),tf$einsum("hj,ij,ij->ih",Lambda,iSigmaObs,S),tf$stack(c(np,nf)))
                     iLiVStlambdaiSigmaS = tfla$triangular_solve(LiVSt, lambdaiSigmaS[,,NULL])
                     eta = tf$squeeze(tfla$triangular_solve(LiVSt, iLiVStlambdaiSigmaS+tf$random$normal(shape=tf$stack(c(np,nf,ic(1))),dtype=tf$float64), adjoint=TRUE),ic(-1))
                     return(eta)
                  }
                  fun_tf = tf_function(fun)
                  assign(sprintf("updateEta_%.3d",r), fun_tf, envir=parent.frame())
               } else{
                  fun_tf = get(sprintf("updateEta_%.3d",r), envir=parent.frame())
               }
               eta_tf = fun_tf(tf$constant(lambda,dtype=tf$float64),tf$constant(iSigma,dtype=tf$float64),tf$constant(S,dtype=tf$float64),tf$constant(ic(np[r]),tf$int32),tf$constant(ic(nf),tf$int32))
               eta = eta_tf$numpy()
               # unLPi = unique(lPi)
               # for(q in 1:np[r]){
               #    rows = which(lPi==unLPi[q])
               #    if(all(Yx[rows,])){
               #       iV = diag(nf) + LamInvSigLam*length(rows)
               #       RiV = chol(iV)
               #       V = chol2inv(RiV)
               #       mu = tcrossprod(apply(S[rows,,drop=FALSE],2,sum), lambda*matrix(iSigma,nf,ns,byrow=TRUE)) %*% V
               #    } else{
               #       LiSL = matrix(0,nf,nf)
               #       for(p in 1:length(rows))
               #          LiSL = LiSL + tcrossprod(lambda*matrix(sqrt(iSigma)*Yx[rows[p],],nf,ns,byrow=TRUE))
               #       iV = diag(nf) + LiSL
               #       RiV = chol(iV)
               #       V = chol2inv(RiV)
               #       mu = colSums(tcrossprod(S[rows,,drop=FALSE]*matrix(iSigma,length(rows),ns,byrow=TRUE)*Yx[rows,], lambda)) %*% V
               #    }
               #    eta[unLPi[q],] = mu + 0*t(backsolve(RiV,rnorm(ny)))
               # }
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
      } else if(inherits(rL[[r]],"HmscSpatialRandomLevel",TRUE)==1){
         eta = matrix(0,np[r],nf)
         alpha = Alpha[[r]]
         iWg = rLPar[[r]]$iWg
         switch(rL[[r]]$spatialMethod,
                "Full" = {
                   iWs = bdiag(lapply(seq_len(nf), function(x) iWg[[alpha[x]]]))
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
      } else if(inherits(rL[[r]],"HmscKroneckerRandomLevel",TRUE)==1){
         alpha = Alpha[[r]]
         iD = !is.na(Y) * matrix(iSigma,ny,ns,byrow=TRUE)
         LamiDLam = tf$scatter_nd(matrix(ic(lPi-1)), tf$einsum("hj,ij,kj->ihk", lambda,iD,lambda), ic(np[r],nf,nf))$numpy()
         fS = tf$scatter_nd(matrix(ic(lPi-1)), tcrossprod(iD*S, lambda) , ic(np[r],nf))$numpy()

         epsRand = rnorm(np[r]*nf)
         dfPiElemLevelList = vector("list", length(rL[[r]]$rLList))
         npElemVec = rep(NA, length(rL[[r]]$rLList))
         for(l in seq_len(length(rL[[r]]$rLList))){
            dfPiElem = as.factor(unlist(lapply(strsplit(as.character(dfPi[,r]), rL[[r]]$sepStr), function(a) a[l])))
            dfPiElemLevelList[[l]] = levels(dfPiElem)
            npElemVec[l] = length(dfPiElemLevelList[[l]])
         }
         dfTmp = expand.grid(rev(dfPiElemLevelList))[,rev(1:length(rL[[r]]$rLList))]
         allUnits = factor(do.call(function(...) paste(..., sep=rL[[r]]$sepStr),dfTmp))
         indKronObs = as.numeric(factor(levels(dfPi[,r]), levels=levels(allUnits)))
         if(rL[[r]]$etaMethod == "R"){
            if(rL[[r]]$sDim == 0){
               iWs = Diagonal(np[r]*nf)
            } else{
               iWsList = vector("list", nf)
               for(h in seq_len(nf)){
                  if(length(allUnits) == np[r]){
                     iWfList = vector("list", length(rL[[r]]$rLList))
                     for(l in seq_len(length(rL[[r]]$rLList))){
                        if(rL[[r]]$rLList[[l]]$sDim == 0){
                           iWfList[[l]] = Diagonal(npElemVec[l])
                        } else{
                           iWfList[[l]] = rLPar[[r]][[l]]$iWg[[alpha[h,l]]]
                        }
                     }
                     iWFull = iWfList[[1]]
                     for(l in 1+seq_len(length(rL[[r]]$rLList)-1)){ #TODO can this be done without loop?
                        iWFull = Matrix::kronecker(iWFull, iWfList[[l]])
                     }
                     iWsList[[h]] = iWFull
                  } else{
                     WfList = vector("list", length(rL[[r]]$rLList))
                     for(l in seq_len(length(rL[[r]]$rLList))){
                        if(rL[[r]]$rLList[[l]]$sDim == 0){
                           WfList[[l]] = Diagonal(npElemVec[l])
                        } else{
                           WfList[[l]] = rLPar[[r]][[l]]$Wg[[alpha[h,l]]]
                        }
                     }
                     WFull = WfList[[1]]
                     for(l in 1+seq_len(length(rL[[r]]$rLList)-1)){ #TODO can this be done without loop?
                        WFull = Matrix::kronecker(WFull, WfList[[l]])
                     }
                     Wf = WFull[indKronObs,indKronObs]
                     iWsList[[h]] = chol2inv(chol(Wf))
                  }
               }
               iWs = sparseMatrix(c(),c(),dims=c(np[r]*nf,np[r]*nf))
               for(h in seq_len(nf))
                  iWs = iWs + Matrix::kronecker(iWsList[[h]], Diagonal(x=c(rep(0,h-1),1,rep(0,nf-h))))
            }
            tmp1 = bdiag(asplit(LamiDLam, 1))
            iUEta = iWs + tmp1
            R = Matrix::chol(iUEta)
            tmp2 = Matrix::solve(t(R), as.vector(t(fS)))
            feta = Matrix::solve(R, tmp2+epsRand)
            eta = matrix(feta,np[r],nf,byrow=TRUE)
         } else{
            nt = npElemVec[1]; nx = npElemVec[2]
            if(np[r] == nx*nt){
               M0Full = fS
               epsFullRand = epsRand
            } else{
               M0Full = matrix(0,nx*nt,nf)
               M0Full[indKronObs,] = fS
               epsFullRand = rnorm(nx*nt*nf)
            }
            if(rL[[r]]$etaMethod == "TF"){
               if(!exists(sprintf("updateEta_kronecker_TF_%.3d",r), envir=parent.frame())){
                  fun = function(alpha,LamiDLam,m0Mat,epsMat){
                     iKtArray = tf$gather(rLPar[[r]][[1]]$iWStack, alpha[,1])
                     iKsArray = tf$gather(rLPar[[r]][[2]]$iWStack, alpha[,2])
                     LamiDLam = tf$reshape(tf$scatter_nd(matrix(ic(indKronObs-1)), LamiDLam, ic(nt*nx,nf,nf)), ic(nt,nx,nf,nf))
                     it = tf$constant(ic(0),tf$int32)
                     diagBlockiKtVec = iKtArray[,2,it]
                     A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
                     # A = A + tf$reshape(tf$transpose(tfla$diag(LamiDLam[,,NULL]*tf$constant(numKronObsMat,tf$float64)[,it]), ic(2,0,3,1)), ic(nf*nx,nf*nx))
                     A = A + tf$reshape(tf$transpose(tfla$diag(tf$transpose(LamiDLam[it,,,], ic(1,2,0))), ic(2,0,3,1)), ic(nx*nf,nx*nf))
                     L = tfla$cholesky(A)
                     iLm0 = tf$TensorArray(tf$float64, size=ic(nt))
                     LArray = tf$TensorArray(tf$float64, size=ic(nt))
                     CArray = tf$TensorArray(tf$float64, size=ic(nt-1))
                     iLm0 = iLm0$write(it, tf$squeeze(tfla$triangular_solve(L, m0Mat[it,,NULL]),ic(-1)))
                     LArray = LArray$write(it, L)
                     for(it in tf$range(ic(1),ic(nt))){
                        diagBlockiKtVec = iKtArray[,2,it]
                        A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
                        A = A + tf$reshape(tf$transpose(tfla$diag(tf$transpose(LamiDLam[it,,,], ic(1,2,0))), ic(2,0,3,1)), ic(nx*nf,nx*nf))
                        offdiagBlockiKtVec = iKtArray[,3,it]
                        B = tf$reshape(tf$transpose(tfla$diag(offdiagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
                        C = tf$transpose(tfla$triangular_solve(LArray$read(it-ic(1)), B))
                        L = tfla$cholesky(A-tf$matmul(C,C,transpose_b=TRUE))
                        v = m0Mat[it,] - tf$squeeze(tf$matmul(C,iLm0$read(it-ic(1))[,NULL]),ic(-1))
                        iLm0 = iLm0$write(it, tf$squeeze(tfla$triangular_solve(L, v[,NULL]),ic(-1)))
                        # iLm0$write(it, v)
                        LArray = LArray$write(it, L)
                        CArray = CArray$write(it-ic(1), C)
                     }
                     LArray = LArray$stack(); CArray = CArray$stack(); iLm0 = iLm0$stack()
                     it = tf$constant(ic(nt-1),tf$int32)
                     vv = iLm0[it,] + epsMat[it,]
                     iLTiLm0 = tf$TensorArray(tf$float64, size=ic(nt))
                     iLTiLm0 = iLTiLm0$write(it, tf$squeeze(tfla$triangular_solve(LArray[it,,], vv[,NULL], adjoint=TRUE),ic(-1)))
                     for(it in tf$range(ic(nt-2),ic(-1),ic(-1))){
                        vv = iLm0[it,] + epsMat[it,] - tf$squeeze(tf$matmul(CArray[it,,],iLTiLm0$read(it+ic(1))[,NULL],transpose_a=TRUE),ic(-1))
                        iLTiLm0 = iLTiLm0$write(it, tf$squeeze(tfla$triangular_solve(LArray[it,,], vv[,NULL], adjoint=TRUE),ic(-1)))
                        # tf$print(tmp)
                        # tf$print(iLTiLm0$stack())
                     }
                     iLTiLm0 = iLTiLm0$stack()
                     return(tf$reshape(iLTiLm0, ic(nt*nx,nf)))
                  }
                  sampleEtaA1_tf_fun = tf_function(fun) #tf_function
                  assign(sprintf("updateEta_kronecker_TF_%.3d",r), sampleEtaA1_tf_fun, envir=parent.frame())
               } else{
                  sampleEtaA1_tf_fun = get(sprintf("updateEta_kronecker_TF_%.3d",r), envir=parent.frame())
               }
               m0Mat = tf$reshape(tf$constant(aperm(array(M0Full,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64), ic(nt,nx*nf))
               epsMat = tf$reshape(tf$constant(aperm(array(epsFullRand,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64), ic(nt,nx*nf))
               EtaFullTf = sampleEtaA1_tf_fun(tf$constant(alpha-1,tf$int32),tf$constant(LamiDLam,tf$float64),m0Mat,epsMat)
               EtaFull = EtaFullTf$numpy()
               if(np[r] == nx*nt){
                  eta = EtaFull
               } else{
                  eta = EtaFull[indKronObs,]
               }
            } else if(rL[[r]]$etaMethod == "TF_krylov"){
               if(np[r] == nx*nt){
                  EtaPrevFull = Eta[[r]]
               } else{
                  EtaPrevFull = matrix(0,nt*nx,nf)
                  EtaPrevFull[indKronObs,] = Eta[[r]]
               }
               if(!exists(sprintf("updateEta_kronecker_TF_krylov_%.3d",r), envir=parent.frame())){
                  fun = function(alpha,LamiDLam,m0Array,epsArray,EtaPrevArray,randMult){
                     cgIterN = rL[[r]]$cgIterN
                     iKtArray = tf$gather(rLPar[[r]][[1]]$iWStack, alpha[,1])
                     eKtMat = tf$gather(rLPar[[r]][[1]]$eWStack, alpha[,1])
                     vKtArray = tf$gather(rLPar[[r]][[1]]$vWStack, alpha[,1])
                     iKsArray = tf$gather(rLPar[[r]][[2]]$iWStack, alpha[,2])
                     eKsMat = tf$gather(rLPar[[r]][[2]]$eWStack, alpha[,2])
                     vKsArray = tf$gather(rLPar[[r]][[2]]$vWStack, alpha[,2])
                     eKArray = tf$einsum("ha,hc->hac",eKtMat,eKsMat)
                     LamiDLamKron = tf$reshape(tf$scatter_nd(matrix(ic(indKronObs-1)), LamiDLam, ic(nt*nx,nf,nf)), ic(nt,nx,nf,nf))
                     LamiDLamDiagMean = tf$reduce_mean(tfla$diag_part(LamiDLamKron), ic(0,1))

                     x = EtaPrevArray
                     iKx = tf$transpose(tfla$tridiagonal_matmul(iKtArray, tf$matmul(tf$transpose(x,ic(2,0,1)), iKsArray)), ic(1,2,0))
                     Bx = tf$einsum("txh,txhk->txk", x, LamiDLamKron)
                     Ax = iKx + Bx
                     r = m0Array - Ax
                     tmp1 = tf$einsum("hab,ach,hcd->hbd", vKtArray,r,vKsArray)
                     tmp2 = tmp1 * (eKArray^-1 + LamiDLamDiagMean[,NULL,NULL])^-1
                     z = tf$einsum("hab,hbd,hcd->ach", vKtArray,tmp2,vKsArray)
                     p = z
                     rTz = tf$reduce_sum(r*z)
                     A1aLoop1Cond = function(cgIter,x,p,r,z,rTz) tfm$less(cgIter, ic(cgIterN))
                     A1aLoop1Body = function(cgIter,x,p,r,z,rTz){
                        iKp = tf$transpose(tfla$tridiagonal_matmul(iKtArray, tf$matmul(tf$transpose(p,ic(2,0,1)), iKsArray)), ic(1,2,0))
                        Bp = tf$einsum("txh,txhk->txk", p, LamiDLamKron)
                        Ap = iKp + Bp
                        pTAp = tf$reduce_sum(p*Ap)
                        a = tfm$divide_no_nan(rTz, pTAp)
                        x = x + a*p
                        rNew = r - a*Ap
                        tmp1 = tf$einsum("hab,ach,hcd->hbd", vKtArray,rNew,vKsArray)
                        tmp2 = tmp1 * (eKArray^-1 + LamiDLamDiagMean[,NULL,NULL])^-1
                        zNew = tf$einsum("hab,hbd,hcd->ach", vKtArray,tmp2,vKsArray)
                        rTrNew = tf$reduce_sum(rNew^2)
                        # tf$print(rTrNew)
                        rTzNew = tf$reduce_sum(rNew*zNew)
                        b = tfm$divide_no_nan(rTzNew, rTz)
                        p = zNew + b*p
                        return(c(cgIter+ic(1),x,p,rNew,zNew,rTzNew))
                     }
                     A1aLoop1Init = c(tf$constant(ic(0),tf$int32),x,p,r,z,rTz)
                     A1aLoop1Res = tf$while_loop(A1aLoop1Cond, A1aLoop1Body, A1aLoop1Init)
                     EtaFullMean = tf$reshape(A1aLoop1Res[[2]],ic(nt*nx,nf))

                     epsArrayNorm = tf$sqrt(tf$reduce_sum(epsArray^2))
                     V = tf$scatter_nd(tf$zeros(ic(1,1),tf$int32), tf$reshape(epsArray,ic(1,nt,nx,nf))/epsArrayNorm, ic(cgIterN+1,nt,nx,nf))
                     alpha = tf$zeros(ic(cgIterN), tf$float64)
                     beta = tf$zeros(ic(cgIterN+1), tf$float64)
                     A1aLoop2Cond = function(j,V,alpha,beta) tfm$less(j, ic(cgIterN))
                     A1aLoop2Body = function(j,V,alpha,beta){
                        tmp1 = tf$einsum("hab,ach,hcd->hbd", vKtArray,V[j,,,],vKsArray)
                        tmp2 = tmp1 * (eKArray^-1 + LamiDLamDiagMean[,NULL,NULL])^-0.5
                        Gvj = tf$einsum("hab,hbd,hcd->ach", vKtArray,tmp2,vKsArray)
                        iKGvj = tf$transpose(tfla$tridiagonal_matmul(iKtArray, tf$matmul(tf$transpose(Gvj,ic(2,0,1)), iKsArray)), ic(1,2,0))
                        BGvj = tf$einsum("txh,txhk->txk", Gvj, LamiDLamKron)
                        AGvj = iKGvj + BGvj
                        tmp1 = tf$einsum("hab,ach,hcd->hbd", vKtArray,AGvj,vKsArray)
                        tmp2 = tmp1 * (eKArray^-1 + LamiDLamDiagMean[,NULL,NULL])^-0.5
                        GAGvj = tf$einsum("hab,hbd,hcd->ach", vKtArray,tmp2,vKsArray)
                        v = GAGvj
                        v = tf$cond(j>ic(0), function() v-beta[j]*V[j-ic(1),,,], function() v)
                        alpha = tf$tensor_scatter_nd_add(alpha, j*tf$ones(ic(1,1),tf$int32), tf$reduce_sum(V[j,,,]*v)[NULL])
                        v = v - alpha[j]*V[j,,,]
                        beta = tf$tensor_scatter_nd_add(beta, (j+ic(1))*tf$ones(ic(1,1),tf$int32), tf$sqrt(tf$reduce_sum(v^2))[NULL])
                        V = tf$tensor_scatter_nd_add(V, (j+ic(1))*tf$ones(ic(1,1),tf$int32), (v/beta[j+ic(1)])[NULL,,,])
                        return(c(j+ic(1),V,alpha,beta))
                     }
                     A1aLoop2Init = c(tf$constant(ic(0),tf$int32),V,alpha,beta)
                     A1aLoop2Res = tf$while_loop(A1aLoop2Cond, A1aLoop2Body, A1aLoop2Init)
                     V=A1aLoop2Res[[2]]; alpha=A1aLoop2Res[[3]]; beta=A1aLoop2Res[[4]]
                     V = tf$reshape(tf$transpose(V[1:cgIterN,,,], ic(1,2,3,0)), ic(nt*nx*nf,cgIterN))
                     #evT = tf.linalg.eigh_tridiagonal(....)# in v2.6.0+ m,m,c(-1,0,1),list(beta[2:m],alpha,beta[2:m]))
                     TM = tfla$LinearOperatorTridiag(list(tf$concat(list(beta[2:cgIterN],tf$zeros(ic(1),tf$float64)),ic(0)), alpha, beta[1:cgIterN]), diagonals_format='sequence')$to_dense()
                     evT = tfla$eigh(TM)
                     eT=evT[[1]]; vT=evT[[2]]
                     e1 = tf$scatter_nd(tf$zeros(ic(1,1),tf$int32), tf$ones(ic(1),tf$float64), ic(cgIterN)*tf$ones(ic(1),tf$int32))
                     iT05e1 = tf$matmul(vT, (eT^-0.5)[,NULL]*tf$matmul(vT,e1[,NULL],transpose_a=TRUE))
                     ViT05e1 = tf$reshape(tf$matmul(V,iT05e1), ic(nt,nx,nf))
                     tmp1 = tf$einsum("hab,ach,hcd->hbd", vKtArray,ViT05e1,vKsArray)
                     # tmp2 = tmp1 * (eKArray^-1 + frac*tfla$diag_part(LamiDLam)[,NULL,NULL])^-0.5
                     tmp2 = tmp1 * (eKArray^-1 + LamiDLamDiagMean[,NULL,NULL])^-0.5
                     GViT05e1 = tf$einsum("hab,hbd,hcd->ach", vKtArray,tmp2,vKsArray)
                     EtaFullRand = epsArrayNorm*tf$reshape(GViT05e1, ic(nt*nx,nf))
                     EtaFull = EtaFullMean + randMult*EtaFullRand
                     return(EtaFull)
                  }
                  sampleEtaA1a_tf_fun = tf_function(fun) #tf_function
                  assign(sprintf("updateEta_kronecker_TF_krylov_%.3d",r), sampleEtaA1a_tf_fun, envir=parent.frame())
               } else{
                  sampleEtaA1a_tf_fun = get(sprintf("updateEta_kronecker_TF_krylov_%.3d",r), envir=parent.frame())
               }
               m0Array = tf$constant(aperm(array(M0Full,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64)
               epsArray = tf$constant(aperm(array(epsFullRand,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64)
               EtaPrevArray = tf$constant(aperm(array(EtaPrevFull,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64)
               EtaFullTf = sampleEtaA1a_tf_fun(tf$constant(alpha-1,tf$int32),LamiDLam,m0Array,epsArray,EtaPrevArray,tf$constant(1,tf$float64))


               # KtArray = tf$gather(rLPar[[r]][[1]]$WStack, ic(alpha[,1]-1))
               # iKtArray = tfla$cholesky_solve(tfla$cholesky(KtArray), tf$eye(ic(nt),batch_shape=ic(nf)*tf$ones(ic(1),tf$int32),dtype=tf$float64))
               # iKsArray = tf$gather(rLPar[[r]][[2]]$iWStack, ic(alpha[,2]-1))
               # iKstArray = tf$reshape(iKtArray[,,NULL,,NULL] * iKsArray[,NULL,,NULL,], ic(nf,nt*nx,nt*nx))
               # iKst = tf$reshape(tf$transpose(tfla$diag(tf$transpose(iKstArray, ic(1,2,0))), ic(0,2,1,3)), ic(nt*nx*nf,nt*nx*nf))
               # LamiDLamKron = tf$scatter_nd(matrix(ic(indKronObs-1)), LamiDLam, ic(nt*nx,nf,nf))
               # W = iKst + tf$reshape(tf$transpose(tfla$diag(tf$transpose(LamiDLamKron, ic(1,2,0))), ic(2,0,3,1)), ic(nt*nx*nf,nt*nx*nf))
               # evW = tfla$eigh(W)
               # eW = evW[[1]]; vW=evW[[2]]
               # iWm0Array = tf$einsum("ij,j,kj,k->i",vW,eW^-1,vW,tf$reshape(m0Array,ic(nt*nx*nf)))
               # iW05Eps = tf$einsum("ij,j,kj,k->i",vW,eW^-0.5,vW,tf$reshape(epsArray,ic(nt*nx*nf)))
               # EtaFullTfExact = iWm0Array + 1*iW05Eps
               # plot(tf$reshape(EtaFullTf,ic(-1))$numpy(), EtaFullTfExact$numpy())


               EtaFull = EtaFullTf$numpy()
               if(np[r] == nx*nt){
                  eta = EtaFull
               } else{
                  eta = EtaFull[indKronObs,]
               }
            }
         }
      }
      rownames(eta)=rnames
      Eta[[r]] = eta
      EtaFullList[[r]] = EtaFull
      if(r < nr){
         if(class(rL[[r]])[1]=="HmscRandomLevel"){
            if(rL[[r]]$xDim == 0){
               LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
            } else{
               LRan[[r]] = matrix(0,ny,ns)
               for(k in 1:rL[[r]]$xDim)
                  LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),r]) %*% Lambda[[r]][,,r]
            }
         } else{
            LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
         }
      }
   }
   return(list(Eta, EtaFullList))
}




# fun = function(alpha,LamiDLam,m0Array,epsArray,EtaPrevArray,randMult){
#    cgIterN = rL[[r]]$cgIterN
#    iKtArray = tf$gather(rLPar[[r]][[1]]$iWStack, alpha[,1])
#    iKsArray = tf$gather(rLPar[[r]][[2]]$iWStack, alpha[,2])
#    LamiDLamKron = tf$reshape(tf$scatter_nd(matrix(ic(indKronObs-1)), LamiDLam, ic(nt*nx,nf,nf)), ic(nt,nx,nf,nf))
#    LamiDLamDiagMean = tf$reduce_mean(tfla$diag_part(LamiDLamKron), ic(0,1))
#    x = EtaPrevArray
#    iKx = tf$transpose(tfla$tridiagonal_matmul(iKtArray, tf$matmul(tf$transpose(x,ic(2,0,1)), iKsArray)), ic(1,2,0))
#    Bx = tf$einsum("txh,txhk->txk", x, LamiDLamKron)
#    Ax = iKx + Bx
#    r = m0Array - Ax
#    p = r
#    rTr = tf$reduce_sum(r*r)
#    A1aLoop1Cond = function(cgIter,x,p,r,rTr) tfm$less(cgIter, ic(cgIterN))
#    A1aLoop1Body = function(cgIter,x,p,r,rTr){
#       iKp = tf$transpose(tfla$tridiagonal_matmul(iKtArray, tf$matmul(tf$transpose(p,ic(2,0,1)), iKsArray)), ic(1,2,0))
#       Bp = tf$einsum("txh,txhk->txk", p, LamiDLamKron)
#       Ap = iKp + Bp
#       pTAp = tf$reduce_sum(p*Ap)
#       a = tfm$divide_no_nan(rTr, pTAp)
#       x = x + a*p
#       rNew = r - a*Ap
#       rTrNew = tf$reduce_sum(rNew^2)
#       tf$print(rTrNew)
#       EPS = 1e-12
#       b = tfm$divide_no_nan(rTrNew, rTr)
#       p = rNew + b*p
#       return(c(cgIter+ic(1),x,p,rNew,rTrNew))
#    }
#    A1aLoop1Init = c(tf$constant(ic(0),tf$int32),x,p,r,rTr)
#    A1aLoop1Res = tf$while_loop(A1aLoop1Cond, A1aLoop1Body, A1aLoop1Init)
#    EtaFullMean = tf$reshape(A1aLoop1Res[[2]],ic(nt*nx,nf))
#    EtaFull = EtaFullMean
#    return(EtaFull)
# }


# if(!exists("updateEta_kronecker_A1a_tf", envir=parent.frame())){
#    fun = function(alpha,LamiDLam,m0Array,epsArray,EtaPrevArray){
#       cgIterN = rL[[r]]$cgIterN
#       tfla = tf$linalg
#       tfm = tf$math
#       frac = tf$constant(sum(numKronObsMat) / (nx*nt), tf$float64)
#       # indL = 2 + (nt+1)*(0:(nt-2))
#       # indD = 1 + (nt+1)*(0:(nt-1))
#       # indU = (nt+1)*(1:(nt-1))
#       # getTridiagonal = function(A) t(matrix(c(A[indU],0,A[indD],0,A[indL]),nrow(A),3))
#       # iKtStack = tf$cast(tf$stack(lapply(rLPar[[r]][[1]]$iWg, getTridiagonal)), tf$float64)
#       # eKtStack = tf$cast(tf$stack(rLPar[[r]][[1]]$eWg), tf$float64)
#       # vKtStack = tf$cast(tf$stack(rLPar[[r]][[1]]$vWg), tf$float64)
#       # iKsStack = tf$cast(tf$stack(lapply(rLPar[[r]][[2]]$iWg, as.matrix)), tf$float64)
#       # eKsStack = tf$cast(tf$stack(rLPar[[r]][[2]]$eWg), tf$float64)
#       # vKsStack = tf$cast(tf$stack(rLPar[[r]][[2]]$vWg), tf$float64)
#       # iKtArray = tf$gather(iKtStack, alpha[,1])
#       # eKtMat = tf$gather(eKtStack, alpha[,1])
#       # vKtArray = tf$gather(vKtStack, alpha[,1])
#       # iKsArray = tf$gather(iKsStack, alpha[,2])
#       # eKsMat = tf$gather(eKsStack, alpha[,2])
#       # vKsArray = tf$gather(vKsStack, alpha[,2])
#       iKtArray = tf$gather(rLPar[[r]][[1]]$iWStack, alpha[,1])
#       eKtMat = tf$gather(rLPar[[r]][[1]]$eWStack, alpha[,1])
#       vKtArray = tf$gather(rLPar[[r]][[1]]$vWStack, alpha[,1])
#       iKsArray = tf$gather(rLPar[[r]][[2]]$iWStack, alpha[,2])
#       eKsMat = tf$gather(rLPar[[r]][[2]]$eWStack, alpha[,2])
#       vKsArray = tf$gather(rLPar[[r]][[2]]$vWStack, alpha[,2])
#    }
# }





# fun = function(iKtArray,iKsArray,m0Mat,LamiDLam,epsMat,alphaTFlag0){
#    conserveMemoryFlag = rL[[r]]$conserveMemoryFlag
#    tfla = tf$linalg
#    tfm = tf$math
#    it = tf$constant(ic(0),tf$int32)
#    diagBlockiKtVec = iKtArray[,2,it]
#    A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
#    A = A + tf$reshape(tf$transpose(tfla$diag(LamiDLam[,,NULL]*tf$constant(numKronObsMat,tf$float64)[,it]), ic(2,0,3,1)), ic(nf*nx,nf*nx))
#    L = tfla$cholesky(A)
#    iLm0 = tf$scatter_nd(it[NULL,NULL], tfla$triangular_solve(L, m0Mat[it,,NULL])[NULL,,1],ic(nt,nx*nf))
#    A1Loop1Cond = function(it,L,C,iLm0) tfm$less(it, ic(nt))
#    conserveMemoryFlagTmp = conserveMemoryFlag
#    conserveMemoryFlag = FALSE
#    if(conserveMemoryFlag==FALSE){ # saving intermediate results for faster backward loop
#       LArray = tf$scatter_nd(it[NULL,NULL], L[NULL,,], ic(nt,nx*nf,nx*nf))
#       CArray = tf$zeros(ic(nt-1,nx*nf,nx*nf), tf$float64)
#    }
#    A1Loop1Body = function(it,L,C,iLm0){
#       if(conserveMemoryFlag==FALSE){
#          LArray = L; CArray = C
#       }
#       diagBlockiKtVec = iKtArray[,2,it]
#       A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
#       A = A + tf$reshape(tf$transpose(tfla$diag(LamiDLam[,,NULL]*tf$constant(numKronObsMat,tf$float64)[,it]), ic(2,0,3,1)), ic(nf*nx,nf*nx))
#       offdiagBlockiKtVec = iKtArray[,3,it]
#       B = tf$reshape(tf$transpose(tfla$diag(offdiagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
#       if(conserveMemoryFlag==FALSE){
#          C = tf$transpose(tfla$triangular_solve(LArray[it-ic(1),,], B))
#       } else{
#          C = tf$transpose(tfla$triangular_solve(L, B))
#       }
#       L = tfla$cholesky(A-tf$matmul(C,C,transpose_b=TRUE))
#       v = m0Mat[it,] - tf$squeeze(tf$matmul(C,iLm0[it-ic(1),,NULL]),-1)
#       iLm0 = tf$tensor_scatter_nd_add(iLm0, it[NULL,NULL], tfla$triangular_solve(L, v[,NULL])[NULL,,1])
#       if(conserveMemoryFlag==FALSE){
#          CArray = tf$tensor_scatter_nd_add(CArray, (it-ic(1))[NULL,NULL], C[NULL,,])
#          LArray = tf$tensor_scatter_nd_add(LArray, it[NULL,NULL], L[NULL,,])
#          return(c(it+ic(1),LArray,CArray,iLm0))
#       } else{
#          return(c(it+ic(1),L,C,iLm0))
#       }
#    }
#    if(conserveMemoryFlag==FALSE){
#       A1Loop1Init = c(tf$constant(ic(1),tf$int32),LArray,CArray,iLm0)
#       A1Loop1Res = tf$while_loop(A1Loop1Cond, A1Loop1Body, A1Loop1Init)
#       LArray = A1Loop1Res[[2]]; CArray = A1Loop1Res[[3]]; iLm0 = A1Loop1Res[[4]]
#       L = LArray[nt,,]
#    } else{
#       A1Loop1Init = c(tf$constant(ic(1),tf$int32),L,tf$zeros(ic(nx*nf,nx*nf),tf$float64),iLm0)
#       A1Loop1Res = tf$while_loop(A1Loop1Cond, A1Loop1Body, A1Loop1Init)
#       L = A1Loop1Res[[2]]; iLm0 = A1Loop1Res[[4]]
#    }
#
#    conserveMemoryFlag = conserveMemoryFlagTmp
#    it = tf$constant(ic(nt-1),tf$int32)
#    v = iLm0[it,] + epsMat[it,]
#    if(conserveMemoryFlag==FALSE){
#       iLTiLm0 = tf$scatter_nd(it[NULL,NULL], tfla$triangular_solve(LArray[it,,], v[,NULL], adjoint=TRUE)[NULL,,1], ic(nt,nx*nf))
#       A1Loop2Cond = function(it,iLTiLm0) tfm$greater_equal(it, ic(0))
#       A1Loop2Body = function(it,iLTiLm0){
#          v = iLm0[it,] + epsMat[it,] - tf$squeeze(tf$matmul(CArray[it,,],iLTiLm0[it+ic(1),,NULL],transpose_a=TRUE),-1)
#          iLTiLm0 = tf$tensor_scatter_nd_add(iLTiLm0, it[NULL,NULL], tfla$triangular_solve(LArray[it,,], v[,NULL], adjoint=TRUE)[NULL,,1])
#          return(c(it-ic(1),iLTiLm0))
#       }
#       A1Loop2Init = c(it-ic(1),iLTiLm0)
#       A1Loop2Res = tf$while_loop(A1Loop2Cond, A1Loop2Body, A1Loop2Init)
#       iLTiLm0 = A1Loop2Res[[2]]
#    } else{ #TODO currently this remains numerically unstable, probably due to accumulation of small errors
#       mult0 = tf$reshape(tf$tile(alphaTFlag0[NULL,], ic(nx,1)), ic(nx*nf))
#       iLTiLm0 = tf$scatter_nd(it[NULL,NULL], tfla$triangular_solve(L, v[,NULL], adjoint=TRUE)[NULL,,1], ic(nt,nx*nf))
#       A1Loop2Cond = function(it,L,iLTiLm0) tfm$greater_equal(it, ic(0))
#       A1Loop2Body = function(it,L,iLTiLm0){
#          tf$print(it)
#          # tf$print("L1-1",L)
#          # tf$print("L1-2",LArray[it+ic(1),,])
#          diagBlockiKtVec = iKtArray[,2,it+ic(1)]
#          A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
#          A = A + tf$reshape(tf$transpose(tfla$diag(LamiDLam[,,NULL]*tf$constant(numKronObsMat,tf$float64)[,it+ic(1)]), ic(2,0,3,1)), ic(nf*nx,nf*nx))
#          offdiagBlockiKtVec = iKtArray[,3,it+ic(1)]
#          B = tf$reshape(tf$transpose(tfla$diag(offdiagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
#          # B = (mult0-1)[,NULL]*B*(mult0-1) + tfla$diag(mult0)
#          H = A - tf$matmul(L,L,transpose_b=TRUE)
#          # H = (mult0-1)[,NULL]*H*(mult0-1) + tfla$diag(mult0)
#          tf$print("H-1",H)
#          tf$print("H-2",tf$matmul(CArray[it,,],CArray[it,,],transpose_b=TRUE))
#          LH = tfla$cholesky(H, name="LH")
#          iLHB = tfla$triangular_solve(LH, B, name="iLHB")
#          QRList = tfla$qr(iLHB)
#          # the diagonal of QRList$r is not necessarily positive
#          # signDiagRVec = tf$sign(tfla$diag_part(QRList$r))
#          L = tf$transpose(QRList$r)
#          C = tf$matmul(LH, QRList$q)
#          # BiHB = tf$matmul(iLHB,iLHB,transpose_a=TRUE)
#          # L = tfla$cholesky(BiHB, name="L")
#          # C = tf$transpose(tfla$triangular_solve(L, B, name="C"))
#          tf$print("L-1",L)
#          tf$print("L-2",LArray[it,,])
#          tf$print("C-1",C)
#          tf$print("C-2",CArray[it,,])
#          v = iLm0[it,] + epsMat[it,] - tf$squeeze(tf$matmul(C,iLTiLm0[it+ic(1),,NULL],transpose_a=TRUE),-1)
#          iLTiLm0 = tf$tensor_scatter_nd_add(iLTiLm0, it[NULL,NULL], tfla$triangular_solve(L, v[,NULL], adjoint=TRUE, name="iLTiLm0")[NULL,,1])
#          return(c(it-ic(1),L,iLTiLm0))
#       }
#       A1Loop2Init = c(it-ic(1),L,iLTiLm0)
#       A1Loop2Res = tf$while_loop(A1Loop2Cond, A1Loop2Body, A1Loop2Init)
#       iLTiLm0 = A1Loop2Res[[3]]
#    }
#    return(tf$reshape(iLTiLm0, ic(nt*nx,nf)))
# }

