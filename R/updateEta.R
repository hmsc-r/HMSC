#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix chol
#' @importFrom SparseM chol backsolve
#' @importFrom plyr mdply
#'
updateEta = function(Y,Z,Beta,iSigma,Eta,Lambda,Alpha, rLPar, X,Pi,dfPi,rL){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   Yx = !is.na(Y)
   ic = function(...){
      return(as.integer(c(...)))
   }

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
         invSigLam = lambda*matrix(iSigma,nf,ns,byrow=TRUE)
         LamInvSigLam = tcrossprod(lambda*matrix(sqrt(iSigma),nf,ns,byrow=TRUE))
         if(np[r] == ny){
            numObsVec = rep(1,ny)
            fS = tcrossprod(S[order(lPi),,drop=FALSE], invSigLam)
         }else{
            P = sparseMatrix(i=1:ny,j=lPi)
            numObsVec = Matrix::colSums(P)
            fS = Matrix::tcrossprod(Matrix::crossprod(P,S), invSigLam)
         }
         epsRand = rnorm(nf*np[r])
         dfPiElemLevelList = vector("list", length(rL[[r]]$rLList))
         npElemVec = rep(NA, length(rL[[r]]$rLList))
         for(l in seq_len(length(rL[[r]]$rLList))){
            dfPiElem = as.factor(unlist(lapply(strsplit(as.character(m$dfPi[,r]), rL[[r]]$sepStr), function(a) a[l])))
            dfPiElemLevelList[[l]] = unique(dfPiElem)
            npElemVec[l] = length(dfPiElemLevelList[[l]])
         }
         dfTmp = expand.grid(rev(dfPiElemLevelList))[,rev(1:length(rL[[r]]$rLList))]
         # allUnits = factor(mdply(dfTmp,paste,sep=rL[[r]]$sepStr,.expand=FALSE)[,2])
         allUnits = factor(do.call(function(...) paste(..., sep=rL[[r]]$sepStr),dfTmp))
         indKronObs = as.numeric(factor(m$dfPi[,r], levels=levels(allUnits)))
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
               tmp1 = Matrix::kronecker(Diagonal(x=numObsVec), LamInvSigLam)
            }
            iUEta = iWs + tmp1
            R = Matrix::chol(iUEta)
            tmp2 = Matrix::solve(t(R), as.vector(t(fS)))
            feta = Matrix::solve(R, tmp2+epsRand)
            eta = matrix(feta,np[r],nf,byrow=TRUE)
         } else if(rL[[r]]$etaMethod == "TF"){
            nt = npElemVec[1]
            nx = npElemVec[2]
            if(np[r] == nx*nt){
               m0Full = as.vector(as.matrix(fS))
               numKronObsMat = matrix(numObsVec,nx,nt)
               epsFullRand = epsRand
            } else{
               M0Full = matrix(0,nx*nt,nf)
               M0Full[indKronObs,] = as.matrix(fS)
               m0Full = as.vector(M0Full)
               numKronObsMat = matrix(0,nx,nt)
               numKronObsMat[indKronObs] = numObsVec
               epsFullRand = rnorm(nx*nt*nf)
            }
            if(!exists("updateEta_kronecker_A1_tf", envir=parent.frame())){
               fun = function(iKtArray,iKsArray,m0Mat,LamiDLam,epsMat){
                  tfla = tf$linalg
                  tfm = tf$math
                  LArray = tf$zeros(ic(nt,nx*nf,nx*nf), tf$float64)
                  CArray = tf$zeros(ic(nt-1,nx*nf,nx*nf), tf$float64)
                  iLm0 = tf$zeros(ic(nt,nx*nf), tf$float64)
                  it = 1
                  diagBlockiKtVec = iKtArray[,2,it]
                  A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
                  A = A + tf$reshape(tf$transpose(tfla$diag(LamiDLam[,,NULL]*tf$constant(numKronObsMat,tf$float64)[,it]), ic(2,0,3,1)), ic(nf*nx,nf*nx))
                  L = tfla$cholesky(A)
                  LArray = tf$tensor_scatter_nd_add(LArray, tf$constant(ic(it-1),tf$int32)[NULL,NULL], L[NULL,,])
                  iLm0 = tf$tensor_scatter_nd_add(iLm0, tf$constant(ic(it-1),tf$int32)[NULL,NULL], tfla$triangular_solve(L, m0Mat[it,,NULL])[NULL,,1])
                  A1Loop1Cond = function(it,LArray,CArray,iLm0) tfm$less_equal(it, as.integer(nt))
                  A1Loop1Body = function(it,LArray,CArray,iLm0){
                     diagBlockiKtVec = iKtArray[,2,it-ic(1)]
                     A = tf$reshape(tf$transpose(tfla$diag(diagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
                     A = A + tf$reshape(tf$transpose(tfla$diag(LamiDLam[,,NULL]*tf$constant(numKronObsMat,tf$float64)[,it-ic(1)]), ic(2,0,3,1)), ic(nf*nx,nf*nx))
                     offdiagBlockiKtVec = iKtArray[,3,it-ic(1)]
                     B = tf$reshape(tf$transpose(tfla$diag(offdiagBlockiKtVec*tf$transpose(iKsArray,ic(1,2,0))),ic(0,2,1,3)), ic(nx*nf,nx*nf))
                     C = tf$transpose(tfla$triangular_solve(LArray[it-ic(2),,], B))
                     CArray = tf$tensor_scatter_nd_add(CArray, (it-ic(2))[NULL,NULL], C[NULL,,])
                     L = tfla$cholesky(A-tf$matmul(C,C,transpose_b=TRUE))
                     LArray = tf$tensor_scatter_nd_add(LArray, (it-ic(1))[NULL,NULL], L[NULL,,])
                     v = m0Mat[it-ic(1),] - tf$squeeze(tf$matmul(C,iLm0[it-ic(2),,NULL]),-1)
                     iLm0 = tf$tensor_scatter_nd_add(iLm0, (it-ic(1))[NULL,NULL], tfla$triangular_solve(L, v[,NULL])[NULL,,1])
                     return(c(it+ic(1),LArray,CArray,iLm0))
                  }
                  A1Loop1Init = c(tf$constant(as.integer(2),tf$int32),LArray,CArray,iLm0)
                  A1Loop1Res = tf$while_loop(A1Loop1Cond, A1Loop1Body, A1Loop1Init)
                  LArray = A1Loop1Res[[2]]; CArray = A1Loop1Res[[3]]; iLm0 = A1Loop1Res[[4]]
                  iLTiLm0 = tf$zeros(ic(nt,nx*nf), tf$float64)
                  it = nt
                  v = iLm0[it,] + epsMat[it,]
                  iLTiLm0 = tf$tensor_scatter_nd_add(iLTiLm0, tf$constant(ic(it-1),tf$int32)[NULL,NULL], tfla$triangular_solve(LArray[it,,], v[,NULL], adjoint=TRUE)[NULL,,1])
                  A1Loop2Cond = function(it,iLTiLm0) tfm$greater_equal(it, as.integer(1))
                  A1Loop2Body = function(it,iLTiLm0){
                     v = iLm0[it-ic(1),] + epsMat[it-ic(1),] - tf$squeeze(tf$matmul(CArray[it-ic(1),,],iLTiLm0[it+ic(0),,NULL],transpose_a=TRUE),-1)
                     iLTiLm0 = tf$tensor_scatter_nd_add(iLTiLm0, (it-ic(1))[NULL,NULL], tfla$triangular_solve(LArray[it-ic(1),,], v[,NULL], adjoint=TRUE)[NULL,,1])
                     return(c(it-ic(1),iLTiLm0))
                  }
                  A1Loop2Init = c(tf$constant(as.integer(nt-1),tf$int32),iLTiLm0)
                  A1Loop2Res = tf$while_loop(A1Loop2Cond, A1Loop2Body, A1Loop2Init)
                  iLTiLm0 = A1Loop2Res[[2]]
                  return(tf$reshape(iLTiLm0, ic(nt*nx,nf)))
               }
               sampleEtaA1_tf_fun = tf_function(fun)
               assign("updateEta_kronecker_A1_tf", sampleEtaA1_tf_fun, envir=parent.frame())
            } else{
               sampleEtaA1_tf_fun = get("updateEta_kronecker_A1_tf", envir=parent.frame())
            }
            m0Mat = tf$reshape(tf$constant(aperm(array(m0Full,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64), ic(nt,nx*nf))
            epsMat = tf$reshape(tf$constant(aperm(array(epsFullRand,c(nx,nt,nf)),c(2,1,3)), dtype=tf$float64), ic(nt,nx*nf))
            # iKtArray = tf$gather(iKtGrid, ic(alpha[h,1]-1))
            # iKsArray = tf$gather(iKsGrid, ic(alpha[h,2]-1))
            # iKtGrid = tf$stack(lapply(iKtGridList, function(A) t(matrix(c(A[indU],0,A[indD],0,A[indL]),nrow(A),3))))
            # iKsGrid = tf$stack(iKsGridList)
            indL = 2 + (nt+1)*(0:(nt-2))
            indD = 1 + (nt+1)*(0:(nt-1))
            indU = (nt+1)*(1:(nt-1))
            getTridiagonal = function(A) t(matrix(c(A[indU],0,A[indD],0,A[indL]),nrow(A),3))
            iKtArray = tf$stack(lapply(rLPar[[r]][[1]]$iWg[alpha[,1]], getTridiagonal))
            iKsArray = tf$stack(lapply(rLPar[[r]][[2]]$iWg[alpha[,2]], as.matrix))
            LamiDLam = tf$constant(LamInvSigLam, dtype=tf$float64)
            EtaFullTf = sampleEtaA1_tf_fun(iKtArray,iKsArray,m0Mat,LamiDLam,epsMat)
            EtaFull = EtaFullTf$numpy()
         }
         if(np[r] == nx*nt){
            eta = EtaFull
         } else{
            eta = EtaFull[indKronObs,]
         }
      }
      rownames(eta)=rnames
      Eta[[r]] = eta
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
   return(Eta)
}

