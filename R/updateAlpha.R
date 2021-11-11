#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#' @importFrom matrixStats colMaxs logSumExp
#' @importFrom pracma ndims
#' @importFrom plyr mdply
#' @importFrom tensorflow tf
#'
updateAlpha = function(Z,Beta,iSigma,Eta,EtaFull,Alpha,Lambda, rLPar, X,Pi,dfPi,rL){
   nr = length(rL)
   ny = nrow(Z)
   ns = ncol(Z)
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
      if(nr > 1){
         S = Z - (LFix + Reduce("+", LRan[setdiff(1:nr, r)]))
      } else{
         S = Z - LFix
      }
      eta = Eta[[r]]
      lambda = Lambda[[r]]
      np = nrow(eta)
      nf = ncol(eta)
      if(inherits(rL[[r]],"HmscRandomLevel",TRUE)==1){
         Alpha[[r]] = rep(1, nf)
      } else if(inherits(rL[[r]],"HmscSpatialRandomLevel",TRUE)==1){
         iWg = rLPar[[r]]$iWg
         RiWg = rLPar[[r]]$RiWg
         detWg = rLPar[[r]]$detWg
         alphapw = rL[[r]]$alphapw
         gN = nrow(alphapw)

         Alpha[[r]] = rep(NA, nf)
         if(rL[[r]]$spatialMethod =="Full" || rL[[r]]$spatialMethod=="NNGP"){
            tmpMat = matrix(NA,nrow=gN,ncol=nf)
            for(g in seq_len(gN)){
               tmp1 = RiWg[[g]]%*%eta
               tmpMat[g,] = colSums(tmp1^2)
            }
            logLikeMat = log(alphapw[,2]) - 0.5*detWg - 0.5*tmpMat
         } else if(rL[[r]]$spatialMethod =="GPP"){ #TODO:GT this part shall be carefully revised, I am not sure of its correctness
            iFg = rLPar[[r]]$iFg
            nK = nrow(iFg)
            idDW12g = rLPar[[r]]$idDW12g
            idDg = rLPar[[r]]$idDg
            detDg = rLPar[[r]]$detDg
            tmpMat2 = array(NA,dim=c(nf,nK,gN))
            tmpMat3 = array(NA,dim=c(nf,nK,gN))
            tmpMat4 = matrix(NA,dim=c(nf,gN))
            for(g in 1:seq_len(gN)){
               tmpMat2[,,g] = crossprod(eta,idDW12g[,,g])
               tmpMat3[,,g] = tmpMat2[,,g]%*%iFg[,,g]
               tmpMat4[,g] = rowSums(tmpMat2[,,g]*tmpMat3[,,g])
            }
            tmpMat = matrix(NA,nrow=gN,ncol=nf)
            for(ag in seq_len(gN)){
               if(alphapw[ag,1] == 0){
                  tmpMat[ag,] = colSums(eta^2)
               } else {
                  tmpMat[ag,] = colSums(eta*(idDg[,ag]*eta)) - tmpMat4[,ag]
               }
            }
            logLikeMat = log(alphapw[,2]) - 0.5*detDg - 0.5*tmpMat #TODO:GT the determinant here seems to be incorrect - shall it not include parts other than diagonal
         }
         likeMat = exp(logLikeMat - matrix(colMaxs(logLikeMat),gN,nf,byrow=TRUE))
         for(h in seq_len(nf)){
            Alpha[[r]][h] = sample.int(gN, size=1, prob=likeMat[,h])
         }
      } else if(inherits(rL[[r]],"HmscKroneckerRandomLevel",TRUE)==1){
         krN = length(rL[[r]]$rLList)
         # Alpha[[r]] = matrix(NA, nf, krN)
         dfPiElemLevelList = vector("list", krN)
         npElemVec = rep(NA, krN)
         for(l in seq_len(krN)){
            dfPiElem = as.factor(unlist(lapply(strsplit(as.character(m$dfPi[,r]), "="), function(a) a[l])))
            dfPiElemLevelList[[l]] = levels(dfPiElem)
            npElemVec[l] = length(dfPiElemLevelList[[l]])
         }
         dfTmp = expand.grid(rev(dfPiElemLevelList))[,rev(1:krN)]
         # allUnits = as.factor(mdply(dfTmp, paste, sep=rL[[r]]$sepStr)[,length(dfPiElemLevelList)+1])
         allUnits = factor(do.call(function(...) paste(..., sep=rL[[r]]$sepStr),dfTmp))
         indKronObs = as.numeric(factor(levels(m$dfPi[,r]), levels=levels(allUnits)))
         gNVec = sapply(rL[[r]]$alphaPrior$alphaGridList, length)
         if(rL[[r]]$alphaMethod=="R" || rL[[r]]$alphaMethod=="R_full"){
            if(rL[[r]]$alphaMethod=="R_full"){
               etaArray = array(EtaFull[[r]],c(rev(npElemVec),nf))
               iWgEtaList = vector("list", krN)
               for(l in seq_len(krN)){
                  if(rL[[r]]$rLList[[l]]$sDim == 0){
                     iWgEtaList[[l]] = array(etaArray, c(1,dim(etaArray)))
                  } else{
                     iWg = tf$stack(lapply(rLPar[[r]][[l]]$iWg, as.matrix))
                     if(l==1){
                        iWgEta = tf$einsum("xtf,gtq->gxqf",etaArray, iWg)
                     }
                     if(l==2){
                        iWgEta = tf$einsum("gyx,xtf->gytf",iWg,etaArray)
                     }
                     # gN = gNVec[l]
                     # iWg = rLPar[[r]][[l]]$iWg
                     # iWgEtaElemList = vector("list",gN)
                     # dimInd = c(krN-l+1,c(1:ndims(etaArray))[-(krN-l+1)])
                     # dimIndOrd = order(dimInd)
                     # etaMatFlat = matrix(aperm(etaArray, dimInd), dim(etaArray)[krN-l+1])
                     # for(ag in seq_len(gN)){
                     #    iWgEtaMatFlat = iWg[[ag]]%*%etaMatFlat
                     #    iWgEtaElemList[[ag]] = array(iWgEtaMatFlat, dim(etaArray)[dimInd])
                     # }
                     # iWgEta2 = aperm(abind(iWgEtaElemList, rev.along=0), c(krN+2,dimIndOrd))
                     iWgEtaList[[l]] = iWgEta
                  }
               }
               #TODO:GT potential for >2 kronecker elements
               gN1 = gNVec[1]
               gN2 = gNVec[2]
               qFArray = tf$linalg$einsum("axtf,bxtf->abf",iWgEtaList[[1]],iWgEtaList[[2]])$numpy()
               tmp1 = matrix(rLPar[[r]][[1]]$detWg, gN1,gN2)
               tmp2 = matrix(rLPar[[r]][[2]]$detWg, gN1,gN2, byrow=TRUE)
               logDetWArray = array(npElemVec[2]*tmp1+npElemVec[1]*tmp2, c(gN1,gN2,nf))
               logPriorArray = array(log(rL[[r]]$alphaPrior$alphaProb), c(gN1,gN2,nf))
               logProbArray = logPriorArray - 0.5*logDetWArray - 0.5*qFArray #dimension-based const coefficients are omitted in this expr
               logProbArray1 = logProbArray
            } else if(rL[[r]]$alphaMethod=="R"){
               ind = as.numeric(factor(m$dfPi[,r], levels=levels(allUnits)))
               qFArray = array(NA,c(nf,gNVec))
               logDetWArray = array(NA,gNVec)
               for(gInd1 in seq_len(gNVec[1])){
                  for(gInd2 in seq_len(gNVec[2])){ #TODO:GT potential for >2 kronecker elements
                     gIndVec = c(gInd1,gInd2)
                     WfList = vector("list", krN)
                     for(l in seq_len(krN)){
                        if(rL[[r]]$rLList[[l]]$sDim == 0){
                           WfList[[l]] = Diagonal(npElemVec[l])
                        } else{
                           WfList[[l]] = rLPar[[r]][[l]]$Wg[[gIndVec[l]]]
                        }
                     }
                     WFull = WfList[[1]]
                     for(l in 1+seq_len(krN-1)){ #TODO can this be done without loop?
                        WFull = Matrix::kronecker(WFull, WfList[[l]])
                     }
                     Wf = WFull[ind,ind]
                     RWf = chol(Wf)
                     logDetWArray[matrix(gIndVec,1,krN)] = 2*sum(log(diag(RWf)))
                     qFArray[cbind(1:nf,matrix(gIndVec,nf,krN,byrow=TRUE))] = colSums(backsolve(RWf, eta, transpose=TRUE)^2)
                  }
               }
               tmp1 = array(log(rL[[r]]$alphaPrior$alphaProb),c(gNVec,nf))
               tmp2 = array(logDetWArray,c(gNVec,nf))
               tmp3 = aperm(qFArray, c(1+seq_len(krN), 1))
               logProbArray =  tmp1 - 0.5*tmp2 - 0.5*tmp3
               logProbArray2 = logProbArray
            }
            for(h in 1:nf){
               logProbVec = as.vector(logProbArray[,,h])
               logProbVec = logProbVec - logSumExp(logProbVec)
               alphaFlatInd = sample(gNVec[1]*gNVec[2], 1, prob=exp(logProbVec))
               alphaInd1 = ((alphaFlatInd-1) %% gNVec[1]) + 1
               alphaInd2 = ((alphaFlatInd-1) %/% gNVec[1]) + 1
               Alpha[[r]][h,] = c(alphaInd1,alphaInd2)
            }
         } else if(rL[[r]]$alphaMethod=="TF_full" || rL[[r]]$alphaMethod=="TF_obs" || rL[[r]]$alphaMethod=="TF_direct_krylov"){
            if(rL[[r]]$alphaMethod=="TF_full"){
               if(!exists("updateAlpha_TF_full", envir=parent.frame())){
                  fun = function(etaArray){
                     iWgEtaList = vector("list", krN)
                     for(l in seq_len(krN)){
                        # iWg = tf$cast(tf$stack(lapply(rLPar[[r]][[l]]$iWg, as.matrix)), tf$float64)
                        iWg = rLPar[[r]][[l]]$iWStack
                        if(l==1){
                           # iWgEta = tf$einsum("xtf,gtq->gxqf",etaArray, iWg)
                           iWgEta = tf$linalg$tridiagonal_matmul(tf$tile(iWg[,NULL,,],ic(1,npElemVec[2],1,1)),tf$tile(etaArray[NULL,,,],ic(gNVec[1],1,1,1)))
                        }
                        if(l==2){
                           iWgEta = tf$einsum("gyx,xtf->gytf",iWg,etaArray)
                        }
                        iWgEtaList[[l]] = iWgEta
                     }
                     gN1 = gNVec[1]
                     gN2 = gNVec[2]
                     qFArray = tf$linalg$einsum("axtf,bxtf->fab",iWgEtaList[[1]],iWgEtaList[[2]])
                     tmp1 = matrix(rLPar[[r]][[1]]$detWg, gN1,gN2)
                     tmp2 = matrix(rLPar[[r]][[2]]$detWg, gN1,gN2, byrow=TRUE)
                     logDetW = npElemVec[2]*tmp1+npElemVec[1]*tmp2
                     logPrior = log(rL[[r]]$alphaPrior$alphaProb)
                     logLikeArray = -tf$constant(0.5,tf$float64)*qFArray - 0.5*logDetW + logPrior #dimension-based const coefficients are omitted in this expr
                     logLikeMat = tf$reshape(logLikeArray, ic(nf,gNVec[1]*gNVec[2]))
                  }
                  logLike_tf_fun = tf_function(fun)
                  assign("updateAlpha_TF_full", logLike_tf_fun, envir=parent.frame())
               } else{
                  logLike_tf_fun = get("updateAlpha_TF_full", envir=parent.frame())
               }
               etaArray = tf$constant(array(EtaFull[[r]],c(rev(npElemVec),nf)), tf$float64)
               logLikeMat = logLike_tf_fun(etaArray)
               logLikeMat1 = logLikeMat
            } else if(rL[[r]]$alphaMethod=="TF_obs"){ #TODO to be checked thoughtfully
               nt = npElemVec[1]
               nx = npElemVec[2]
               P = sparseMatrix(i=1:ny,j=Pi[,r])
               numObsVec = Matrix::colSums(P)
               # Y1 = as.matrix(Matrix::crossprod(P, (S - LRan[[r]] + tf$einsum("ih,hj->hij", eta[Pi[,r],], lambda))))
               Y1 = S - LRan[[r]] + tf$einsum("ih,hj->hij", eta[Pi[,r],], lambda)
               iDY1 = Y1 * matrix(iSigma,ny,ns,byrow=TRUE)
               LambdaiDY1 = tf$einsum("hj,hij->hi",lambda,iDY1)
               LambdaiDY1 = tf$transpose(tf$scatter_nd(tf$constant(matrix(ic(Pi[,r]-1),ny,1)), tf$transpose(LambdaiDY1), tf$constant(ic(np,nf))))
               if(np == nt*nx){
                  m0Mat = LambdaiDY1
                  numKronObsMat = matrix(numObsVec,nx,nt)
               } else{
                  m0Mat = matrix(0,nf,nt*nx)
                  m0Mat[,indKronObs] = LambdaiDY1$numpy()
                  numKronObsMat = matrix(0,nx,nt)
                  numKronObsMat[indKronObs] = numObsVec
                  m0Mat = tf$constant(m0Mat, tf$float64)
               }
               m0Mat = tf$reshape(m0Mat, ic(nf,nt,nx))
               y1iDy1 = tf$reduce_sum(Y1**2 * tf$constant(iSigma,tf$float64), ic(1,2))
               logDetD = -ny*sum(log(iSigma))
               lidl = tf$einsum("hj,j,hj->h",lambda,tf$constant(iSigma,tf$float64),lambda)
               if(!exists("updateAlpha_kronecker_B2_tf", envir=parent.frame())){
                  fun = function(m0Mat,lidl){
                     tfla = tf$linalg
                     tfm = tf$math
                     tAlphaGridN = gNVec[1]
                     sAlphaGridN = gNVec[2]
                     indL = 2 + (nt+1)*(0:(nt-2))
                     indD = 1 + (nt+1)*(0:(nt-1))
                     indU = (nt+1)*(1:(nt-1))
                     getTridiagonal = function(A) t(matrix(c(A[indU],0,A[indD],0,A[indL]),nrow(A),3))
                     iKtGrid = tf$stack(lapply(rLPar[[r]][[1]]$iWg, getTridiagonal))
                     iKsGrid = tf$stack(lapply(rLPar[[r]][[2]]$iWg, as.matrix))
                     iKsArray = tf$gather(iKsGrid, rep(0:(sAlphaGridN-1),tAlphaGridN), axis=ic(0))
                     iKtArray = tf$gather(iKtGrid, rep(0:(tAlphaGridN-1),each=sAlphaGridN), axis=ic(0))
                     ldLoopCond = function(it,logDet,L,iLm0,qF) tfm$less_equal(it, as.integer(nt))
                     ldLoopBody = function(it,logDet,L,iLm0,qF){
                        A = (iKsArray*iKtArray[,2,it-ic(1)][,NULL,NULL])[,NULL,,] + tfla$diag(lidl[,NULL]*tf$constant(numKronObsMat,tf$float64)[,it-ic(1)])
                        B = (iKsArray*iKtArray[,1,it-ic(2)][,NULL,NULL])[,NULL,,]
                        iLB = tfla$triangular_solve(L, B)
                        L = tfla$cholesky(A - tf$matmul(iLB,iLB,transpose_a=TRUE))
                        logDetV = tf$constant(2,tf$float64)*tf$reduce_sum(tfm$log(tfla$diag_part(L)), ic(-1))
                        iLm0 = tf$squeeze(tfla$triangular_solve(L, m0Mat[,it-ic(1),,NULL] - tf$matmul(iLB,iLm0[,,,NULL],transpose_a=TRUE)),ic(-1))
                        qF = qF + tf$reduce_sum(iLm0**2, ic(-1))
                        return(c(it+ic(1),logDet+logDetV,L,iLm0,qF))
                     }
                     it = 1
                     A = (iKsArray*iKtArray[,2,it,drop=FALSE])[,NULL,,] + tfla$diag(lidl[,NULL]*tf$constant(numKronObsMat,tf$float64)[,it])
                     LInit = tfla$cholesky(A)
                     iLm0Init = tf$squeeze(tfla$triangular_solve(LInit, m0Mat[,it,,NULL]))
                     qFInit = tf$reduce_sum(iLm0Init**2, ic(-1))
                     logDetInit = tf$constant(2,tf$float64)*tf$reduce_sum(tfm$log(tfla$diag_part(LInit)), ic(-1))
                     ldLoopInit = c(tf$constant(as.integer(2),tf$int32),logDetInit,LInit,iLm0Init,qFInit)
                     ldLoopRes = tf$while_loop(ldLoopCond, ldLoopBody, ldLoopInit)
                     logDetW = ldLoopRes[[2]]
                     qF2 = ldLoopRes[[5]]
                     return(c(qF2,logDetW))
                  }
                  qF2_logDetW_B2_tf_fun = tf_function(fun)
                  assign("updateAlpha_kronecker_B2_tf", qF2_logDetW_B2_tf_fun, envir=parent.frame())
               } else{
                  qF2_logDetW_B2_tf_fun = get("updateAlpha_kronecker_B2_tf", envir=parent.frame())
               }
               res = qF2_logDetW_B2_tf_fun(m0Mat,lidl)
               qF2 = res[[1]]
               logDetW = res[[2]]
               qF = y1iDy1 - qF2
               logDetKstMat = matrix(nt*rLPar[[r]][[2]]$detWg,gNVec[2],gNVec[1]) + matrix(nx*rLPar[[r]][[1]]$detWg,gNVec[2],gNVec[1],byrow=TRUE)
               logDet = logDetD + matrix(logDetKstMat,gNVec[2]*gNVec[1],1) + logDetW
               logLikeMat = tf$transpose(-0.5*qF -0.5*logDet)
            } else if(rL[[r]]$alphaMethod=="TF_direct_krylov"){
               nt = npElemVec[1]
               nx = npElemVec[2]
               if(np == nt*nx){
                  vMat = eta
               } else{
                  vMat = matrix(0,nt*nx,nf)
                  vMat[indKronObs,] = eta
               }
               vArray = tf$transpose(tf$constant(array(vMat,c(nx,nt,nf)),tf$float64), ic(2,1,0))
               if(!exists("updateAlpha_TF_direct_krylov", envir=parent.frame())){
                  fun = function(vArray){
                     rInd = r
                     tfla = tf$linalg
                     tfm = tf$math
                     obsMat = matrix(0,nx,nt)
                     obsMat[indKronObs] = 1
                     obsMat = tf$transpose(tf$constant(obsMat, tf$float64))
                     KtSt = tf$stack(rLPar[[r]][[1]]$Wg)
                     KsSt = tf$stack(rLPar[[r]][[2]]$Wg)
                     iKtSt = tfla$inv(KtSt)
                     iKsSt = tfla$inv(KsSt) # tf$stack(rLPar[[r]][[2]]$iWg)

                     x = tf$zeros(ic(gNVec[1],gNVec[2],nf,nt,nx),tf$float64)
                     r = tf$tile(vArray[NULL,NULL,,,], ic(gNVec[1],gNVec[2],1,1,1))
                     z = obsMat*tfla$matmul(iKtSt[,NULL,NULL,,], tf$matmul(r, iKsSt[NULL,,NULL,,]))
                     p = z
                     rTz = tf$reduce_sum(r*z, ic(-1,-2))
                     # for(cgIter in 1:rL[[rInd]]$cgIterN){
                     #    Ap = obsMat*tfla$matmul(KtSt[,NULL,NULL,,], tf$matmul(p, KsSt[NULL,,NULL,,]))
                     #    pTAp = tf$reduce_sum(p*Ap, ic(-1,-2))
                     #    a = tfm$divide_no_nan(rTz, pTAp)
                     #    x = x + a[,,,NULL,NULL]*p
                     #    rNew = r - a[,,,NULL,NULL]*Ap
                     #    zNew = obsMat*tfla$matmul(iKtSt[,NULL,NULL,,], tf$matmul(rNew, iKsSt[NULL,,NULL,,]))
                     #    rTrNew = tf$reduce_sum(rNew^2, ic(-1,-2))
                     #    rTzNew = tf$reduce_sum(rNew*zNew, ic(-1,-2))
                     #    EPS = 1e-12
                     #    b = tfm$multiply_no_nan((rTzNew/rTz), tf$cast(rTrNew>=EPS,tf$float64))
                     #    p = zNew + b[,,,NULL,NULL]*p
                     #    r=rNew; z=zNew; rTz=rTzNew
                     # }
                     cgLoopCond = function(cgIter,x,r,z,p,rTz) tfm$less_equal(cgIter, as.integer(rL[[rInd]]$cgIterN))
                     cgLoopBody = function(cgIter,x,r,z,p,rTz){
                        Ap = obsMat*tfla$matmul(KtSt[,NULL,NULL,,], tf$matmul(p, KsSt[NULL,,NULL,,]))
                        pTAp = tf$reduce_sum(p*Ap, ic(-1,-2))
                        a = tfm$divide_no_nan(rTz, pTAp)
                        x = x + a[,,,NULL,NULL]*p
                        rNew = r - a[,,,NULL,NULL]*Ap
                        zNew = obsMat*tfla$matmul(iKtSt[,NULL,NULL,,], tf$matmul(rNew, iKsSt[NULL,,NULL,,]))
                        rTrNew = tf$reduce_sum(rNew^2, ic(-1,-2))
                        rTzNew = tf$reduce_sum(rNew*zNew, ic(-1,-2))
                        EPS = 1e-12
                        b = tfm$multiply_no_nan((rTzNew/rTz), tf$cast(rTrNew>=EPS,tf$float64))
                        p = zNew + b[,,,NULL,NULL]*p
                        return(list(cgIter+ic(1),x,rNew,zNew,p,rTzNew))
                     }
                     cgLoopInit = list(tf$constant(ic(1),tf$int32),x,r,z,p,rTz)
                     cgLoopRes = tf$while_loop(cgLoopCond, cgLoopBody, cgLoopInit)
                     x = cgLoopRes[[2]]
                     qF = tf$reduce_sum(vArray*x, ic(-1,-2))
                     r = rInd
                     return(qF)
                  }
                  qF_B3_tf_fun = tf_function(fun)
                  assign("updateAlpha_TF_direct_krylov", qF_B3_tf_fun, envir=parent.frame())
               } else{
                  qF_B3_tf_fun = get("updateAlpha_TF_direct_krylov", envir=parent.frame())
               }
               qF = qF_B3_tf_fun(vArray)
               logLikeArray = -0.5*tf$transpose(qF,ic(2,0,1)) - 0.5*rLPar[[r]]$logDetK + log(rL[[r]]$alphaPrior$alphaProb)
               logLikeMat = tf$reshape(logLikeArray, ic(nf,gNVec[1]*gNVec[2]))
            }
            logLikeMat = logLikeMat - tf$reduce_logsumexp(logLikeMat, ic(-1), keepdims=TRUE)
            sampleInd = tf$random$categorical(logLikeMat, ic(1))$numpy() + 1
            Alpha[[r]][,2] = (sampleInd-1)%%gNVec[2]+1
            Alpha[[r]][,1] = (sampleInd-1)%/%gNVec[2]+1
         } else if(rL[[r]]$alphaMethod=="TF_direct_local_krylov"){
            rad = 1
            nt = npElemVec[1]
            nx = npElemVec[2]
            if(np == nt*nx){
               vMat = eta
            } else{
               vMat = matrix(0,nt*nx,nf)
               vMat[indKronObs,] = eta
            }
            vArray = tf$transpose(tf$constant(array(vMat,c(nx,nt,nf)),tf$float64), ic(2,1,0))
            gN = (2*rad+1)^2
            if(!exists("updateAlpha_TF_direct_local_krylov", envir=parent.frame())){
               fun = function(vArray, tInd, sInd){
                  rInd = r
                  tfla = tf$linalg
                  tfm = tf$math
                  obsMat = matrix(0,nx,nt)
                  obsMat[indKronObs] = 1
                  obsMat = tf$transpose(tf$constant(obsMat, tf$float64))
                  KtSt = tf$stack(rLPar[[r]][[1]]$Wg)
                  KsSt = tf$stack(rLPar[[r]][[2]]$Wg)
                  iKtSt = tfla$inv(KtSt)
                  iKsSt = tfla$inv(KsSt)

                  Kt = tf$gather(KtSt, tInd)
                  Ks = tf$gather(KsSt, sInd)
                  iKt = tf$gather(iKtSt, tInd)
                  iKs = tf$gather(iKsSt, sInd)
                  x = tf$zeros(ic(gN,nf,nt,nx),tf$float64)
                  r = tf$tile(vArray[NULL,,,], ic(gN,1,1,1))
                  z = obsMat*tfla$matmul(iKt, tf$matmul(r, iKs))
                  p = z
                  rTz = tf$reduce_sum(r*z, ic(-1,-2))
                  cgLoopCond = function(cgIter,x,r,z,p,rTz) tfm$less_equal(cgIter, as.integer(rL[[rInd]]$cgIterN))
                  cgLoopBody = function(cgIter,x,r,z,p,rTz){
                     Ap = obsMat*tfla$matmul(Kt, tf$matmul(p, Ks))
                     pTAp = tf$reduce_sum(p*Ap, ic(-1,-2))
                     a = tfm$divide_no_nan(rTz, pTAp)
                     x = x + a[,,NULL,NULL]*p
                     rNew = r - a[,,NULL,NULL]*Ap
                     zNew = obsMat*tfla$matmul(iKt, tf$matmul(rNew, iKs))
                     rTrNew = tf$reduce_sum(rNew^2, ic(-1,-2))
                     rTzNew = tf$reduce_sum(rNew*zNew, ic(-1,-2))
                     EPS = 1e-12
                     b = tfm$multiply_no_nan((rTzNew/rTz), tf$cast(rTrNew>=EPS,tf$float64))
                     p = zNew + b[,,NULL,NULL]*p
                     return(list(cgIter+ic(1),x,rNew,zNew,p,rTzNew))
                  }
                  cgLoopInit = list(tf$constant(ic(1),tf$int32),x,r,z,p,rTz)
                  cgLoopRes = tf$while_loop(cgLoopCond, cgLoopBody, cgLoopInit)
                  x = cgLoopRes[[2]]
                  qF = tf$reduce_sum(vArray*x, ic(-1,-2))
                  r = rInd
                  return(qF)
               }
               qF_C1_tf_fun = tf_function(fun)
               assign("updateAlpha_TF_direct_local_krylov", qF_C1_tf_fun, envir=parent.frame())
            } else{
               qF_C1_tf_fun = get("updateAlpha_TF_direct_local_krylov", envir=parent.frame())
            }
            tInd = matrix(rep(-rad:rad,each=2*rad+1),gN,nf) + matrix(Alpha[[r]][,1],gN,nf,byrow=TRUE) - 1
            sInd = matrix(rep(-rad:rad,2*rad+1),gN,nf) + matrix(Alpha[[r]][,2],gN,nf,byrow=TRUE) - 1
            remFlag = tf$constant(((tInd<0) | (tInd>=gNVec[1]) | (sInd<0) | (sInd>=gNVec[2])), tf$float64)
            tInd = (1-remFlag)*tInd
            sInd = (1-remFlag)*sInd
            qF = qF_C1_tf_fun(vArray, tf$cast(tInd,tf$int32), tf$cast(sInd,tf$int32))
            indMat = tf$cast(tf$stack(list(tInd,sInd),ic(-1)),tf$int32)
            logDetK = tf$gather_nd(tf$constant(rLPar[[r]]$logDetK,tf$float64), indMat)
            logPriorProb = tf$math$log(tf$gather_nd(tf$constant(rL[[r]]$alphaPrior$alphaProb,tf$float64), indMat))
            logLikeMat = -0.5*qF - 0.5*logDetK + logPriorProb

            logLikeMat = logLikeMat*(1-remFlag) - tf$math$multiply_no_nan(tf$constant(Inf,tf$float64),remFlag)
            logLikeMat = logLikeMat - tf$reduce_logsumexp(logLikeMat, ic(0), keepdims=TRUE)
            sampleInd = tf$squeeze(tf$random$categorical(tf$transpose(logLikeMat), ic(1), tf$int32), ic(-1))
            sampleIndMat = tf$stack(list(sampleInd,tf$constant(0:(nf-1),tf$int32)), ic(-1))
            Alpha[[r]][,2] = tf$gather_nd(sInd,sampleIndMat)$numpy() + 1
            Alpha[[r]][,1] = tf$gather_nd(tInd,sampleIndMat)$numpy() + 1
         }
      }
   }
   return(Alpha)
}


