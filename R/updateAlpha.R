#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#' @importFrom matrixStats colMaxs logSumExp
#' @importFrom pracma ndims
#' @importFrom plyr mdply
#' @importFrom tensorflow tf
#'
updateAlpha = function(Z,Beta,iSigma,Eta,Lambda, rLPar, X,Pi,dfPi,rL){
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
   Alpha = vector("list", nr)
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
         Alpha[[r]] = matrix(NA, nf, krN)
         dfPiElemLevelList = vector("list", krN)
         npElemVec = rep(NA, krN)
         for(l in seq_len(krN)){
            dfPiElem = as.factor(unlist(lapply(strsplit(as.character(m$dfPi[,r]), "="), function(a) a[l])))
            dfPiElemLevelList[[l]] = levels(dfPiElem)
            npElemVec[l] = length(dfPiElemLevelList[[l]])
         }
         dfTmp = expand.grid(rev(dfPiElemLevelList))[rev(1:krN)]
         # allUnits = as.factor(mdply(dfTmp, paste, sep=rL[[r]]$sepStr)[,length(dfPiElemLevelList)+1])
         allUnits = factor(do.call(function(...) paste(..., sep=rL[[r]]$sepStr),dfTmp))
         indKronObs = as.numeric(factor(m$dfPi[,r], levels=levels(allUnits)))
         gNVec = sapply(rL[[r]]$alphaPrior$alphaGridList, length)
         if(rL[[r]]$alphaMethod == "R"){
            if(length(allUnits) == np[r]){
               etaArray = array(eta,c(rev(npElemVec),nf))
               iWgEtaList = vector("list", krN)
               for(l in seq_len(krN)){
                  if(rL[[r]]$rLList[[l]]$sDim == 0){
                     iWgEtaList[[l]] = array(etaArray, c(1,dim(etaArray)))
                  } else{
                     gN = gNVec[l]
                     iWg = rLPar[[r]][[l]]$iWg
                     # iWgEta = array(NA, c(gN,dim(etaArray)))
                     iWgEtaElemList = vector("list",gN)
                     dimInd = c(krN-l+1,c(1:ndims(etaArray))[-(krN-l+1)])
                     dimIndOrd = order(dimInd)
                     etaMatFlat = matrix(aperm(etaArray, dimInd), dim(etaArray)[krN-l+1])
                     for(ag in seq_len(gN)){
                        iWgEtaMatFlat = iWg[[ag]]%*%etaMatFlat
                        # agInd = ag + (1:length(etaArray) - 1)*gN
                        # iWgEta[agInd] = aperm(array(iWgEtaMatFlat, dim(etaArray)[dimInd]), dimIndOrd)
                        iWgEtaElemList[[ag]] = array(iWgEtaMatFlat, dim(etaArray)[dimInd])
                     }
                     # abind(rLPar[[r]][[l]]$iWg, along=0)
                     iWgEta = aperm(abind(iWgEtaElemList, rev.along=0), c(krN+2,dimIndOrd))
                     iWgEtaList[[l]] = iWgEta
                  }
               }
               #TODO:GT potential for >2 kronecker elements
               gN1 = gNVec[1]
               gN2 = gNVec[2]
               # tmp1 = aperm(array(rep(iWgEtaList[[1]], each=gN2), c(gN2,gN1,dim(etaArray))), c(2,1,2+1:ndims(etaArray)))
               # tmp2 = array(rep(iWgEtaList[[2]], each=gN1), c(gN1,gN2,dim(etaArray)))
               # EtaiWEta = tmp1 * tmp2
               # qFArray = apply(EtaiWEta, c(1,2,2+ndims(etaArray)), sum)
               qFArray = tf$linalg$einsum("axtf,bxtf->abf",iWgEtaList[[1]],iWgEtaList[[2]])$numpy()
               tmp1 = matrix(rLPar[[r]][[1]]$detWg, gN1,gN2)
               tmp2 = matrix(rLPar[[r]][[2]]$detWg, gN1,gN2, byrow=TRUE)
               logDetWArray = array(npElemVec[2]*tmp1+npElemVec[1]*tmp2, c(gN1,gN2,nf))
               logPriorArray = array(log(rL[[r]]$alphaPrior$alphaProb), c(gN1,gN2,nf))
               logProbArray = logPriorArray - 0.5*logDetWArray - 0.5*qFArray #dimension-based const coefficients are omitted in this expr
               logProbArray1 = logProbArray
            } else{
               # stop("updateAlpha: not implemented yet")
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
         } else if(rL[[r]]$alphaMethod == "TF"){
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
            logLikeMat = logLikeMat - tf$reduce_logsumexp(logLikeMat, ic(-1), keepdims=TRUE)
            sampleInd = tf$random$categorical(logLikeMat, ic(1))$numpy() + 1
            Alpha[[r]][,2] = (sampleInd-1)%%gNVec[2]+1
            Alpha[[r]][,1] = (sampleInd-1)%/%gNVec[2]+1
         }
      }
   }
   return(Alpha)
}


