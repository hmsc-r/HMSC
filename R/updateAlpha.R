#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#' @importFrom matrixStats colMaxs logSumExp
#' @importFrom pracma ndims
#' @importFrom plyr mdply
#' @importFrom tensorflow tf
#'
updateAlpha = function(Eta,EtaFull,Alpha, rLPar, dfPi,rL){
   tfla = tf$linalg
   tfm = tf$math
   ic = function(...){
      return(as.integer(c(...)))
   }
   nr = length(rL)

   for(r in seq_len(nr)){
      eta = Eta[[r]]
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
         dfPiElemLevelList = vector("list", krN)
         npElemVec = rep(NA, krN)
         for(l in seq_len(krN)){
            dfPiElem = as.factor(unlist(lapply(strsplit(as.character(dfPi[,r]), rL[[r]]$sepStr), function(a) a[l])))
            dfPiElemLevelList[[l]] = levels(dfPiElem)
            npElemVec[l] = length(dfPiElemLevelList[[l]])
         }
         dfTmp = expand.grid(rev(dfPiElemLevelList))[,rev(1:krN)]
         allUnits = factor(do.call(function(...) paste(..., sep=rL[[r]]$sepStr),dfTmp))
         indKronObs = as.numeric(factor(levels(dfPi[,r]), levels=levels(allUnits)))
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
               ind = as.numeric(factor(dfPi[,r], levels=levels(allUnits)))
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
         } else if(rL[[r]]$alphaMethod=="TF_full" || rL[[r]]$alphaMethod=="TF_direct_krylov"){
            if(rL[[r]]$alphaMethod=="TF_full"){
               if(!exists("updateAlpha_TF_full", envir=parent.frame())){
                  fun = function(etaArray){
                     iWgEtaList = vector("list", krN)
                     for(l in seq_len(krN)){
                        iWg = rLPar[[r]][[l]]$iWStack
                        if(l==1){
                           iWgEta = tf$linalg$tridiagonal_matmul(tf$tile(iWg[,NULL,,],ic(1,npElemVec[2],1,1)),tf$tile(etaArray[NULL,,,],ic(gNVec[1],1,1,1)))
                        }
                        if(l==2){
                           iWgEta = tf$einsum("gyx,xtf->gytf",iWg,etaArray)
                        }
                        iWgEtaList[[l]] = iWgEta
                     }
                     qFArray = tf$linalg$einsum("axtf,bxtf->fab",iWgEtaList[[1]],iWgEtaList[[2]])
                     tmp1 = matrix(rLPar[[r]][[1]]$detWg, gNVec[1],gNVec[2])
                     tmp2 = matrix(rLPar[[r]][[2]]$detWg, gNVec[1],gNVec[2], byrow=TRUE)
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
            } else if(rL[[r]]$alphaMethod=="TF_direct_krylov"){
               nt = npElemVec[1]; nx = npElemVec[2]
               if(np == nt*nx){
                  vMat = eta
               } else{
                  vMat = matrix(0,nt*nx,nf)
                  vMat[indKronObs,] = eta
               }
               vArray = tf$transpose(tf$constant(array(vMat,c(nx,nt,nf)),tf$float64), ic(2,1,0))
               if(!exists("updateAlpha_TF_direct_krylov", envir=parent.frame())){
                  x_var = r_var = z_var = p_var = rTz_var = NULL
                  fun = function(vArray){
                     rInd = r
                     obsMat = matrix(0,nx,nt)
                     obsMat[indKronObs] = 1
                     obsMat = tf$transpose(tf$constant(obsMat, tf$float64))
                     KtSt = rLPar[[rInd]][[1]]$WStack
                     KsSt = rLPar[[rInd]][[2]]$WStack
                     iKtSt = rLPar[[rInd]][[1]]$iWStack
                     iKsSt = rLPar[[rInd]][[2]]$iWStack

                     x = tf$zeros(ic(gNVec[1],gNVec[2],nf,nt,nx),tf$float64)
                     r = tf$tile(vArray[NULL,NULL,,,], ic(gNVec[1],gNVec[2],1,1,1))
                     z = obsMat*tfla$tridiagonal_matmul(tf$tile(iKtSt[,NULL,NULL,,],ic(1,gNVec[2],nf,1,1)),
                                                        tf$matmul(r, iKsSt[NULL,,NULL,,]))
                     p = z
                     rTz = tf$reduce_sum(r*z, ic(-1,-2))
                     if(is.null(x_var)) x_var<<-tf$Variable(x) else x_var$assign(x)
                     if(is.null(r_var)) r_var<<-tf$Variable(r) else r_var$assign(r)
                     if(is.null(z_var)) z_var<<-tf$Variable(z) else z_var$assign(z)
                     if(is.null(p_var)) p_var<<-tf$Variable(p) else p_var$assign(p)
                     if(is.null(rTz_var)) rTz_var<<-tf$Variable(rTz) else rTz_var$assign(rTz)
                     for(cgIter in tf$range(ic(rL[[rInd]]$cgIterN))){
                        x=x_var; r=r_var; z=z_var; p=p_var; rTz=rTz_var
                        Ap = obsMat*tfla$matmul(KtSt[,NULL,NULL,,], tf$matmul(p, KsSt[NULL,,NULL,,]))
                        pTAp = tf$reduce_sum(p*Ap, ic(-1,-2))
                        a = tfm$divide_no_nan(rTz, pTAp)
                        x = x + a[,,,NULL,NULL]*p
                        rNew = r - a[,,,NULL,NULL]*Ap
                        zNew = obsMat*tfla$tridiagonal_matmul(tf$tile(iKtSt[,NULL,NULL,,],ic(1,gNVec[2],nf,1,1)),
                                                              tf$matmul(rNew, iKsSt[NULL,,NULL,,]))
                        rTrNew = tf$reduce_sum(rNew^2, ic(-1,-2))
                        rTzNew = tf$reduce_sum(rNew*zNew, ic(-1,-2))
                        EPS = 1e-12
                        b = tfm$multiply_no_nan((rTzNew/rTz), tf$cast(rTrNew>=EPS,tf$float64))
                        p = zNew + b[,,,NULL,NULL]*p
                        r=rNew; z=zNew; rTz=rTzNew
                        x_var$assign(x); r_var$assign(r); z_var$assign(z); p_var$assign(p); rTz_var$assign(rTz)
                        # tf$print(tf$reduce_sum(vArray*x_var, ic(-1,-2))[ic(2),ic(2),])
                     }
                     qF = tf$reduce_sum(vArray*x_var, ic(-1,-2))
                     r = rInd
                     return(qF)
                  }
                  qF_B3_tf_fun = tf_function(fun)
                  assign("updateAlpha_TF_direct_krylov", qF_B3_tf_fun, envir=parent.frame())
               } else{
                  qF_B3_tf_fun = get("updateAlpha_TF_direct_krylov", envir=parent.frame())
               }
               qF = qF_B3_tf_fun(vArray)
               logLikeArray = -0.5*tf$transpose(qF,ic(2,0,1)) - 0.5*rLPar[[r]]$logDetK
               logPostArray = logLikeArray + log(rL[[r]]$alphaPrior$alphaProb)
               logPostMat = tf$reshape(logPostArray, ic(nf,gNVec[1]*gNVec[2]))
            }
            logPostMat = logPostMat - tf$reduce_logsumexp(logPostMat, ic(-1), keepdims=TRUE)
            sampleInd = tf$random$categorical(logPostMat, ic(1))$numpy() + 1
            Alpha[[r]][,2] = (sampleInd-1)%%gNVec[2]+1
            Alpha[[r]][,1] = (sampleInd-1)%/%gNVec[2]+1
         }
      }
   }
   return(Alpha)
}


# cgLoopCond = function(cgIter,x,r,z,p,rTz) tfm$less_equal(cgIter, as.integer(rL[[rInd]]$cgIterN))
# cgLoopBody = function(cgIter,x,r,z,p,rTz){
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
#    return(list(cgIter+ic(1),x,rNew,zNew,p,rTzNew))
# }
# cgLoopInit = list(tf$constant(ic(1),tf$int32),x,r,z,p,rTz)
# cgLoopRes = tf$while_loop(cgLoopCond, cgLoopBody, cgLoopInit)
# x = cgLoopRes[[2]]



# else{
#    stop("Local samplers are currently incorrect. Logic of local samplers must be reconsidered and fixed.")
#    rad = 1
#    nt = npElemVec[1]; nx = npElemVec[2]
#    gN = (2*rad+1)^2
#    if(rL[[r]]$alphaMethod=="TF_direct_local_krylov"){
#       if(np == nt*nx){
#          vMat = eta
#       } else{
#          vMat = matrix(0,nt*nx,nf)
#          vMat[indKronObs,] = eta
#       }
#       if(!exists("updateAlpha_TF_direct_local_krylov", envir=parent.frame())){
#          fun = function(vArray, tInd, sInd){
#             rInd = r
#             tfla = tf$linalg
#             tfm = tf$math
#             obsMat = matrix(0,nx,nt)
#             obsMat[indKronObs] = 1
#             obsMat = tf$transpose(tf$constant(obsMat, tf$float64))
#             KtSt = rLPar[[r]][[1]]$WStack
#             KsSt = rLPar[[r]][[2]]$WStack
#             iKtSt = rLPar[[r]][[1]]$iWStack
#             iKsSt = rLPar[[r]][[2]]$iWStack
#
#             Kt = tf$gather(KtSt, tInd)
#             Ks = tf$gather(KsSt, sInd)
#             iKt = tf$gather(iKtSt, tInd)
#             iKs = tf$gather(iKsSt, sInd)
#             x = tf$zeros(ic(gN,nf,nt,nx),tf$float64)
#             r = tf$tile(vArray[NULL,,,], ic(gN,1,1,1))
#             z = obsMat*tfla$tridiagonal_matmul(iKt, tf$matmul(r, iKs))
#             p = z
#             rTz = tf$reduce_sum(r*z, ic(-1,-2))
#             cgLoopCond = function(cgIter,x,r,z,p,rTz) tfm$less_equal(cgIter, as.integer(rL[[rInd]]$cgIterN))
#             cgLoopBody = function(cgIter,x,r,z,p,rTz){
#                Ap = obsMat*tfla$matmul(Kt, tf$matmul(p, Ks))
#                pTAp = tf$reduce_sum(p*Ap, ic(-1,-2))
#                a = tfm$divide_no_nan(rTz, pTAp)
#                x = x + a[,,NULL,NULL]*p
#                rNew = r - a[,,NULL,NULL]*Ap
#                zNew = obsMat*tfla$tridiagonal_matmul(iKt, tf$matmul(rNew, iKs))
#                rTrNew = tf$reduce_sum(rNew^2, ic(-1,-2))
#                rTzNew = tf$reduce_sum(rNew*zNew, ic(-1,-2))
#                EPS = 1e-12
#                b = tfm$multiply_no_nan((rTzNew/rTz), tf$cast(rTrNew>=EPS,tf$float64))
#                p = zNew + b[,,NULL,NULL]*p
#                return(list(cgIter+ic(1),x,rNew,zNew,p,rTzNew))
#             }
#             cgLoopInit = list(tf$constant(ic(1),tf$int32),x,r,z,p,rTz)
#             cgLoopRes = tf$while_loop(cgLoopCond, cgLoopBody, cgLoopInit)
#             x = cgLoopRes[[2]]
#             qF = tf$reduce_sum(vArray*x, ic(-1,-2))
#             r = rInd
#             return(qF)
#          }
#          qF_C1_tf_fun = tf_function(fun)
#          assign("updateAlpha_TF_direct_local_krylov", qF_C1_tf_fun, envir=parent.frame())
#       } else{
#          qF_C1_tf_fun = get("updateAlpha_TF_direct_local_krylov", envir=parent.frame())
#       }
#       tInd = matrix(rep(-rad:rad,each=2*rad+1),gN,nf) + matrix(Alpha[[r]][,1],gN,nf,byrow=TRUE) - 1
#       sInd = matrix(rep(-rad:rad,2*rad+1),gN,nf) + matrix(Alpha[[r]][,2],gN,nf,byrow=TRUE) - 1
#       remFlag = tf$constant(((tInd<0) | (tInd>=gNVec[1]) | (sInd<0) | (sInd>=gNVec[2])), tf$float64)
#       tInd = (1-remFlag)*tInd
#       sInd = (1-remFlag)*sInd
#       vArray = tf$transpose(tf$constant(array(vMat,c(nx,nt,nf)),tf$float64), ic(2,1,0))
#       qF = qF_C1_tf_fun(vArray, tf$cast(tInd,tf$int32), tf$cast(sInd,tf$int32))
#       indMat = tf$cast(tf$stack(list(tInd,sInd),ic(-1)),tf$int32)
#       logDetK = tf$gather_nd(tf$constant(rLPar[[r]]$logDetK,tf$float64), indMat)
#    }
#    logPriorProb = tf$math$log(tf$gather_nd(tf$constant(rL[[r]]$alphaPrior$alphaProb,tf$float64), indMat))
#    logLikeMat = -0.5*qF - 0.5*logDetK + logPriorProb
#
#    logLikeMat = logLikeMat*(1-remFlag) - tf$math$multiply_no_nan(tf$constant(Inf,tf$float64),remFlag)
#    logLikeMat = logLikeMat - tf$reduce_logsumexp(logLikeMat, ic(0), keepdims=TRUE)
#    sampleInd = tf$squeeze(tf$random$categorical(tf$transpose(logLikeMat), ic(1), tf$int32), ic(-1))
#    sampleIndMat = tf$stack(list(sampleInd,tf$constant(0:(nf-1),tf$int32)), ic(-1))
#    # print(vMat)
#    # print(tf$round(logLikeMat,2))
#    # print(sampleIndMat)
#    Alpha[[r]][,2] = tf$gather_nd(sInd,sampleIndMat)$numpy() + 1
#    Alpha[[r]][,1] = tf$gather_nd(tInd,sampleIndMat)$numpy() + 1
# }

