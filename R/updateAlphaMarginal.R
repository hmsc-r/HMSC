#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#' @importFrom matrixStats colMaxs logSumExp
#' @importFrom pracma ndims
#' @importFrom plyr mdply
#' @importFrom tensorflow tf
#'
updateAlphaMarginal = function(iter,Z,Beta,iSigma,Eta,EtaFull,Alpha,Lambda, rLPar, Y,X,Pi,dfPi,rL){
   tfla = tf$linalg
   tfm = tf$math
   ic = function(...){
      return(as.integer(c(...)))
   }
   nr = length(rL)
   ny = nrow(Z)
   ns = ncol(Z)
   iD = (!is.na(Y)) * matrix(iSigma,ny,ns,byrow=TRUE)
   logDetD = sum((!is.na(Y)) * matrix((-log(iSigma)),ny,ns,byrow=TRUE))

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
      eta = Eta[[r]]
      lambda = Lambda[[r]]
      np = nrow(eta)
      nf = ncol(eta)
      if(inherits(rL[[r]],"HmscRandomLevel",TRUE)==1){
         # do nothing
      } else if(inherits(rL[[r]],"HmscSpatialRandomLevel",TRUE)==1){
         #TODO:GT probably worth to implement this later
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
         nt = npElemVec[1]; nx = npElemVec[2]
         logDetKt = tf$constant(rLPar[[r]][[1]]$detWg, tf$float64)
         logDetKs = tf$constant(rLPar[[r]][[2]]$detWg, tf$float64)
         logDetKstMat = nx*logDetKt[,NULL] + nt*logDetKs[NULL,]

         for(h in ((iter-1)%%nf + 1)){ #1:nf
            S = Z - (LFix + Reduce("+",LRan)) + tcrossprod(eta[Pi[,r],h], lambda[h,])
            LambdaiDS = as.vector((iD*S) %*% lambda[h,])
            m0Mat = tf$reshape(tf$scatter_nd(matrix(ic(indKronObs[Pi[,r]]-1)), LambdaiDS, ic(nt*nx)*tf$ones(ic(1),tf$int32)), ic(nt,nx))
            m0Mat = tf$cast(m0Mat, tf$float64)
            tmp = tf$einsum("j,ij,j->i",tf$constant(lambda[h,]),iD,tf$constant(lambda[h,]))
            LamiDLam = tf$reshape(tf$scatter_nd(matrix(ic(indKronObs[Pi[,r]]-1)), tmp, ic(nt*nx)*tf$ones(ic(1),tf$int32)), ic(nt,nx))
            LamiDLam = tf$cast(LamiDLam, tf$float64)
            SiDS = tf$reduce_sum(iD * S**2)

            if(rL[[r]]$alphaMarginalMethod=="TF_obs"){ #TODO to be checked thoughtfully
               if(!exists("updateAlpha_kronecker_TF_obs", envir=parent.frame())){
                  L_var = iLm0_var = qF_var = logDet_var = NULL
                  fun = function(m0Mat,LamiDLam){
                     iKtArray = tf$gather(rLPar[[r]][[1]]$iWStack, rep(0:(gNVec[1]-1),each=gNVec[2]), axis=ic(0))
                     iKsArray = tf$gather(rLPar[[r]][[2]]$iWStack, rep(0:(gNVec[2]-1),gNVec[1]), axis=ic(0))
                     it = tf$constant(ic(0),tf$int32)
                     A = iKsArray*iKtArray[,2,it,NULL,NULL] + tfla$diag(LamiDLam[it,])
                     LInit = tfla$cholesky(A)
                     iLm0Init = tf$squeeze(tfla$triangular_solve(LInit, m0Mat[it,,NULL]), ic(-1))
                     qFInit = tf$reduce_sum(iLm0Init^2, ic(-1))
                     logDetInit = tf$constant(2,tf$float64)*tf$reduce_sum(tfm$log(tfla$diag_part(LInit)), ic(-1))
                     if(is.null(L_var)) L_var<<-tf$Variable(LInit) else L_var$assign(LInit)
                     if(is.null(iLm0_var)) iLm0_var<<-tf$Variable(iLm0Init) else iLm0_var$assign(iLm0Init)
                     if(is.null(qF_var)) qF_var<<-tf$Variable(qFInit) else qF_var$assign(qFInit)
                     if(is.null(logDet_var)) logDet_var<<-tf$Variable(logDetInit) else logDet_var$assign(logDetInit)
                     for(it in tf$range(ic(1),ic(nt))){
                        L=L_var; iLm0=iLm0_var
                        A = iKsArray*iKtArray[,2,it,NULL,NULL] + tfla$diag(LamiDLam[it,])
                        B = iKsArray*iKtArray[,1,it-ic(1),NULL,NULL]
                        iLB = tfla$triangular_solve(L, B)
                        L = tfla$cholesky(A - tf$matmul(iLB,iLB,transpose_a=TRUE))
                        logDetV = tf$constant(2,tf$float64)*tf$reduce_sum(tfm$log(tfla$diag_part(L)), ic(-1))
                        iLm0 = tf$squeeze(tfla$triangular_solve(L, m0Mat[it,,NULL] - tf$matmul(iLB,iLm0[,,NULL],transpose_a=TRUE)),ic(-1))
                        L_var$assign(L); iLm0_var$assign(iLm0)
                        qF_var$assign_add(tf$reduce_sum(iLm0^2, ic(-1)))
                        logDet_var$assign_add(logDetV)
                     }
                     return(c(qF_var,logDet_var))
                  }
                  qF2_logDetW_obs_tf_fun = tf_function(fun) #tf_function
                  assign("updateAlpha_kronecker_TF_obs", qF2_logDetW_obs_tf_fun, envir=parent.frame())
               } else{
                  qF2_logDetW_obs_tf_fun = get("updateAlpha_kronecker_TF_obs", envir=parent.frame())
               }
               res = qF2_logDetW_obs_tf_fun(m0Mat,LamiDLam)
               qF2 = tf$reshape(res[[1]],ic(gNVec))
               logDetW = tf$reshape(res[[2]],ic(gNVec))
               qF = SiDS - qF2
               logDet = logDetD + logDetKstMat + logDetW
               logLikeMat = -0.5*qF -0.5*logDet
               logPostVec = tf$reshape(logLikeMat + log(rL[[r]]$alphaPrior$alphaProb), ic(prod(gNVec)))
               logPostVec = logPostVec - tf$reduce_logsumexp(logPostVec)
               sampleInd = tf$squeeze(tf$random$categorical(logPostVec[NULL,], ic(1)), ic(0))$numpy() + 1
               Alpha[[r]][h,2] = (sampleInd-1)%%gNVec[2]+1
               Alpha[[r]][h,1] = (sampleInd-1)%/%gNVec[2]+1
            }
         }
      }
   }
   return(list(Alpha=Alpha,Eta=Eta))
}



# else if(rL[[r]]$alphaMarginalMethod=="TF_obs_local"){ #TODO:GT unfinished part
#    rad = 1
#    nt = npElemVec[1]
#    nx = npElemVec[2]
#    gN = (2*rad+1)^2
#    P = sparseMatrix(i=1:ny,j=Pi[,r])
#    numObsVec = Matrix::colSums(P)
#    Y1 = S - LRan[[r]] + tf$einsum("ih,hj->hij", eta[Pi[,r],], lambda)
#    iDY1 = Y1 * matrix(iSigma,ny,ns,byrow=TRUE)
#    LambdaiDY1 = tf$einsum("hj,hij->hi",lambda,iDY1)
#    LambdaiDY1 = tf$transpose(tf$scatter_nd(tf$constant(matrix(ic(Pi[,r]-1),ny,1)), tf$transpose(LambdaiDY1), tf$constant(ic(np,nf))))
#    if(np == nt*nx){
#       m0Mat = LambdaiDY1
#       numKronObsMat = matrix(numObsVec,nx,nt)
#    } else{
#       m0Mat = matrix(0,nf,nt*nx)
#       m0Mat[,indKronObs] = LambdaiDY1$numpy()
#       numKronObsMat = matrix(0,nx,nt)
#       numKronObsMat[indKronObs] = numObsVec
#       m0Mat = tf$constant(m0Mat, tf$float64)
#    }
#    m0Mat = tf$reshape(m0Mat, ic(nf,nt,nx))
#    y1iDy1 = tf$reduce_sum(Y1**2 * tf$constant(iSigma,tf$float64), ic(1,2))
#    logDetD = -ny*sum(log(iSigma))
#    lidl = tf$einsum("hj,j,hj->h",lambda,tf$constant(iSigma,tf$float64),lambda)
#
#    logPriorProb = tf$math$log(tf$gather_nd(tf$constant(rL[[r]]$alphaPrior$alphaProb,tf$float64), indMat))
#    logLikeMat = -0.5*qF - 0.5*logDetK + logPriorProb
#    logLikeMat = logLikeMat*(1-remFlag) - tf$math$multiply_no_nan(tf$constant(Inf,tf$float64),remFlag)
#    logLikeMat = logLikeMat - tf$reduce_logsumexp(logLikeMat, ic(0), keepdims=TRUE)
#    sampleInd = tf$squeeze(tf$random$categorical(tf$transpose(logLikeMat), ic(1), tf$int32), ic(-1))
#    sampleIndMat = tf$stack(list(sampleInd,tf$constant(0:(nf-1),tf$int32)), ic(-1))
#    Alpha[[r]][,2] = tf$gather_nd(sInd,sampleIndMat)$numpy() + 1
#    Alpha[[r]][,1] = tf$gather_nd(tInd,sampleIndMat)$numpy() + 1
# }

