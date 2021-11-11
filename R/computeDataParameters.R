#' @title computeDataParameters
#'
#' @description Computes initial values before the sampling starts
#'
#' @param hM a fitted \code{Hmsc} model object
#'
#' @return a list including pre-computed matrix inverses and determinants (for phylogenetic and spatial random effects) needed in MCMC sampling
#'
#' @importFrom stats dist
#' @importFrom methods is
#' @importFrom sp spDists
#' @importFrom FNN get.knn
#' @importFrom Matrix .sparseDiagonal t solve band
#' @importFrom tensorflow tf
#'
#' @export


computeDataParameters = function(hM){
   tfla = tf$linalg
   tfm = tf$math
   ic = function(...){
      return(as.integer(c(...)))
   }
   parList = list()

   if(!is.null(hM$C)){
      Qg = array(NA, c(hM$ns,hM$ns,nrow(hM$rhopw)))
      iQg = array(NA, c(hM$ns,hM$ns,nrow(hM$rhopw)))
      RQg = array(NA, c(hM$ns,hM$ns,nrow(hM$rhopw)))
      detQg = rep(NA, nrow(hM$rhopw))
      if(any(hM$rhopw[,1] < 0))
         iC = chol2inv(chol(hM$C))
      for(rg in 1:nrow(hM$rhopw)){
         rho = hM$rhopw[rg,1]
         if(rho >= 0){
            rhoC = rho*hM$C;
         } else{
            rhoC = (-rho)*iC;
         }
         Q = rhoC + (1-abs(rho))*diag(hM$ns)
         Qg[,,rg] = Q
         RQ = chol(Q);
         iQg[,,rg] = chol2inv(RQ)
         RQg[,,rg] = RQ
         detQg[rg] = 2*sum(log(diag(RQ)))
      }
   } else{
      Qg = array(diag(hM$ns), c(hM$ns,hM$ns,1))
      iQg = array(diag(hM$ns), c(hM$ns,hM$ns,1))
      detQg = 0
      RQg = array(diag(hM$ns), c(hM$ns,hM$ns,1))
   }

   rLPar = vector("list", hM$nr)
   for(r in seq_len(hM$nr)){
      # print(class(hM$rL[[r]]))
      if(inherits(hM$rL[[r]],"HmscSpatialRandomLevel",TRUE)==1){
         if(hM$rL[[r]]$sDim > 0){
            alphapw = hM$rL[[r]]$alphapw
            np = hM$np[r]
            alphaN = nrow(alphapw)
            switch(hM$rL[[r]]$spatialMethod,
                   "Full" = {
                      if(is.null(hM$rL[[r]]$distMat)){
                         s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                         if (is(s, "Spatial"))
                            distance <- spDists(s)
                         else
                            distance = as.matrix(dist(s))
                      } else{
                         distance = hM$rL[[r]]$distMat[levels(hM$dfPi[,r]),levels(hM$dfPi[,r])]
                      }
                      Wg = vector("list", alphaN)
                      iWg = vector("list", alphaN)
                      RiWg = vector("list", alphaN)
                      detWg = rep(NA, alphaN)
                      for(ag in 1:alphaN){
                         alpha = alphapw[ag,1]
                         if(alpha==0){
                            W = diag(np)
                         } else{
                            W = exp(-distance/alpha) #TODO if the spatial dimensionality is 1, then the sparsity in iWg shall be accounted for
                         }
                         RW = chol(W)
                         iW = chol2inv(RW)
                         if(hM$rL[[r]]$sDim > 0){
                            iW = Matrix(band(iW,-1,1), sparse=TRUE)
                         }
                         Wg[[ag]] = W
                         iWg[[ag]] = iW
                         RiWg[[ag]] = chol(iW)
                         detWg[ag] = 2*sum(log(diag(RW)))
                      }
                      rLPar[[r]] = list(Wg=Wg, iWg=iWg, RiWg=RiWg, detWg=detWg)
                   },
                   "NNGP" = {
                      if(is.null(hM$rL[[r]]$nNeighbours)){
                         hM$rL[[r]]$nNeighbours = 10
                      }
                      if(!is.null(hM$rL[[r]]$distMat)){
                         dnam <- levels(hM$dfPi[,r])
                         distMat <- hM$rL[[r]]$distMat[dnam, dnam]
                      }
                      ## SpatialPoints are non-Euclidean, and we need
                      ## distances to get the nearest neighbours
                      if(is(hM$rL[[r]]$s, "Spatial"))
                         distMat <- spDists(hM$rL[[r]]$s[levels(hM$dfPi[,r]),])
                      iWg = vector("list", alphaN)
                      RiWg = vector("list", alphaN)
                      detWg = rep(NA, alphaN)
                      ## get nNeighbours nearest neighbours. If we
                      ## already have distances, we use them, but
                      ## otherwise use a fast method finding nearest
                      ## (Euclidean) neighbours and find their distances
                      if (exists("distMat", inherits = FALSE)) {
                         k <- seq_len(hM$rL[[r]]$nNeighbours)
                         diag(distMat) <- Inf
                         indNN <- t(apply(distMat, 2, function(x)
                            sort(x, index.return = TRUE)$ix[k]))
                         diag(distMat) <- 0
                      } else { # fast with coordinate data and FNN package
                         s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                         indNN = get.knn(s,k=hM$rL[[r]]$nNeighbours)[[1]]
                      }
                      indNN = t(apply(indNN,1,sort,decreasing=FALSE))
                      indices = list()
                      distList = list()
                      for(i in 2:np){
                         ind = indNN[i,]
                         ind = ind[ind<i]
                         if(!is.na(ind[1])){
                            indices[[i]] = rbind(i*rep(1,length(ind)),ind)
                            if (exists("distMat", inherits = FALSE))
                               distList[[i]] <- distMat[c(ind,i), c(ind,i)]
                            else
                               distList[[i]] = as.matrix(dist(s[c(ind,i),]))
                         }
                      }
                      for (ag in 1:alphaN){
                         alpha = alphapw[ag,1]
                         if(alpha==0){
                            iW = .sparseDiagonal(np)
                            RiW = .sparseDiagonal(np)
                            detW = 0
                         } else{
                            D = rep(0,np)
                            D[1] = 1
                            values = list()
                            for (i in 2:np){
                               if(!is.null(indices[[i]][2])){
                                  Kp = exp(-distList[[i]]/alpha)
                                  values[[i]] = solve(Kp[1:(nrow(Kp)-1),1:(ncol(Kp)-1)],Kp[1:(nrow(Kp)-1),ncol(Kp)])
                                  D[i] = Kp[nrow(Kp),ncol(Kp)] - Kp[nrow(Kp),1:(ncol(Kp)-1)]%*%values[[i]]
                               } else{
                                  D[i] = 1
                               }
                            }
                            A = Matrix(0,nrow=np, ncol=np,sparse=TRUE)
                            A[t(matrix(unlist(indices),nrow=2))] = unlist(values)
                            B =.sparseDiagonal(np) - A
                            RiW = (.sparseDiagonal(np)*D^-0.5) %*% B
                            iW = Matrix::t(RiW) %*% RiW
                            detW = sum(log(D))
                         }
                         iWg[[ag]] = iW
                         RiWg[[ag]] = RiW
                         detWg[ag] = detW
                      }
                      rLPar[[r]] = list(iWg=iWg, RiWg=RiWg, detWg=detWg)

                   },
                   "GPP" = {
                      if(!is.null(hM$rL[[r]]$distMat)){
                         stop("predictive gaussian process not available for distance matrices")
                      }
                      s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                      sKnot = hM$rL[[r]]$sKnot
                      if (is(s, "Spatial")) {
                         dim <- ncol(coordinates(s))
                         nKnots <- nrow(coordinates(sKnot))
                         di12 <- spDists(s, sKnot)
                         di22 <- spDists(sKnot)
                      } else {
                         dim <- ncol(s)
                         nKnots <- nrow(sKnot)
                         di12 <- sqrt(Reduce("+",
                                             Map(function(i)
                                                outer(s[,i], sKnot[,i], "-")^2,
                                                seq_len(dim))))
                         di22 = as.matrix(dist(sKnot))
                      }
                      idDg = vector("list", alphaN)
                      idDW12g = vector("list", alphaN)
                      Fg = vector("list", alphaN)
                      iFg = vector("list", alphaN)
                      detDg = rep(NA,alphaN)

                      for(ag in 1:alphaN){
                         alpha = alphapw[ag,1]
                         if(alpha==0){
                            W22 = diag(nKnots)
                            W12 = matrix(0,nrow=np,ncol=nKnots)
                         } else{
                            W22 = exp(-di22/alpha)
                            W12 = exp(-di12/alpha)
                         }
                         iW22 = solve(W22)
                         D = W12%*%iW22%*%t(W12)
                         dD = 1 - diag(D)

                         liW22 = t(chol(iW22))
                         idD = 1/dD

                         tmp0 = matrix(rep(idD,nKnots),ncol=nKnots)
                         idDW12 = tmp0*W12
                         FMat = W22 + t(W12)%*%idDW12
                         iF = solve(FMat)
                         tmp2 = W12%*%liW22
                         DS = t(tmp2)%*%(tmp0*tmp2) + diag(nKnots)
                         LDS = t(chol(DS))
                         detD = sum(log(dD)) + 2*sum(log(diag(LDS)))
                         idDg[[ag]] = idD
                         idDW12g[[ag]] = idDW12
                         Fg[[ag]] = FMat
                         iFg[[ag]] = iF
                         detDg[ag] = detD
                      }
                      rLPar[[r]] = list(idDg=idDg, idDW12g=idDW12g, Fg=Fg, iFg=iFg, detDg=detDg)
                   }
            )
         }
      } else if(inherits(hM$rL[[r]],"HmscKroneckerRandomLevel",TRUE)==1){
         rLPar[[r]] = vector(mode ="list", length=length(hM$rL[[r]]$rLList))
         alphaPrior = hM$rL[[r]]$alphaPrior
         for(l in seq_len(length(hM$rL[[r]]$rLList))){
            rL = hM$rL[[r]]$rLList[[l]]
            if(rL$sDim > 0){
               alphaGrid = alphaPrior$alphaGridList[[l]]
               dfPiElem = as.factor(unlist(lapply(strsplit(as.character(hM$dfPi[,r]), hM$rL[[r]]$sepStr), function(a) a[l])))
               np = length(levels(dfPiElem))
               alphaN = length(alphaGrid)
               if(rL$spatialMethod=="Full"){
                  if(is.null(rL$distMat)){
                     s = rL$s[levels(dfPiElem),]
                     if (is(s, "Spatial"))
                        distance <- spDists(s)
                     else
                        distance = as.matrix(dist(s))
                  } else{
                     distance = rL$distMat[levels(dfPiElem),levels(dfPiElem)]
                  }
                  # Wg = vector("list", alphaN)
                  # iWg = vector("list", alphaN)
                  # RiWg = vector("list", alphaN)
                  # detWg = rep(NA, alphaN)
                  # for(ag in 1:alphaN){
                  #    alpha = alphaGrid[ag]
                  #    if(alpha==0){
                  #       W = diag(np)
                  #    } else{
                  #       W = exp(-distance/alpha)
                  #    }
                  #    RW = chol(W)
                  #    iW = chol2inv(RW)
                  #    if(rL$sDim == 1){
                  #       iW = Matrix(band(iW,-1,1), sparse=TRUE)
                  #    }
                  #    Wg[[ag]] = W
                  #    iWg[[ag]] = iW
                  #    RiWg[[ag]] = chol(iW)
                  #    detWg[ag] = 2*sum(log(diag(RW)))
                  # }
                  Wg = tf$exp(-tf$math$multiply_no_nan(tf$constant(alphaGrid^-1,tf$float64)[,NULL,NULL],tf$constant(distance,tf$float64)))
                  LWg = tfla$cholesky(Wg)
                  iWg = tfla$cholesky_solve(LWg, tf$eye(ic(nrow(distance)),batch_shape=ic(alphaN)*tf$ones(ic(1),tf$int32),dtype=tf$float64))
                  detWg = 2*tf$reduce_sum(tfm$log(tfla$diag_part(LWg)), ic(-1))
                  if(hM$rL[[r]]$etaMethod=="TF_krylov" || hM$rL[[r]]$alphaMethod=="TF_full" || hM$rL[[r]]$alphaMethod=="TF_direct_krylov"){
                     WStack = Wg #tf$cast(tf$stack(Wg), tf$float64)
                     if(l==1){
                        # nt = np
                        # indL = 2 + (nt+1)*(0:(nt-2)) - 1
                        # indD = 1 + (nt+1)*(0:(nt-1)) - 1
                        # indU = (nt+1)*(1:(nt-1)) - 1
                        # getTridiagonal = function(A) t(matrix(c(A[indU],0,A[indD],0,A[indL]),nrow(A),3))
                        # iWStack = tf$cast(tf$stack(lapply(iWg, getTridiagonal)), tf$float64)
                        iWStack = tf$reverse(tf$linalg$diag_part(iWg, k=ic(-1,1)), ic(-2)*tf$ones(ic(1),tf$int32)) # works correctly due to symmetric property
                     } else{
                        iWStack = iWg #tf$cast(tf$stack(lapply(iWg, as.matrix)), tf$float64)
                     }
                  } else{
                     WStack = iWStack = NULL
                  }
                  if(hM$rL[[r]]$etaMethod=="TF_krylov"){
                     evWStack = tf$linalg$eigh(WStack)
                     eWStack = evWStack[[1]]; vWStack = evWStack[[2]]
                  } else{
                     eWStack = eWStack = NULL
                  }
                  Wg = lapply(tf$split(Wg, ic(alphaN)), function(e) tf$squeeze(e,ic(0))$numpy())
                  RiWg = lapply(tf$split(tf$transpose(tfla$cholesky(iWg),ic(0,2,1)), ic(alphaN)), function(e) tf$squeeze(e,ic(0))$numpy())
                  iWg = lapply(tf$split(iWg, ic(alphaN)), function(e) tf$squeeze(e,ic(0))$numpy())
                  detWg = detWg$numpy()
                  if(rL$sDim == 1){
                     for(ag in 1:alphaN){
                        iWg[[ag]] = Matrix(band(iWg[[ag]],-1,1), sparse=TRUE)
                     }
                  }
                  rLPar[[r]][[l]] = list(Wg=Wg, iWg=iWg, RiWg=RiWg, detWg=detWg,
                                         WStack=WStack, iWStack=iWStack, eWStack=eWStack, vWStack=vWStack)
               } else if(rL$spatialMethod=="NNGP"){
                  stop("Hmsc.computeDataParameters: kronecker random levels with NNGP are not implemented")
               } else if(rL$spatialMethod=="NNGP"){
                  stop("Hmsc.computeDataParameters: kronecker random levels with GPP are not implemented")
               }
            }
         }
         if(hM$rL[[r]]$alphaMethod=="TF_direct_krylov" || hM$rL[[r]]$alphaMethod=="TF_direct_local_krylov"){
            krN = length(hM$rL[[r]]$rLList)
            dfPiElemLevelList = vector("list", krN)
            npElemVec = rep(NA, krN)
            for(l in seq_len(krN)){
               dfPiElem = as.factor(unlist(lapply(strsplit(as.character(hM$dfPi[,r]), hM$rL[[r]]$sepStr), function(a) a[l])))
               dfPiElemLevelList[[l]] = levels(dfPiElem)
               npElemVec[l] = length(dfPiElemLevelList[[l]])
            }
            dfTmp = expand.grid(rev(dfPiElemLevelList))[rev(1:krN)]
            allUnits = factor(do.call(function(...) paste(..., sep=hM$rL[[r]]$sepStr),dfTmp))
            indKronObs = as.numeric(factor(hM$dfPi[,r], levels=levels(allUnits)))
            gNVec = sapply(hM$rL[[r]]$alphaPrior$alphaGridList, length)
            tAlphaGridN = gNVec[1]
            sAlphaGridN = gNVec[2]
            nt = npElemVec[1]
            nx = npElemVec[2]
            obsMat = tf$transpose(tf$constant(matrix(table(c(indKronObs,1:(nt*nx)))-1,nx,nt), tf$float64))
            iKtSt = rLPar[[r]][[1]]$iWStack
            iKsSt = rLPar[[r]][[2]]$iWStack
            batchSize = hM$rL[[r]]$logDetKstBatchSize
            if(batchSize==0){
               A = (1-obsMat[1,,NULL])*(iKtSt[,NULL,2,1,NULL,NULL]*iKsSt[NULL,,,])*(1-obsMat[1,]) + tfla$diag(obsMat[1,])
               L = tfla$cholesky(A)
               logDet = 2*tf$reduce_sum(tfm$log(tfla$diag_part(L)),ic(-1))
               for(i in 2:nt){
                  print(i)
                  B = (1-obsMat[i,,NULL])*(iKtSt[,NULL,3,i,NULL,NULL]*iKsSt[NULL,,,])*(1-obsMat[i-1,])
                  CT = tfla$triangular_solve(L, tf$transpose(B,ic(0,1,3,2)))
                  A = (1-obsMat[i,,NULL])*(iKtSt[,NULL,2,i,NULL,NULL]*iKsSt[NULL,,,])*(1-obsMat[i,]) + tfla$diag(obsMat[i,])
                  H = A - tf$matmul(CT,CT,transpose_a=TRUE)
                  L = tfla$cholesky(H)
                  logDet = logDet + 2*tf$reduce_sum(tfm$log(tfla$diag_part(L)),ic(-1))
                  gc()
               }
            } else{
               batchN = ceiling((tAlphaGridN*sAlphaGridN)/batchSize)
               logDet = tf$zeros(ic(tAlphaGridN*sAlphaGridN), tf$float64)
               # outerLoopCond = function(batch,logDet) tfm$less(batch, ic(batchN))
               # outerLoopBody = function(batch,logDet){
               tmp_fun = function(batch){
                  tf$print(batch)
                  batchLen = tfm$minimum((batch+ic(1))*ic(batchSize), ic(tAlphaGridN*sAlphaGridN)) - batch*ic(batchSize)
                  batchInd = tf$range(batch*ic(batchSize), batch*ic(batchSize)+batchLen)
                  # batchInd = tf$reshape(batchInd, batchInd$shape[1])
                  tGridInd = tfm$floordiv(batchInd, ic(sAlphaGridN))
                  sGridInd = tfm$floormod(batchInd, ic(sAlphaGridN))
                  iKt = tf$gather(iKtSt, tGridInd)
                  iKs = tf$gather(iKsSt, sGridInd)
                  i = tf$constant(ic(0), tf$int32)
                  A = (tf$constant(1,tf$float64)-obsMat[i,,NULL])*(iKt[,tf$constant(ic(1),tf$int32),i,NULL,NULL]*iKs)*(tf$constant(1,tf$float64)-obsMat[i,]) + tfla$diag(obsMat[i,])
                  L = tfla$cholesky(A)
                  # L = tf$Variable(tfla$cholesky(A))
                  logDetBatch = tf$constant(2,tf$float64)*tf$reduce_sum(tfm$log(tfla$diag_part(L)),ic(-1))
                  innerLoopCond = function(i,logDetBatch,L) tfm$less(i, ic(nt))
                  innerLoopBody = function(i,logDetBatch,L){
                     B = (tf$constant(1,tf$float64)-obsMat[i,,NULL])*(iKt[,tf$constant(ic(2),tf$int32),i,NULL,NULL]*iKs)*(tf$constant(1,tf$float64)-obsMat[i-ic(1),])
                     CT = tfla$triangular_solve(L, tf$transpose(B,ic(0,2,1)))
                     A = (tf$constant(1,tf$float64)-obsMat[i,,NULL])*(iKt[,tf$constant(ic(1),tf$int32),i,NULL,NULL]*iKs)*(tf$constant(1,tf$float64)-obsMat[i,]) + tfla$diag(obsMat[i,])
                     H = A - tf$matmul(CT,CT,transpose_a=TRUE)
                     L = tfla$cholesky(H)
                     # L$assign(tfla$cholesky(H))
                     logDetBatch = logDetBatch + tf$constant(2,tf$float64)*tf$reduce_sum(tfm$log(tfla$diag_part(L)),ic(-1))
                     return(c(i+ic(1),logDetBatch,L))
                  }
                  innerLoopInit = c(i+ic(1),logDetBatch,L)
                  innerLoopRes = tf$while_loop(innerLoopCond, innerLoopBody, innerLoopInit)
                  logDetBatch = innerLoopRes[[2]]
                  return(logDetBatch)
               }
               tmp_fun_tf = tf_function(tmp_fun)
               for(batch_r in 1:batchN){
                  batch = tf$constant(ic(batch_r-1), tf$int32)
                  logDetBatch = tmp_fun_tf(batch)
                  batchLen = tfm$minimum((batch+ic(1))*ic(batchSize), ic(tAlphaGridN*sAlphaGridN)) - batch*ic(batchSize)
                  batchInd = tf$range(batch*ic(batchSize), batch*ic(batchSize)+batchLen)
                  logDet = tf$tensor_scatter_nd_add(logDet, batchInd[,NULL], logDetBatch)
                  # return(c(batch+ic(1),logDet))
               }
               # outerLoopInit = c(tf$constant(ic(0),tf$int32), logDet)
               # outerLoopRes = tf$while_loop(outerLoopCond, outerLoopBody, outerLoopInit)
               # logDet = outerLoopRes[[2]]
               logDet = tf$reshape(logDet, ic(tAlphaGridN,sAlphaGridN))
            }
            logDetKt = tf$constant(rLPar[[r]][[1]]$detWg, tf$float64)
            logDetKs = tf$constant(rLPar[[r]][[2]]$detWg, tf$float64)
            rLPar[[r]]$logDetK = (nx*logDetKt[,NULL] + nt*logDetKs[NULL,] + logDet)$numpy()
         }
      }
   }
   parList$Qg = Qg
   parList$iQg = iQg
   parList$RQg = RQg
   parList$detQg = detQg
   parList$rLPar = rLPar

   # print(rLPar)
   return(parList)
}

