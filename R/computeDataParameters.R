#' @title computeDataParameters
#'
#' @description Computes initial values before the sampling starts
#'
#' @param hM a fitted \code{Hmsc} model object
#'
#' @return a list including pre-computed matrix inverses and determinants (for phylogenetic and spatial random effects) needed in MCMC sampling
#'
#' @importFrom stats dist
#' @importFrom FNN get.knn
#' @importFrom Matrix .sparseDiagonal t solve
#'
#' @export


computeDataParameters = function(hM){
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
      if(hM$rL[[r]]$sDim > 0){
         alphapw = hM$rL[[r]]$alphapw
         np = hM$np[r]
         alphaN = nrow(alphapw)
         switch(hM$rL[[r]]$spatialMethod,
                "Full" = {
                   if(is.null(hM$rL[[r]]$distMat)){
                      s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                      distance = as.matrix(dist(s))
                   } else{
                      distance = hM$rL[[r]]$distMat[levels(hM$dfPi[,r]),levels(hM$dfPi[,r])]
                   }
                   Wg = array(NA, c(np,np,alphaN))
                   iWg = array(NA, c(np,np,alphaN))
                   RiWg = array(NA, c(np,np,alphaN))
                   detWg = rep(NA, alphaN)
                   for(ag in 1:alphaN){
                      alpha = alphapw[ag,1]
                      if(alpha==0){
                         W = diag(np)
                      } else{
                         W = exp(-distance/alpha)
                      }
                      RW = chol(W)
                      iW = chol2inv(RW)

                      Wg[,,ag] = W
                      iWg[,,ag] = iW
                      RiWg[,,ag] = chol(iW)
                      detWg[ag] = 2*sum(log(diag(RW)))
                   }
                   rLPar[[r]] = list(Wg=Wg, iWg=iWg, RiWg=RiWg, detWg=detWg)
                },
                "NNGP" = {
                   if(is.null(hM$rL[[r]]$nNeighbours)){
                      hM$rL[[r]]$nNeighbours = 10
                   }
                   if(!is.null(hM$rL[[r]]$distMat)){
                      stop("computeDataParameters: Nearest neighbours not available for distance matrices")
                   }
                   iWg = list()
                   RiWg = list()
                   detWg = rep(NA,alphaN)
                   s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                   indNN = get.knn(s,k=hM$rL[[r]]$nNeighbours)[[1]] #This uses FNN package, may be updated later to incorporate distance matrices
                   indNN = t(apply(indNN,1,sort,decreasing=FALSE))
                   indices = list()
                   distList = list()
                   for(i in 2:np){
                      ind = indNN[i,]
                      ind = ind[ind<i]
                      if(!is.na(ind[1])){
                         indices[[i]] = rbind(i*rep(1,length(ind)),ind)
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
                      stop("computeDataParameters: predictive gaussian process not available for distance matrices")
                   }
                   s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                   sKnot = hM$rL[[r]]$sKnot
                   dim = ncol(s)
                   nKnots = nrow(sKnot)
                   di = matrix(0,nrow=np,ncol=nKnots)
                   for(i in 1:dim){
                      xx1 = matrix(rep(s[,i],nKnots),ncol=nKnots)
                      xx2 = matrix(rep(sKnot[,i],np),ncol=np)
                      dx = xx1 - t(xx2)
                      di = di+dx^2
                   }
                   di12 = sqrt(di)

                   di22 = as.matrix(dist(sKnot))

                   idDg = matrix(NA,nrow=np,ncol=alphaN)
                   idDW12g = array(NA, c(np,nKnots,alphaN))
                   Fg = array(NA, c(nKnots,nKnots,alphaN))
                   iFg = array(NA, c(nKnots,nKnots,alphaN))
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
                      idDg[,ag] = idD
                      idDW12g[,,ag] = idDW12
                      Fg[,,ag] = FMat
                      iFg[,,ag] = iF
                      detDg[ag] = detD
                   }
                   rLPar[[r]] = list(idDg=idDg, idDW12g=idDW12g, Fg=Fg, iFg=iFg, detDg=detDg)
                }
                )
      }
   }
   parList$Qg = Qg
   parList$iQg = iQg
   parList$RQg = RQg
   parList$detQg = detQg
   parList$rLPar = rLPar

   return(parList)
}

