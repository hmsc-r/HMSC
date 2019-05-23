#' @title computeInitialParameters
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'
#' @return
#'
#' @importFrom stats dist
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
                   Wg = list()
                   iWg = list() #slow??
                   RiWg = list()
                   detWg = rep(NA,alphaN)
                   s = hM$rL[[r]]$s[levels(hM$dfPi[,r]),]
                   indNN = get.knn(s,k=hM$rL[[r]]$nNeighbours)[[1]] #This uses FNN package, may be updated later to incorporate distance matrices
                   indNN = t(apply(indNN,1,sort,decreasing=FALSE))
                   indices = list() #slow?
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
                         detW = 1
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
                      Wg[[ag]] = Matrix::solve(iW,.sparseDiagonal(np))
                      iWg[[ag]] = iW
                      RiWg[[ag]] = RiW
                      detWg[ag] = detW
                   }
                   rLPar[[r]] = list(Wg = Wg, iWg=iWg, RiWg=RiWg, detWg=detWg)
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

