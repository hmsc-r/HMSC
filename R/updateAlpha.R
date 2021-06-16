#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#' @importFrom matrixStats colMaxs
#' @importFrom pracma ndims
#' @importFrom plyr mdply
#'
updateAlpha = function(Eta ,rL, rLPar){
   nr = length(rL)
   Alpha = vector("list", nr)
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
         Alpha[[r]] = matrix(NA, nf, length(rL[[r]]$rLList))
         dfPiElemLevelList = vector("list", length(rL[[r]]$rLList))
         npElemVec = rep(NA, length(rL[[r]]$rLList))
         for(l in seq_len(length(rL[[r]]$rLList))){
            dfPiElem = as.factor(unlist(lapply(strsplit(as.character(m$dfPi[,r]), "="), function(a) a[l])))
            dfPiElemLevelList[[l]] = levels(dfPiElem)
            npElemVec[l] = length(dfPiElemLevelList[[l]])
         }
         dfTemp = expand.grid(rev(dfPiElemLevelList))[rev(1:length(rL[[r]]$rLList))]
         allUnits = as.factor(mdply(dfTemp, paste, sep=rL[[r]]$sepStr)[,length(dfPiElemLevelList)+1])
         ind = as.numeric(factor(m$dfPi[,r], levels=levels(allUnits)))
         if(length(allUnits) == np[r]){
            etaArray = array(eta,c(npElemVec,nf))
            iWgEtaList = vector("list", length(rL[[r]]$rLList))
            for(l in seq_len(length(rL[[r]]$rLList))){
               if(rL[[r]]$rLList[[l]]$sDim == 0){
                  iWgEtaList[[l]] = array(etaArray, c(1,dim(etaArray)))
               } else{
                  gN = length(rL[[r]]$alphaPrior$alphaGridList[[l]])
                  iWg = rLPar[[r]][[l]]$iWg
                  iWgEta = array(NA, c(gN,dim(etaArray)))
                  dimInd = c(l,c(1:ndims(etaArray))[-l])
                  dimIndOrd = order(dimInd)
                  etaMatFlat = matrix(aperm(etaArray, dimInd), dim(etaArray)[l])
                  for(ag in seq_len(gN)){
                     iWgEtaMatFlat = iWg[[ag]]%*%etaMatFlat
                     iWgEta[((ag-1)*length(etaArray)+1):(ag*length(etaArray))] = aperm(array(iWgEtaMatFlat, dim(etaArray)[dimInd]), dimIndOrd)
                  }
                  iWgEtaList[[l]] = iWgEta
               }
            }
            #TODO:GT potential for >2 kronecker elements
            gN1 = length(rL[[r]]$alphaPrior$alphaGridList[[1]])
            gN2 = length(rL[[r]]$alphaPrior$alphaGridList[[2]])
            tmp1 = aperm(array(rep(iWgEtaList[[1]], each=gN2), c(gN2,gN1,dim(etaArray))), c(2,1,2+1:ndims(etaArray)))
            tmp2 = array(rep(iWgEtaList[[2]], each=gN1), c(gN1,gN2,dim(etaArray)))
            EtaiWEta = tmp1 * tmp2
            qFArray = apply(EtaiWEta, c(1,2,2+ndims(etaArray)), sum)
            tmp1 = matrix(rLPar[[r]][[1]]$detWg, gN1,gN2)
            tmp2 = matrix(rLPar[[r]][[2]]$detWg, gN1,gN2, byrow=TRUE)
            logDetWArray = array(tmp1+tmp2, c(gN1,gN2,nf))
            logPriorArray = array(log(rL[[r]]$alphaPrior$alphaProb), c(gN1,gN2,nf))
            logProbArray = logPriorArray - 0.5*logDetWArray - 0.5*qFArray #dimension-based const coefficients are omitted in this expr
            for(h in 1:nf){
               logProbVec = as.vector(logProbArray[,,h])
               logProbVec = logProbVec - max(logProbVec)#log(sum(exp(logProbVec)))
               alphaFlatInd = sample(gN1*gN2, 1, prob=exp(logProbVec))
               alphaInd1 = ((alphaFlatInd-1) %% gN1) + 1
               alphaInd2 = ((alphaFlatInd-1) %/% gN1) + 1
               Alpha[[r]][h,] = c(alphaInd1,alphaInd2)
            }
         } else{
            stop("updateAlpha: not implemented yet")
         }
      }
   }
   return(Alpha)
}


