#' @importFrom Matrix bdiag Diagonal sparseMatrix t Matrix
#'
updateAlpha = function(Eta ,rL, rLPar){
   nr = length(rL)

   Alpha = vector("list", nr)
   for(r in seq_len(nr)){
      eta = Eta[[r]]
      np = nrow(eta)
      nf = ncol(eta)
      if(rL[[r]]$sDim > 0){
         iWg = rLPar[[r]]$iWg
         RiWg = rLPar[[r]]$RiWg
         detWg = rLPar[[r]]$detWg
         alphapw = rL[[r]]$alphapw
         gN = nrow(alphapw)

         Alpha[[r]] = rep(NA, nf)

         switch(rL[[r]]$spatialMethod,
                'Full' = {
                   tmp1 = array(NA,dim=c(np,nf,gN))
                   for(g in 1:gN){ # this cycle should be replaced with tensor operation
                      tmp1[,,g] = RiWg[,,g]%*%eta
                   }
                   tmpMat = t(colSums(tmp1^2))
                },
                'NNGP' = {
                   tmpMat = matrix(NA,nrow=gN,ncol=nf)
                   for(g in 1:gN){
                      tmp1 = RiWg[[g]]%*%eta
                      tmpMat[g,] = Matrix::colSums(tmp1^2)
                   }
                },
                'GPP' = {
                   iFg = rLPar[[r]]$iFg
                   nK = nrow(iFg)
                   idDW12g = rLPar[[r]]$idDW12g
                   idDg = rLPar[[r]]$idDg
                   detDg = rLPar[[r]]$detDg
                   tmpMat2 = array(NA,dim=c(nf,nK,gN))
                   tmpMat3 = array(NA,dim=c(nf,nK,gN))
                   tmpMat4 = array(NA,dim=c(nf,nf,gN))
                   for(g in 1:gN){
                      tmpMat2[,,g] = crossprod(eta,idDW12g[,,g])
                      tmpMat3[,,g] = tmpMat2[,,g]%*%iFg[,,g]
                      tmpMat4[,,g] = matrix(tmpMat3[,,g],ncol=nK,nrow=nf)%*%t(matrix(tmpMat2[,,g],,ncol=nK,nrow=nf)) #This extra matrix wrap is neccessary in case there is only one LV
                   }
                }
         )
         # tmp1 = mul.tensor(to.tensor(RiWg),"I2", to.tensor(rep(eta,gN),c(I1=np,I2=nf,I3=gN)),"I1", by="I3")

         # tmpMat = t(apply(tmp1^2, c(2,3), sum))
         for(h in 1:nf){
            switch(rL[[r]]$spatialMethod,
                   'Full' = {
                      v = tmpMat[,h]
                      like = log(alphapw[,2]) - 0.5*detWg - 0.5*v
                   },
                   'NNGP' = {
                      v = tmpMat[,h]
                      like = log(alphapw[,2]) - 0.5*detWg - 0.5*v
                   },
                   'GPP' = {
                      tmp = rep(NA,gN)
                      for(ag in 1:gN){
                         if(alphapw[ag,1] == 0){
                            tmp[ag] = t(eta[,h])%*%eta[,h]
                         } else {
                            tmp1 = t(eta[,h])%*%(idDg[,ag]*eta[,h])
                            tmp[ag] = tmp1 - tmpMat4[h,h,ag]
                         }
                      }
                      like = log(alphapw[,2]) - 0.5*detDg - 0.5*tmp
                   }
            )
            like = exp(like - max(like))
            like = like / sum(like)
            Alpha[[r]][h] = sample.int(gN, size=1, prob=like)
         }
      } else{
         Alpha[[r]] = rep(1, nf)
      }
   }
   return(Alpha)
}


