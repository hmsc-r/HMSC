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
                   for(g in 1:gN){ # this cycle should be replaced with tensor operation
                      tmp1 = RiWg[[g]]%*%eta
                      tmpMat[g,] = Matrix::colSums(tmp1^2)
                   }
                }
         )
         # tmp1 = mul.tensor(to.tensor(RiWg),"I2", to.tensor(rep(eta,gN),c(I1=np,I2=nf,I3=gN)),"I1", by="I3")

         # tmpMat = t(apply(tmp1^2, c(2,3), sum))
         for(h in 1:nf){
            v = tmpMat[,h]
            like = log(alphapw[,2]) - 0.5*detWg - 0.5*v
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


