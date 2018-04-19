updateAlpha = function(Eta ,rL, rLPar){
   nr = length(Eta)

   Alpha = vector("list", nr)
   for(r in 1:nr){
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

         tmp1 = array(NA,dim=c(np,nf,gN))
         for(g in 1:gN){ # this cycle should be replaced with tensor operation
            tmp1[,,g] = RiWg[,,g]%*%eta
         }
         tmpMat = t(apply(tmp1^2, c(2,3), sum))
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


