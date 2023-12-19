computePredictor.HmscRandomLevel = function(rL,Eta,Lambda,pi=c(1:nrow(Eta)),dfpi=NULL){
   if(rL$xDim == 0){
      LRan = Eta[pi,]%*%Lambda
   } else{
      LRan = matrix(0,length(pi),ncol(Lambda))
      for(k in 1:rL$xDim)
         LRan = LRan + (Eta[pi,]*rL$x[as.character(dfpi),k]) %*% Lambda[,,k]
   }
   return(LRan)
}
