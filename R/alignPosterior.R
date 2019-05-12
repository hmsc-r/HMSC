#' @title alignPosterior
#'
#' @description Aligns posterior in terms of variable susceptible for label switching
#'
#' @export

alignPosterior=function(hM){
   bind0 = function(...){
      abind(...,along=0)
   }
   ncRRR = hM$ncRRR
   ncNRRR = hM$ncNRRR
   if(ncRRR>0){
      for (i in 1:length(hM$postList)){
         cpL= hM$postList[[i]]
         valList=lapply(cpL, function(a) a[["wRRR"]])
         if (i==1){
            val = do.call(bind0, valList)
            valMean = colMeans(val)
         }
         s=NULL
         for(k in 1:ncRRR){
            s = rbind(s,sign(abind(lapply(valList, function(a) cor(as.vector(valMean[k,]),as.vector(a[k,]))))))
         }
         for(j in 1:length(cpL)){
            for(k in 1:ncRRR){
               if(s[k,j]<0){
                  cpL[[j]]$wRRR[k,] = -cpL[[j]]$wRRR[k,]
                  cpL[[j]]$Beta[ncNRRR+k,] = -cpL[[j]]$Beta[ncNRRR+k,]
                  cpL[[j]]$Gamma[ncNRRR+k,] = -cpL[[j]]$Gamma[ncNRRR+k,]
                  cpL[[j]]$V[ncNRRR+k,] = -cpL[[j]]$V[ncNRRR+k,]
                  cpL[[j]]$V[,ncNRRR+k] = -cpL[[j]]$V[,ncNRRR+k]
               }
            }
         }
         hM$postList[[i]] = cpL
      }
   }
   return(hM)
}
