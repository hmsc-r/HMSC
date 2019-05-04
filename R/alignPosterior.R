#' @title alignPosterior
#'
#' @description Aligns posterior in terms of variable susceptible for label switching
#'
#' @export

alignPosterior=function(hM){
   bind0 = function(...){
      abind(...,along=0)
   }
   ncDR = hM$ncDR
   ncNDR = hM$ncNDR
   if(ncDR>0){
      for (i in 1:length(hM$postList)){
         cpL= hM$postList[[i]]
         valList=lapply(cpL, function(a) a[["wDR"]])
         if (i==1){
            val = do.call(bind0, valList)
            valMean = colMeans(val)
         }
         s=NULL
         for(k in 1:ncDR){
            s = rbind(s,sign(abind(lapply(valList, function(a) cor(as.vector(valMean[k,]),as.vector(a[k,]))))))
         }
         for(j in 1:length(cpL)){
            for(k in 1:ncDR){
               if(s[k,j]<0){
                  cpL[[j]]$wDR[k,] = -cpL[[j]]$wDR[k,]
                  cpL[[j]]$Beta[ncNDR+k,] = -cpL[[j]]$Beta[ncNDR+k,]
                  cpL[[j]]$Gamma[ncNDR+k,] = -cpL[[j]]$Gamma[ncNDR+k,]
                  cpL[[j]]$V[ncNDR+k,] = -cpL[[j]]$V[ncNDR+k,]
                  cpL[[j]]$V[,ncNDR+k] = -cpL[[j]]$V[,ncNDR+k]
               }
            }
         }
         hM$postList[[i]] = cpL
      }
   }
   return(hM)
}
