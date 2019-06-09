#' @title alignPosterior
#'
#' @description Aligns posterior in terms of variable susceptible for label switching
#'
#' @importFrom stats cor
#' @importFrom abind abind
#' @export

alignPosterior=function(hM){
   bind0 = function(...){
      abind(...,along=0)
   }
   ncRRR = hM$ncRRR
   ncNRRR = hM$ncNRRR
   nr = hM$nr

   if(nr>0){
      for (r in 1:nr){
          for (i in 1:length(hM$postList)){
            cpL= hM$postList[[i]]
            valList=lapply(cpL, function(a) a[["Lambda"]][[1]])
            if (i==1){
               val = do.call(bind0, valList)
               valMean = colMeans(val)
            }
            s=NULL
            for(k in 1:dim(valMean)[1]){
               s=rbind(s,sign(abind(lapply(valList, function(a) cor(as.vector(valMean[k,]),as.vector(a[k,]))))))
            }
            for(j in 1:length(cpL)){
               for(k in 1:dim(valMean)[1]){
                  if(s[k,j]<0){
                     cpL[[j]]$Lambda[[1]][k,] = -cpL[[j]]$Lambda[[1]][k,]
                     cpL[[j]]$Eta[[1]][,k] = -cpL[[j]]$Eta[[1]][,k]
                  }
               }
            }
            hM$postList[[i]] = cpL
         }
      }
   }
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
