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

   for(r in seq_len(nr)){
      nfVec = unlist(lapply(hM$postList, function(postChain) dim(postChain[[1]][["Lambda"]][[r]])[1]))
      nfMax = max(nfVec)
      cIndNfMax = which.max(nfVec)
      cpL = hM$postList[[cIndNfMax]]
      LambdaPostList = lapply(cpL, function(a) a[["Lambda"]][[r]])
      LambdaPostArray = do.call(bind0, LambdaPostList)
      LambdaPostMean = colMeans(LambdaPostArray)
      for(cInd in 1:length(hM$postList)){
         cpL = hM$postList[[cInd]]
         LambdaPostList = lapply(cpL, function(a) a[["Lambda"]][[r]])
         nf = nrow(LambdaPostList[[1]])
         s = NULL
         for(k in 1:nf){
            if(hM$ns > 1 || hM$rL[[r]]$xDim > 1){
               s = rbind(s,sign(abind(lapply(LambdaPostList, function(a) cor(as.vector(LambdaPostMean[k,]),as.vector(a[k,]))))))
            } else
               s = rbind(s,abind(lapply(LambdaPostList, function(a) sign(LambdaPostMean[k,])*sign(a[k,]))))
         }
         for(j in 1:length(cpL)){
            for(k in 1:nf){
               if(s[k,j]<0){
                  cpL[[j]]$Lambda[[r]][k,] = -cpL[[j]]$Lambda[[r]][k,]
                  cpL[[j]]$Eta[[r]][,k] = -cpL[[j]]$Eta[[r]][,k]
               }
            }
         }
         hM$postList[[cInd]] = cpL
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
