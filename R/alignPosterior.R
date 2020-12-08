#' @title alignPosterior
#'
#' @description Aligns posterior in terms of variables susceptible to label switching
#'
#' @param hM a fitted \code{Hmsc} model object
#'
#' @return an \code{Hmsc} model object that is identical to the input except the posterior being aligned
#'
#' @examples
#' # Align the posterior for a previously fitted HMSC model
#' m = alignPosterior(TD$m)
#'
#' @importFrom stats cor sd
#' @importFrom abind abind
#'
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
               s = rbind(s,sign(abind(lapply(LambdaPostList, function(a){
                     if(sd(as.vector(LambdaPostMean[k,]))>0 && sd(as.vector(a[k,]))>0){
                        return(cor(as.vector(LambdaPostMean[k,]),as.vector(a[k,])))
                     } else
                        return(0)
                  }))))
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
            if(nf < nfMax){
               LambdaAddDim = dim(cpL[[j]]$Lambda[[r]])
               LambdaAddDim[1] = nfMax-nf
               cpL[[j]]$Lambda[[r]] = abind(cpL[[j]]$Lambda[[r]], array(0,LambdaAddDim), along=1)
               cpL[[j]]$Psi[[r]] = abind(cpL[[j]]$Psi[[r]], array(0,LambdaAddDim), along=1)
               DeltaAddDim = dim(cpL[[j]]$Delta[[r]])
               DeltaAddDim[1] = nfMax-nf
               cpL[[j]]$Delta[[r]] = abind(cpL[[j]]$Delta[[r]], array((hM$rL[[r]]$a2-1)/hM$rL[[r]]$b2,DeltaAddDim), along=1)
               cpL[[j]]$Eta[[r]] = abind(cpL[[j]]$Eta[[r]], matrix(0,nrow(cpL[[j]]$Eta[[r]]),nfMax-nf), along=2)
               if(hM$rL[[r]]$sDim > 0)
                  cpL[[j]]$Alpha[[r]] = abind(cpL[[j]]$Alpha[[r]], rep(1,nfMax-nf), along=1)
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
