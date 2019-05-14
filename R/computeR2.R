#' @title Hmsc$computeR2
#'
#' @description Computes explainatory R2 measure based on predicted species abundance matrix
#' @param predY vector of predictions from the trained model
#'
#' @return
#'
#' @examples
#'
#' @importFrom stats cor
#' @export

computeR2 = function(hM, predY){
   R2 = matrix(NA,hM$ns,1)
   for (i in 1:hM$ns){
      if (hM$distr[i,1]==1){
         R2[i,]=(cor(predY[,i],hM$Y[,i]))^2
      }
      if (hM$distr[i,1]==2){
         R2[i,]=mean(predY[hM$Y[,i]==1,i])-mean(predY[hM$Y[,i]==0,i])
      }
   }
   return(R2)
}
