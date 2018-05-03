#' @title Hmsc$computeR2
#'
#' @description Computes explainatory R2 measure based on predicted species abundance matrix
#' @param predY
#'
#' @examples
#'

computeR2 = function(predY){
   m = self
   R2 = matrix(NA,m$ns,1)
   for (i in 1:m$ns){
      R2[i,]=mean(predY[m$Y[,i]==1,i])-mean(predY[m$Y[,i]==0,i])
   }
   return(R2)
}

Hmsc$set("public", "computeR2", computeR2, overwrite=TRUE)
