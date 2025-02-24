#' @title coralSplitToBatches
#'
#' @description Splits Hmsc data into common and rare parts for CORAL analysis
#'
#' @param Y fitted \code{Hmsc}-class object
#' @param TrData arg2
#' @param C.common.rare arg3
#' @param batchN arg4
#'
#' @return
#' A list containing three lists with splits of the input arguments
#'#'
#' @export

coralSplitToBatches = function(Y, TrData, C.common.rare, batchN=1){
   spNames = colnames(Y)
   ns = length(spNames)
   spNamesList = split(spNames, ceiling(seq_along(spNames)/ceiling(ns/batchN)))
   YList = vector("list", batchN)
   for(bInd in seq_len(batchN)){
      YList[[bInd]] = Y[,match(spNamesList[[bInd]], colnames(Y))]
   }
   TrDataList = vector("list", batchN)
   if(!is.null(TrData)){
      for(bInd in seq_len(batchN)){
         TrDataList[[bInd]] = TrData[match(spNamesList[[bInd]], rownames(TrData)), ]
      }
   }
   CcrList = vector("list", batchN)
   if(!is.null(C.common.rare)){
      for(bInd in seq_len(batchN)){
         CcrList[[bInd]] = C.common.rare[, match(spNamesList[[bInd]], colnames(C.common.rare))]
      }
   }
   return(list(Y=YList, TrData=TrDataList, C.common.rare=CcrList))
}
