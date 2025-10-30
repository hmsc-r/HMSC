#' @title coralSplitToBatches
#'
#' @description Splits CORAL data objects that contain all rare species into smaller batches for distributed computation
#'
#' @param Y community matrix of rare species
#' @param TrData trait matrix (dataframe) for rare species
#' @param C.common.rare phylogeny similarity matrix between common and rare species
#' @param batchN target number of batches
#'
#' @return
#' A named list of three elements \code{Y}, \code{TrData} and \code{C.common.rare},
#' each of which is a list of length \code{batchN} with splits of the corresponding input argument
#'
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
