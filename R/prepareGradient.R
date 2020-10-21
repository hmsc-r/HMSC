#' @title prepareGradient
#'
#' @description prepares a user-made environmental and/or spatial gradient to be used for prediction
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param XDataNew a dataframe of the new XData
#' @param sDataNew a list of the new sData
#'
#'
#' @return a named list with members \code{XDataNew}, \code{studyDesignNew} and \code{rLNew}
#'
#' @details
#' The dataframe \code{XDataNew} is the for output as for input. The main purpose of this function is to prepare
#' the study design and random levels so that predictions can be made with the \code{\link{predict}} function.
#' Note that the difference between \code{constructGradient} and \code{prepareGradient} is that while
#' \code{prepareGradient} takes as input the new environmental and spatial data, \code{constructGradient}
#' generates those data to represent a new environmental gradient.
#'
#'
#' @seealso
#' \code{\link{constructGradient}}, \code{\link{predict}}
#'
#' @export

prepareGradient = function(hM, XDataNew, sDataNew){
   #xDataNew needs to be added and the function needs to be more extensively tested
   nyNew = NROW(XDataNew)
   dfPiNew = matrix(NA,nyNew,hM$nr)
   colnames(dfPiNew) = hM$rLNames
   rLNew = vector("list", hM$nr)
   names(rLNew) = hM$rLNames
   for (r in seq_len(hM$nr)){
      rL1 = hM$rL[[r]]
      if (rL1$sDim==0){
         dfPiNew[,r] = sprintf('new_unit',1:(nyNew))
         unitsAll = c(rL1$pi,dfPiNew[,r])
         rL1$pi = unitsAll
         rL1$N = rL1$N+1 # '+1' - shouldn't this be length(unitsAll)?
      } else {
         index = which(names(sDataNew)==hM$rLNames)
         xyNew = sDataNew[[index]]
         nxyNew = NROW(xyNew)
         if (nxyNew != nyNew) # or fails mystically later
            stop("new XData and sData must have same numbers of rows")
         dfPiNew[,r] = sprintf('new_spatial_unit%.6d',1:nxyNew)
         unitsAll = c(rL1$pi,dfPiNew[,r])
         rL1$pi = unitsAll
         row.names(xyNew) = dfPiNew[,r]
         xyOld = rL1$s
         rL1$s = rbind(xyOld, xyNew)
      }
      rLNew[[r]] = rL1
   }
   studyDesignNew = as.data.frame(dfPiNew, stringsAsFactors = TRUE)
   Gradient = list(XDataNew=XDataNew, studyDesignNew=studyDesignNew, rLNew=rLNew)
   return(Gradient)
}
