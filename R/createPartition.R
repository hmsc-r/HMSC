#' @title createPartition
#'
#' @description Constructs a partition vector given the number of folds and column of study design
#' @param nfolds number of cross-validation folds
#' @param column name or index of the column in studyDesign matrix, corresponding to the level, units of which are splitted to folds
#'
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export

createPartition = function(hM, nfolds=10, column=NULL){
   if(ncol(hM$studyDesign)>0 && !is.null(column)){
      np = length(unique(hM$studyDesign[,column]))
      if(np < nfolds){
         stop("HMSC.createPartition: nfolds must be no bigger than number of units in specified random level")
      }
      tmp1 = data.frame(col=unique(hM$studyDesign[,column]), partition=sample(rep(1:nfolds,ceiling(np/nfolds)),np))
      colnames(tmp1)[1] = colnames(hM$studyDesign[,column,drop=FALSE])
      tmp2 = merge(hM$studyDesign[,column,drop=FALSE], tmp1, all.x=TRUE, sort=FALSE)
      part = tmp2$partition
   }
   else{
      part=sample(rep(1:nfolds,ceiling(hM$ny/nfolds)),hM$ny)
   }
   return(part)
}
