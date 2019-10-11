#' @title createPartition
#'
#' @description Constructs a partition vector given the number of folds and column of study design
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param nfolds number of cross-validation folds
#' @param column name or index of the column in the \code{studyDesign} matrix, corresponding to the level for which units are splitted into folds
#'
#' @return a vector describing the fold of each sampling unit
#'
#' @examples
#' # Create 3 folds for a HMSC object
#' partition = createPartition(TD$m, nfolds = 3)
#' @export

createPartition = function(hM, nfolds=10, column=NULL){
   if(ncol(hM$studyDesign)>0 && !is.null(column)){
      np = length(unique(hM$studyDesign[,column]))
      if(np < nfolds){
         stop("HMSC.createPartition: nfolds cannot exceed the number of units in the specified random level")
      }
      level = hM$studyDesign[,column]
      levels = unique(hM$studyDesign[,column])
      partition1 = sample(rep(1:nfolds,ceiling(np/nfolds)),np)
      part= rep(NA,hM$ny)
      for (i in 1:nfolds){
         sel = levels[partition1==i]
         for (lev in sel){
            part[level==lev] = i
         }
      }
   }
   else{
      part=sample(rep(1:nfolds,ceiling(hM$ny/nfolds)),hM$ny)
   }
   return(part)
}
