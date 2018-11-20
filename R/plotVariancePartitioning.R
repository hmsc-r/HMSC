#' @title plotVariancePartitioning
#'
#' @description Plots ...
#' @param VP
#' @param ... other parameters passed to the barplot function
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

plotVariancePartitioning = function(hM, VP, ...){

   ng = dim(VP$vals)[1]
   leg = VP$groupnames
   for(r in 1:hM$nr){
      leg = c(leg, paste("random: ", hM$levelNames[r], sep=""))
   }
   means = round(100*rowMeans(VP$vals), 1)
   for(i in 1:ng){
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]), ")", sep="")
   }
   # par(mfrow=c(1,1), mar=c(5,5,4,8)) what is the point of this?
   mainTitle = paste("Variance partitioning. R2(traits) = ", round(VP$traitR2,2), ".", sep="")
   barplot(VP$vals, main=mainTitle, xlab="Species", ylab="Variance proportion", legend=leg,
     col=heat.colors(ng,alpha=1), ...)
}
