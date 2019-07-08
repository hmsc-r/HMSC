#' @title plotVariancePartitioning
#'
#' @description Plots the results of variance partitioning of a \code{Hmsc} model produced by
#' \code{\link{computeVariancePartitioning}} as a barplot
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param VP a Hmsc variance partitioning object produced by \code{\link{computeVariancePartitioning}}
#' @param cols colors of the barplot
#' @param ... additional parameters passed to the barplot function
#'
#' @examples
#' # Plot how the explained variance of a previously fitted model is partitioned
#' VP = computeVariancePartitioning(TD$m)
#' plotVariancePartitioning(TD$m, VP)
#'
#' @importFrom graphics barplot
#' @importFrom grDevices heat.colors
#'
#' @export

plotVariancePartitioning=function (hM, VP, cols=NULL, ...)
{
   ng = dim(VP$vals)[1]
   if(is.null(cols)){
      cols = heat.colors(ng, alpha = 1)
   }
   leg = VP$groupnames
   for (r in 1:hM$nr) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
   }
   means = round(100 * rowMeans(VP$vals), 1)
   for (i in 1:ng) {
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]),
                     ")", sep = "")
   }

   mainTitle = "Variance partitioning"
   barplot(VP$vals, main = mainTitle, xlab= "Species", ylab = "Variance proportion", las = 1,
           legend = leg, col = cols,...)
#   mtext("Species", 1,line = 1)
}
