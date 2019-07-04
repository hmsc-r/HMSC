#' @title plotVariancePartitioning
#'
#' @description Plots the results of variance partitioning of a HMSC model produced by `HMSC::computeVariancePartitioning` as a barplot
#' @param hM A fitted HMSC model
#' @param VP A HMSC variance partitioning produced by `HMSC::computeVariancePartitioning`
#' @param ... other parameters passed to the barplot function
#'
#' @examples
#' # Plot how the explained variance of a previously fitted model is partitioned
#' VP = computeVariancePartitioning(TD$m)
#' plotVariancePartitioning(TD$m,VP)
#'
#' @importFrom graphics barplot
#' @importFrom grDevices heat.colors
#'
#' @export

plotVariancePartitioning=function (hM,VP, cols=NULL,...)
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

   mainTitle=substitute("Variance partitioning")
   barplot(VP$vals, main = mainTitle, xlab= "Species", ylab = "Variance proportion", las = 1,
           legend = leg, col = cols,...)
#   mtext("Species", 1,line = 1)
}
