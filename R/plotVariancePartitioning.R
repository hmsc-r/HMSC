#' @title plotVariancePartitioning
#'
#' @description Plots the results of variance partitioning of a HMSC model produced by `HMSC::computeVariancePartitioning` as a barplot
#' @param hM A fitted HMSC model
#' @param VP A HMSC variance partitioning produced by `HMSC::computeVariancePartitioning`
#' @param ... other parameters passed to the barplot function
#'
#'
#' @examples
#'
#' @export

plotVariancePartitioning=function (hM,VP, ...)
{
   ng = dim(VP$vals)[1]
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
           legend = leg, col = heat.colors(ng, alpha = 1),...)
#   mtext("Species", 1,line = 1)
}
