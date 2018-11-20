#' @title Hmsc$plotMarginalEffects
#'
#' @description Plots ...
#'
#' @param pred
#' @param covariate
#' @param measure
#' @param index
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

plotMarginalEffects = function(hM, model, pred, covariate, measure, index=1){

   if (measure == "S"){
      val = t(pred[[covariate]]$S)
      ylab = "Species richness"
   }

   if (measure == "Y"){
      val = t(pred[[covariate]]$Y[,,index])
      ylab = hM$spNames[[index]]
   }
   if (measure == "T"){
      val = t(pred[[covariate]]$Tr[,,index])
      ylab = hM$trNames[[index]]
   }

   xlab = colnames(val)[[1]]
   plot(val[,1], val[,3], ylim = c(min(val[,2]),max(val[,4])), type = "l", xlab = xlab, ylab = ylab)
   polygon(c(val[,1],rev(val[,1])),c(val[,2],rev(val[,4])),col = "grey75", border = FALSE)
   lines(val[,1], val[,3], lwd = 2)
   lines(val[,1], val[,2], col="red",lty=2)
   lines(val[,1], val[,4], col="red",lty=2)
}
