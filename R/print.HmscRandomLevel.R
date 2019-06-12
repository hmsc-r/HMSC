#' @export

print.HmscRandomLevel = function(x, ...){
   if ((x$sDim)<Inf){
      cat(sprintf("Hmsc random level object with %d units. Spatial dimensionality is %d and number of covariates is %d.\n",
               x$N, x$sDim, x$xDim))
   } else
   {
      cat(sprintf("Hmsc random level object with %d units. Spatiality defined through a distance matrix and number of covariates is %d.\n",
                  x$N, x$xDim))
   }
}
