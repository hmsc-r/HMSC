#' @title print.HmscRandomLevel
#'
#' @description print method for HmscRandomLevel
#' @param rL HmscRandomLevel object
#'
#' @examples
#'
#' @export

print.HmscRandomLevel = function(rL){
   if ((rL$sDim)<Inf){
      cat(sprintf("Hmsc random level object with %d units. Spatial dimensionality is %d and number of covariates is %d.\n",
               rL$N,rL$sDim,rL$xDim))
   } else
   {
      cat(sprintf("Hmsc random level object with %d units. Spatiality defined through a distance matrix and number of covariates is %d.\n",
                  rL$N,rL$xDim))
   }
}
