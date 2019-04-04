#' @title print.HmscRandomLevel
#'
#' @description print method for HmscRandomLevel
#' @param rL HmscRandomLevel object
#'
#' @examples
#'
#' @export

print.HmscRandomLevel = function(rL){
   cat(sprintf("Hmsc random level object with %d units. Spatial dimentionality is %d and number of covariates is %d.\n",
               rL$N,rL$sDim,rL$xDim))
}
