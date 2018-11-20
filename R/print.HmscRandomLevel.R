#' @title print.HmscRandomLevel
#'
#' @description print method for HmscRandomLevel
#' @param rL HmscRandomLevel object
#'
#'
#'
#'
#' @examples
#'

print.HmscRandomLevel = function(rL){
   cat(sprintf("Hmsc %s random level object with %d units\n", if(rL$sDim>0) "spatial" else "non-spatial",rL$N))
}
