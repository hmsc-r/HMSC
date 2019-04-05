#' @title print.Hmsc
#'
#' @description print method for Hmsc
#' @param hM Hmsc object
#'
#'
#' @examples
#'
#' @export

print.Hmsc = function(hM){
   cat(sprintf("Hmsc object with %d units, %d species, %d covariates, %d traits and %d random levels\n", hM$ny, hM$ns, hM$nc, hM$nt, hM$nr))
}
