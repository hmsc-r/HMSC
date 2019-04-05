#' @title setPriors
#'
#' @description Sets or resets priors to objects
#'
#' @param x Hmsc or HmscRandolLevel object
#'
#' @return Object of same type as first input
#'
#' @seealso setPriors.Hmsc, setPriors.HmscRandomLevel
#'
#' @export

setPriors <- function(x, ...) {
   UseMethod("setPriors")
}
