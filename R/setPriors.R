#' @title setPriors
#'
#' @description Sets or resets priors to objects
#'
#' @param \dots Hmsc or HmscRandolLevel object and other arguments.
#'
#' @return Object of same type as first input
#'
#' @seealso setPriors.Hmsc, setPriors.HmscRandomLevel
#'
#' @export

setPriors <- function(...) {
   UseMethod("setPriors")
}
