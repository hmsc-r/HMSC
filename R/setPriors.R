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
#' @examples
#' # Set priors for random level so that there is minimum of 2 latent factors and maximum of 3
#' rL1 = HmscRandomLevel(units=TD$studyDesign$plot)
#' rL1 = setPriors(rL1, nfMax=3, nfMin=2)
#'
#' # Set shrinkage parameters for priors of random level
#' rL1 = HmscRandomLevel(units=TD$studyDesign$plot)
#' rL1 = setPriors(rL1, a1=10, a2=10, b1=1, b2=1)
#'
#' @export

setPriors <- function(...) {
   UseMethod("setPriors")
}
