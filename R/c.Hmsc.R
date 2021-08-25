#' Combine Posterior Samples of Several Hmsc Models
#'
#' Function combines posterior samples of several sampled
#' \code{\link{Hmsc}} models (see \code{\link{sampleMcmc}}) as new
#' chains in the first fitted model. The combined models must be
#' identical, but this is not checked in the function and is left as
#' user responsibility.
#'
#' @param ... Sampled \code{Hmsc} models with posterior samples that
#'     will be added as new chaings in the first listed model.
#' @return An \code{\link{Hmsc}} model with chains of posterior
#'     samples.
#' @examples
#' ## Fitted model with two chains
#' TD$m
#' ## New chains: check carefully that these are sampled exactly like
#' ## the previous model
#' m2 <- sampleMcmc(TD$m, nChains=2, samples=100, thin=1, transient=50,
#'     verbose=0)
#' ## Now four chains
#' mnew <- c(TD$m, m2)
#' mnew
#'
#' @export
`c.Hmsc` <-
    function(...)
{
    ## get models
    hMList <- list(...)
    ## extract postLists
    pLists <- lapply(hMList, function(x) x$postList)
    ## combine postLists
    pLists <- do.call(c, pLists)
    ## replace postList of the first model with the combined list
    hMList[[1]]$postList <- pLists
    hMList[[1]]
}
