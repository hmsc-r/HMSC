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
    ## Check inputs:
    ## all elements are Hmsc objects?
    if (!all(sapply(hMList, inherits, what = "Hmsc")))
        stop("all elements should be Hmsc objects")
    ## all Hmsc object have identical Calls (can give false alarms)
    tmp <- sapply(hMList, getCall)
    if (!all(sapply(tmp, identical, y = tmp[[1]])))
        warning("not all elements have identical function Call")
    ## all chains should be same size
    tmp <- sapply(hMList, function(x) x$sample)
    if (!all(tmp[1] == tmp))
        warning("chains had different lengths: ", paste(tmp, collapse = ", "))
    ## all chains should have same thins
    tmp <- sapply(hMList, function(x) x$thin)
    if (!all(tmp[1] == tmp))
        warning("chains had different thins: ", paste(tmp, collapse = ", "))
    tmp <- sapply(hMList, function(x) x$transient)
    if (!all(tmp[1] == tmp))
        warning("chains had different transients: ", paste(tmp, collapse = ", "))
    ## extract postLists
    pLists <- lapply(hMList, function(x) x$postList)
    ## combine postLists
    pLists <- do.call(c, pLists)
    ## replace postList of the first model with the combined list
    hMList[[1]]$postList <- pLists
    hMList[[1]]
}
