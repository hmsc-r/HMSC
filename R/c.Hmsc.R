#' @title c.Hmsc
#'
#' @description
#' Combine Posterior Samples of Several Hmsc Models
#'
#' Functions can be used to add more samples to existing sampled
#' \code{Hmsc} models. New chains can be added with \code{c} and old
#' chains can be continued with \code{merge}. In general, adding new
#' chains is less risky and continuing chains more prone to errors.
#'
#' @details
#' Function \code{c} combines posterior samples of several sampled
#' \code{\link{Hmsc}} models (see \code{\link{sampleMcmc}}) as new
#' chains in the first fitted model. The combined models must be
#' comparable, and there are some tests for detecting non-equal
#' models. These tests will only give warning, and it is at user
#' deliberation to decide which models and which posterior samples can
#' be combined.  You should be careful not start two models from the
#' same random number seed, because these will only duplicate your
#' data instead of providing new independent samples.
#'
#' Function \code{merge} adds posterior samples of second model
#' (\code{y}) to the first one (\code{x}). If the second model was
#' started from the last sampled values of the first one, this is the
#' same as continuing the chain (and other use cases are
#' dubious). Function \code{getLastPar} extracts the last parameters
#' of a sampled model, and its output can be used as \code{initPar} in
#' \code{sampleMcmc} to continue sampling. Users should be very
#' careful that both models are defined and sampled in the same way,
#' except that it makes no sense to have \code{transient} when
#' sampling is started from the last values of the first model. Adding
#' new samples really makes sense only when sampling is started from
#' the end of the first model, and other cases should rather be
#' handled as adding new chains with \code{c}.
#'
#' @param ... Sampled \code{Hmsc} models with posterior samples that
#'     will be added as new chaings in the first listed model, or
#'     ignored in \code{merge}.
#' @return An \code{\link{Hmsc}} model with chains of posterior
#'     samples.
#' @examples
#' ## Fit a toy model with two chains
#' m1 <- sampleMcmc(TD$m, samples=10, transient=5, nChains=2, verbose=0)
#' ## Need more data? Add chains: check carefully that these are
#' ## sampled exactly like the previous model
#' m2 <- sampleMcmc(TD$m, nChains=2, samples=10, transient=5, verbose=0)
#' ## Now four chains
#' m4 <- c(m1, m2)
#' m4
#'
#' @export
`c.Hmsc` <-
    function(...)
{
    ## get models
    hMList <- list(...)
    ## Check inputs
    ## all elements are Hmsc objects?
    if (!all(sapply(hMList, inherits, what = "Hmsc")))
        stop("all elements should be Hmsc objects")
    ## all Hmsc objects are equal. This does not check elements set by
    ## sampleMcmc - sampling features are checked separately
    ## before. This will give false positives and therefore we just
    ## warn. As of now, we do not check elements that are changed in
    ## sampleMcmc (but we have separate checks for some of these):
    ## samples, transient, thin, verbose, adaptNf, randSeed, postList,
    ## HmscVersion; we do not check formulae but only model.matrix
    ## procuded (XFormula, XRRRFormula, TrFormula) as typographically
    ## different formulae can define identical models; we do not check
    ## phyloTree (but only C); call. What about *Names, data.frames
    ## (built to model.matrix)?
    checkItems <-
        c("Y", "XData", "X", "XScaled", "XRRRData", "XRRRScaled",
          "YScaled", "XInterceptInd", "studyDesign", "ranLevels",
          "ranLevelsUsed", "dfPi", "rL", "Pi", "TrData","Tr",
          "TrScaled", "TrInterceptInd", "C", "distr", "ny", "ns",
          "nc", "ncNRRR", "ncRRR", "ncORRR", "ncsel", "nr", "nt",
          "nf", "ncr", "ncs", "np", "spNames", "covNames", "trNames",
          "rLNames", "XScalePar", "XRRRScalePar", "YScalePar",
          "TrScalePar", "V0", "f0", "mGamma", "UGamma", "aSigma",
          "bSigma", "nu", "a1", "b1", "a2", "b2", "rhopw", "nuRRR",
          "a1RRR", "b1RRR", "a2RRR", "b2RRR",  "initPar", "repN")
    tmp <- hMList[[1]]
    objCheck <- lapply(hMList[-1],
                       function(z) all.equal(tmp[checkItems], z[checkItems],
                                             check.attributes=FALSE))
    if(!all(equalObjs <- sapply(objCheck, isTRUE))) {
        warning("some objects differ from the first: may be false alarm, but check")
        pick <- which(!equalObjs)
        cat("Some objects differ when compared to the first object\n")
        cat("object(s) ", paste0(pick+1, collapse=", "), ":\n\n", sep="")
        print(objCheck[pick])
    }
    ## chains should not start from the same random seed
    tmp <- do.call(rbind, lapply(hMList, function(x) x$randSeed))
    if (any(duplicated(tmp)))
        warning("some chains start from the same random seed and may give duplicated samples")
    ## all chains should be same size
    tmp <- unlist(lapply(hMList, function(x) x$samples))
    if (!all(tmp[1] == tmp))
        warning("chains had different lengths: ", paste(tmp, collapse = ", "))
    ## all chains should have same thins
    tmp <- unlist(lapply(hMList, function(x) x$thin))
    if (!all(tmp[1] == tmp))
        warning("chains had different thins: ", paste(tmp, collapse = ", "))
    tmp <- unlist(lapply(hMList, function(x) x$transient))
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


### Add samples from MCMC chain y to MCMC chain x using
### base::merge(). The call is dictated by merge defintion in base

#' @param x,y Hmsc objects: posterior samples of \code{y} are added to
#'     the samples of \code{x}.
#'
#' @rdname c.Hmsc
#' @export
`merge.Hmsc` <-
    function(x, y, ...)
{
    pl1 <- x$postList
    pl2 <- y$postList
    nl <- length(pl1)
    if (nl != length(pl2))
        stop("models have different numbers of chains")
    for(i in seq_len(nl))
        pl1[[i]] <- append(pl1[[i]], pl2[[i]])
    x$postList <- pl1
    x$samples <- x$samples + y$samples
    x
}
### get last posterior sample for use as initPar

#' @param hM Sampled Hmsc object.
#'
#' @rdname c.Hmsc
#' @export
`getLastPar` <-
    function(hM)
{
    if (!inherits(hM, "Hmsc"))
        stop("function can be applied only to sampled Hmsc objects")
    lastone <- hM$samples
    if (is.null(lastone) || is.null(hM$postList))
        stop("object has no posterior samples")
    lapply(hM$postList, function(z) z[[lastone]])
}
