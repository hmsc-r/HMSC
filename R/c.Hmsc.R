#' Combine Posterior Samples of Several Hmsc Models
#'
#' Function combines posterior samples of several sampled
#' \code{\link{Hmsc}} models (see \code{\link{sampleMcmc}}) as new
#' chains in the first fitted model. The combined models must be
#' comparable, and there are some tests for detecting non-equal
#' models. These tests will only give warning, and it is at user
#' deliberation to decide which models and which posterior samples can
#' be combined.  You should be careful not start two models from the
#' same random number seed, because these will only duplicate your
#' data instead of providing new independent samples.
#'
#' @param ... Sampled \code{Hmsc} models with posterior samples that
#'     will be added as new chaings in the first listed model.
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
