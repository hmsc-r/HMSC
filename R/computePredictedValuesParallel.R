### Complete re-write of computePredictedValues for parallel processing

#' @importFrom predict predict
#' @importFrom abind abind

#' @export
`pcomputePredictedValues` <-
    function(hM, partition=NULL, partition.sp=NULL, start=1, thin=1,
             Yc=NULL, mcmcStep=1, expected=TRUE, initPar=NULL,
             nParallel=1, clusterType = "fork",
             nChains = length(hM$postList), updater=list(),
             verbose = hM$verbose, alignPost = TRUE)
{
    ## No cross-validation, just return what we got and bail out
    if (is.null(partition)) {
        postList <- poolMcmcChains(hM$postList, start=start, thin=thin)
        pred <- predict(hM, post=postList, Yc=Yc, mcmcStep=1, expected=expected)
        ## Done! Return and exit
        return(abind(pred, along=3))
    }
    ## We have partitions: start with the the simple case without
    ## species partitions and implement first only fork clusters
    if (!missing(clusterType))
        .NotYetUsed("clusterType", error = FALSE)
    if (!missing(partition.sp))
        .NotYetUsed("partition.sp", error = TRUE)
    ## STAGE 1: Basic housekeeping
    nfolds <- length(unique(partition))
    if (thin > 1 || start > 1)
        postN <- sum(sapply(hM$postList, function(z)
            length(seq(from=start, to=length(z), by=thin))))
    else
        postN <- Reduce(sum, lapply(hM$postList, length))
    ## output array
    predArray <- array(NA, c(hM$ny, hM$ns, postN))
    ## STEP 1: define new Hmsc model for each nfolds partition
    ## STEP 2: sample Hmsc models for each nfolds * nChains case
    ## STEP 3: combine predictions
    ## OUT
    predArray
}
