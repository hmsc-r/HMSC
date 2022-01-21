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
    parts <- unique(partition)
    nfolds <- length(parts)
    if (thin > 1 || start > 1)
        postN <- sum(sapply(hM$postList, function(z)
            length(seq(from=start, to=length(z), by=thin))))
    else
        postN <- Reduce(sum, lapply(hM$postList, length))
    ## output array
    predArray <- array(NA, c(hM$ny, hM$ns, postN))
    ## STEP 1: define new Hmsc model for each nfolds partition
    ## only implement for the simple case first
    if (is.list(hM$X))
        stop("not yet implemented for a list of model matrices")
    if (hM$ncRRR > 0)
        stop("not yet implemendted for RRR models")
    ## no need to parallelism in defining model
    setHmsc <- function(k, hM) {
        Hmsc(Y = hM$Y[k == partition, , drop=FALSE],
             X = hM$X[k == partition,, drop=FALSE],
             distr = hM$distr,
             studyDesign = hM$dfPi[k == partition,, drop=FALSE],
             Tr = hM$Tr, C = hM$C, ranLevels = hM$rL)
    }
    hM1 <- lapply(parts, function(k) setHmsc(k, hM))

    ## STEP 2: sample Hmsc models for each nfolds * nChains case
    ## STEP 3: combine predictions
    ## OUT
    hM1
}
