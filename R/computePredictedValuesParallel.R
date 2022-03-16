### Complete re-write of computePredictedValues for parallel processing

#' @importFrom stats predict
#' @importFrom abind abind
#' @importFrom parallel mclapply detectCores

#' @rdname computePredictedValues
#' @export
`pcomputePredictedValues` <-
    function(hM, partition=NULL, partition.sp=NULL, start=1, thin=1,
             Yc=NULL, mcmcStep=1, expected=TRUE, initPar=NULL,
             nParallel=1, clusterType = "fork",
             nChains = length(hM$postList), updater=list(),
             verbose = nParallel == 1, alignPost = TRUE)
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
    clusterType <- match.arg(clusterType)
    if (nParallel > 1 && .Platform$OS.type == "windows")
        stop("parallel processing not yet implemented for Windows: use Linux or Mac")
    if (!missing(partition.sp))
        .NotYetUsed("partition.sp", error = TRUE)
    ## STAGE 1: Basic housekeeping
    parts <- sort(unique(partition))
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
    ## Pack setting training model into one function
    setHmsc <- function(k, hM) {
        train <- k != partition
        m <- Hmsc(Y = hM$Y[train, , drop=FALSE],
                  X = hM$X[train, , drop=FALSE],
                  XRRR = hM$XRRR[train, , drop=FALSE],
                  ncRRR = hM$ncRRR, XSelect = hM$XSelect,
                  distr = hM$distr,
                  studyDesign = droplevels(hM$dfPi[train,, drop=FALSE]),
                  Tr = hM$Tr, C = hM$C, ranLevels = hM$rL)
        ## old code calls here setPriors, but that does nothing as its
        ## result is not saved, and it is currently skipped: CHECK
        ## THIS!
        m$YScalePar <- hM$YScalePar
        m$YScaled <- scale(m$Y, m$YScalePar[1,], m$YScalePar[2,])
        m$XInterceptInd <- hM$XInterceptInd
        m$XScalePar <- hM$XScalePar
        m$XScaled <- scale(m$X, m$XScalePar[1,], m$XScalePar[2,])
        m$TrInterceptInd <- hM$TrInterceptInd
        m$TrScalePar <- hM$TrScalePar
        m$TrScaled <- scale(m$Tr, m$TrScalePar[1,], m$TrScalePar[2,])
        m
    }
    ## parallelism would only slow-down (due to start-off time)
    hM1 <- lapply(parts, function(k) setHmsc(k, hM))
    ## STEP 2: sample Hmsc models for each nfolds * nChains case
    chains <- length(hM$postList)
    threads <- nfolds * chains
    idfold <- rep(parts, each = chains)
    seeds <- sample.int(.Machine$integer.max, threads)
    ## to be called in parallel for each chain x fold, and
    ## therefore we set nChains=1, nParallel=1 within
    getSample <- function(i, ...) {
        set.seed(seeds[i])
        k <- idfold[i]
        message("starting thread ", i, "/", threads)
        m <- sampleMcmc(hM1[[k]], samples = hM$samples, thin = hM$thin,
                        transient = hM$transient, adaptNf = hM$adaptNf,
                        initPar = initPar, nChains = 1, nParallel = 1,
                        updater = updater, verbose = verbose,
                        alignPost = alignPost)
        message("finished thread ", i, "/", threads)
        attr(m, "fold") <- k
        m
    }
    ## the next call can be made parallel
    if (nParallel <= 1) # serial
        mods <- lapply(seq_len(threads), function(i, hM, hM1)
            getSample(i, hM, hM1))
    else { # parallel
        Ncores <- min(nParallel, threads, detectCores())
        message("using ", Ncores, " cores")
        mods <- mclapply(seq_len(threads), function(i, hM, hM1)
            getSample(i, hM, hM1),
            mc.cores = Ncores, mc.preschedule = FALSE)
    }
    ## STEP 3: combine predictions: this is still a loop
    idfold <- sapply(mods, attr, which = "fold")
    for (p in parts) {
        val <- partition == p
        m <- do.call(c.Hmsc, mods[which(idfold == p)])
        m <- alignPosterior(m)
        postList <- poolMcmcChains(m$postList, start=start, thin=thin)
        dfPi <- droplevels(hM$dfPi[val,, drop=FALSE])
        pred1 <- predict(m, post=postList, X = hM$X[val,, drop=FALSE],
                         Yc = Yc[val,, drop=FALSE], studyDesign = dfPi,
                         mcmcStep = mcmcStep, expected = expected)
        predArray[val,,] <- abind(pred1, along=3)
    }
    ## OUT
    predArray
}
