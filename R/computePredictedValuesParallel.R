### Complete re-write of computePredictedValues for parallel processing

#' @importFrom stats predict
#' @importFrom parallel mclapply detectCores makeCluster clusterExport
#'     clusterEvalQ clusterApplyLB stopCluster

## Most parameters are the same as in computePredictedValues: document
## only the new one
#' @param useSocket (logical) use socket clusters in parallel processing;
#'     these can be used in all operating systems, but they are
#'     usually slower than forking which can only be used
#'     in non-Windows operating systems (macOS, Linux, unix-like
#'     systems).

#' @rdname computePredictedValues
#' @export
`pcomputePredictedValues` <-
    function(hM, partition=NULL, partition.sp=NULL, start=1, thin=1,
             Yc=NULL, mcmcStep=1, expected=TRUE, initPar=NULL,
             nParallel=1, useSocket = .Platform$OS.type == "windows",
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
    if (nParallel > 1) {
        if (.Platform$OS.type == "windows" && !useSocket) {
            useSocket <- TRUE
            message("setting useSocket=TRUE: the only choice in Windows")
        }
    }
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
        ## X can be a list of matrices
        XTrain <- if(is.matrix(hM$X))
                      hM$X[train, , drop=FALSE]
                  else if(is.list(hM$X))
                      lapply(hM$X, function(a) a[train, , drop=FALSE])
        m <- Hmsc(Y = hM$Y[train, , drop=FALSE], X = XTrain,
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
    if (nParallel < 2) # serial
        mods <- lapply(seq_len(threads), function(i, hM, hM1)
            getSample(i, hM, hM1))
    else { # parallel
        Ncores <- min(nParallel, threads, detectCores())
        message("using ", Ncores, " cores")
        if (!useSocket) { # everywhere except Windows
        mods <- mclapply(seq_len(threads), function(i, hM, hM1)
            getSample(i, hM, hM1),
            mc.cores = Ncores, mc.preschedule = FALSE)
        } else { # socket cluster, works everywhere, incl Windows
            cl <- makeCluster(Ncores)
            clusterExport(cl, "getSample", envir = environment())
            clusterEvalQ(cl, {
                library(BayesLogit);
                library(MCMCpack);
                library(truncnorm);
                library(Matrix);
                library(abind);
                library(Hmsc)})
            mods <- clusterApplyLB(cl, seq_len(threads),
                                   function(i, hM, hM1) getSample(i, hM, hM1))
            stopCluster(cl)
        }
    }
    ## STEP 3: combine predictions: this is still a loop
    idfold <- sapply(mods, attr, which = "fold")
    for (p in parts) {
        message("predictions for partition ", p)
        val <- partition == p
        m <- do.call(c.Hmsc, mods[which(idfold == p)])
        m <- alignPosterior(m)
        postList <- poolMcmcChains(m$postList, start=start, thin=thin)
        dfPi <- droplevels(hM$dfPi[val,, drop=FALSE])
        Xval <- if (is.matrix(hM$X)) hM$X[val, , drop=FALSE]
                else lapply(hM$X, function(a) a[val, , drop=FALSE])
        pred1 <- if (is.null(partition.sp)) {
                     predict(m, post=postList, X = Xval,
                             XRRR = hM$XRRR[val,, drop=FALSE],
                             Yc = Yc[val,, drop=FALSE], studyDesign = dfPi,
                             mcmcStep = mcmcStep, expected = expected)
                 } else {
                     getSpeciesFoldPrediction(hM, m, val, postList, dfPi,
                                              partition.sp = partition.sp,
                                              mcmcStep = mcmcStep,
                                              expected = expected,
                                              nParallel = nParallel,
                                              useSocket = useSocket)
                 }
        predArray[val,,] <- simplify2array(pred1)
    }
    ## OUT
    predArray
}

### Non-exported function to be called in STEP 3 of CV cycle if
### partition.sp was defined. This not parallelized (yet?). NB, it may
### be better to parallelize predict.Hmsc by posterior samples, since
### that is the function that takes time.
#
# @param hM Original non-cv Hmsc model
# @param val Units in this CV partition (logical)
# @param postList Current partition pooled postList
# @param dfPi Current partition random level data frame
# @param partition.sp Partitioning vector for species
# @param mcmcStep Parameter passed to predict
# @param expected Parameter passed to predict
#
#  @return Predictions for current partition
##
`getSpeciesFoldPrediction` <-
    function(hM, hM1, val, postList, dfPi, partition.sp = NULL, mcmcStep,
             expected = expected, nParallel = nParallel,
             useSocket = useSocket)
{
    if (is.null(partition.sp))
        return(NULL)
    nfolds.sp <- length(unique(partition.sp))
    Ncores <- min(nParallel, detectCores())
    message("species prediction using ", Ncores, " cores")
    ## update Hmsc object
    if (is.null(hM$rLPar)) {
        hM$rLPar <- computeDataParameters(hM)$rLPar
    }
    for (r in seq_len(hM$nr)) {
        postEta <- lapply(postList, function(c) c$Eta[[r]])
        postAlpha <- lapply(postList, function(c) c$Alpha[[r]])
        predPostEta <-
            predictLatentFactor(unitsPred = levels(hM$dfPi[,r]),
                                units = levels(hM1$dfPi[,r]),
                                postEta = postEta, postAlpha = postAlpha,
                                rL = hM$rL[[r]])
        for(i in seq_len(length(postList))) {
            postList[[i]]$Eta[[r]] <- predPostEta[[i]]
        }
    }
    ## Handle each species partition
    predArray <- array(dim = c(sum(val), hM$ns, length(postList)))
    for (i in seq_len(nfolds.sp)) {
        val.sp <- partition.sp == i
        YcFull <- hM$Y
        YcFull[val, val.sp] <- NA
        pred <- predict(hM, post = postList, X = hM$X, XRRR = hM$XRRR,
                        studyDesign = hM$studyDesign, Yc = YcFull,
                        mcmcStep = mcmcStep, expected = expected,
                        nParallel = Ncores, useSocket = useSocket)
        pred <- simplify2array(pred)
        predArray[, val.sp, ] <- pred[val, val.sp, ]
        message("finished species fold ", i, "/", nfolds.sp)
    }
    predArray
}
