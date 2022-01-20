#' @title computePredictedValues
#'
#' @description Computes predicted values from the fitted \code{Hmsc} model
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param partition partition vector for cross-validation created by
#'     \code{\link{createPartition}}
#' @param partition.sp species partitioning vector for conditional
#'     cross-validation
#' @param start index of first MCMC sample included
#' @param thin thinning interval of posterior distribution
#' @param Yc response matrix on which the predictions are to be
#'     conditioned
#' @param mcmcStep number of MCMC steps used to make conditional
#'     predictions
#' @param expected whether expected values (TRUE) or realizations
#'     (FALSE) are to be predicted
#' @param initPar a named list of parameter values used for
#'     initialization of MCMC states
#' @param nParallel number of parallel processes by which the chains
#'     are executed, or alternatively a predefined socket cluster
#'     (passed to \code{\link{sampleMcmcm}})
#' @param clusterType cluster type in parallel processing; socket
#'     clusters are also used if supplied in \code{nParallel} and
#'     always in Windows (passed to \code{\link{sampleMcmc}}
#' @param nChains number of independent MCMC chains to be run
#' @param updater a named list, specifying which conditional updaters
#'     should be omitted
#' @param verbose the interval between MCMC steps printed to the
#'     console
#' @param alignPost boolean flag indicating whether the posterior of
#'     each chains should be aligned
#'
#' @return an array of model predictions, made for each posterior sample
#'
#' @details If the option \code{partition} is not used, the posterior predictive distribution is based on the model
#' fitted to the full data. If the option \code{partition} is used but \code{partition.sp} is not used, the posterior predictive distribution
#' is based on cross-validation over the sampling units. If \code{partition.sp} is additionally used, then, when predictions are made for
#' each fold of the sampling units, the predictions are done separately for each fold of species. When making the predictions
#' for one fold of species, the predictions are conditional on known occurrences (those observed in the data) of the species
#' belonging to the other folds. If \code{partition.sp} is used, the parameter \code{mcmcStep} should be set high enough to obtain
#' appropriate conditional predictions. The option \code{Yc} can be used alternatively to \code{partition.sp} if the conditioning is to be done
#' based on a fixed set of data (independently of which sampling unit and species the predictions are made for).
#'
#'
#' @seealso \code{\link{predict.Hmsc}}
#'
#' @examples
#' # Compute predicted values using a previously fitted HMSC model
#' preds = computePredictedValues(TD$m)
#'
#' \dontrun{
#' # Compute predicted values for a previously fitted HMSC model using 2 folds
#' partition = createPartition(TD$m, nfolds = 2)
#' predsCV1 = computePredictedValues(TD$m,partition=partition)
#'
#' # Compute conditional predictions for a previously fitted HMSC model using 2 folds
#' partition = createPartition(TD$m, nfolds = 2)
#' predsCV2 = computePredictedValues(TD$m, partition = partition,
#' partition.sp = 1:TD$m$ns, mcmcStep = 100)
#' }
#'
#' @importFrom stats predict
#' @importFrom abind abind
#' @importFrom rlang duplicate
#' @export

computePredictedValues =
    function(hM, partition=NULL, partition.sp=NULL, start=1, thin=1,
             Yc=NULL, mcmcStep=1, expected=TRUE, initPar=NULL,
             nParallel=1, clusterType = c("socket","fork"),
             nChains = length(hM$postList), updater=list(),
             verbose = hM$verbose, alignPost = TRUE)
{
    clusterType <- match.arg(clusterType)
    if(is.null(partition)){
        postList = poolMcmcChains(hM$postList, start=start, thin=thin)
        pred = predict(hM, post=postList, Yc=Yc, mcmcStep=1, expected=expected)
        predArray = abind(pred, along=3)
    } else {
        if(length(partition) != hM$ny){
            stop("partition parameter must be a vector of length ny")
        }
        nfolds = length(unique(partition))
        if (thin > 1 || start > 1)
            postN = sum(sapply(hM$postList, function(z)
                length(seq(from=start, to=length(z), by=thin))))
        else
            postN = Reduce(sum, lapply(hM$postList, length))
        predArray = array(NA,c(hM$ny,hM$ns,postN))
        ## collect models for each fold
        hM1 <- list()
        for (k in 1:nfolds){
            cat(sprintf("Cross-validation, fold %d out of %d\n", k, nfolds))
            train = (partition!=k)
            val = (partition==k)
            ## stringsAsFactors probably not needed below
            dfPi = as.data.frame(matrix(NA,sum(train),hM$nr),
                                 stringsAsFactors = TRUE)
            colnames(dfPi) = hM$rLNames
            for(r in seq_len(hM$nr)){
                dfPi[,r] = factor(hM$dfPi[train,r])
            }
            switch(class(hM$X)[1L],
                   matrix = {
                XTrain = hM$X[train,,drop=FALSE]
                XVal = hM$X[val,,drop=FALSE]
            },
            list = {
                XTrain = lapply(hM$X, function(a) a[train,,drop=FALSE])
                XVal = lapply(hM$X, function(a) a[val,,drop=FALSE])
            }
            )
            if(hM$ncRRR>0){
                XRRRTrain = hM$XRRR[train,,drop=FALSE]
                XRRRVal = hM$XRRR[val,,drop=FALSE]
            } else {
                XRRRTrain=NULL
                XRRRVal=NULL
            }
            ## NB, if ranLevels was NULL, rL (and ranLevels) do not exist
            ## in hM (Hmsc() removed them), and partial matching gives
            ## the value (character(0)) of hM$rLNames to ranLevels
            ## below. If this disturbs, it should probably be fixed in
            ## Hmsc(), but currently ranLevels=character(0) seems to
            ## work.
            hM1[[k]] = Hmsc(Y=hM$Y[train,,drop=FALSE], X=XTrain,
                            XRRR=XRRRTrain, ncRRR = hM$ncRRR,
                            XSelect = hM$XSelect,
                            distr=hM$distr, studyDesign=dfPi, Tr=hM$Tr, C=hM$C,
                            ranLevels=hM$rL)
            setPriors(hM1[[k]], V0=hM$V0, f0=hM$f0, mGamma=hM$mGamma,
                      UGamma=hM$UGamma, aSigma=hM$aSigma, bSigma=hM$bSigma,
                      nu=hM$nu, a1=hM$a1, b1=hM$b1, a2=hM$a2, b2=hM$b2,
                      rhopw=hM$rhowp)
            hM1[[k]]$YScalePar = hM$YScalePar
            hM1[[k]]$YScaled = (hM1[[k]]$Y -
                                matrix(hM1[[k]]$YScalePar[1,], hM1[[k]]$ny,
                                       hM1[[k]]$ns, byrow=TRUE)) /
                matrix(hM1[[k]]$YScalePar[2,], hM1[[k]]$ny, hM1[[k]]$ns,
                       byrow=TRUE)
            hM1[[k]]$XInterceptInd = hM$XInterceptInd
            hM1[[k]]$XScalePar = hM$XScalePar
            switch(class(hM$X)[1L],
                   matrix = {
                hM1[[k]]$XScaled = (hM1[[k]]$X -
                                    matrix(hM1[[k]]$XScalePar[1,],
                                           hM1[[k]]$ny, hM1[[k]]$ncNRRR,
                                           byrow=TRUE)) /
                    matrix(hM1[[k]]$XScalePar[2,], hM1[[k]]$ny, hM1[[k]]$ncNRRR,
                           byrow=TRUE)
            },
            list = {
                hM1[[k]]$XScaled = list()
                for (zz in seq_len(length(hM1[[k]]$X))) {
                    hM1[[k]]$XScaled[[zz]] = (hM1[[k]]$X[[zz]] -
                                              matrix(hM1[[k]]$XScalePar[1,],
                                                     hM1[[k]]$ny,
                                                     hM1[[k]]$ncNRRR,
                                                     byrow=TRUE)) /
                        matrix(hM1[[k]]$XScalePar[2,], hM1[[k]]$ny,
                               hM1[[k]]$ncNRRR,
                               byrow=TRUE)
                }
            }
            )
            if(hM1[[k]]$ncRRR>0) {
                hM1[[k]]$XRRRScalePar = hM$XRRRScalePar
                hM1[[k]]$XRRRScaled = (hM1[[k]]$XRRR -
                                       matrix(hM1[[k]]$XRRRScalePar[1,],
                                              hM1[[k]]$ny, hM1[[k]]$ncORRR,
                                              byrow=TRUE)) /
                    matrix(hM1[[k]]$XRRRScalePar[2,], hM1[[k]]$ny,
                           hM1[[k]]$ncORRR, byrow=TRUE)
            }
            hM1[[k]]$TrInterceptInd = hM$TrInterceptInd
            hM1[[k]]$TrScalePar = hM$TrScalePar
            hM1[[k]]$TrScaled = (hM1[[k]]$Tr -
                                 matrix(hM1[[k]]$TrScalePar[1,], hM1[[k]]$ns,
                                        hM1[[k]]$nt, byrow=TRUE)) /
                matrix(hM1[[k]]$TrScalePar[2,], hM1[[k]]$ns, hM1[[k]]$nt,
                       byrow=TRUE)
        }
        for (k in 1:nfolds) {
            val = partition == k
            hM1[[k]] = sampleMcmc(hM1[[k]], samples=hM$samples, thin=hM$thin,
                                  transient=hM$transient, adaptNf=hM$adaptNf,
                                  initPar=initPar, nChains=nChains,
                                  nParallel=nParallel,
                                  clusterType = clusterType, updater = updater,
                                  verbose = verbose, alignPost=alignPost)
            postList = poolMcmcChains(hM1[[k]]$postList, start=start,
                                      thin = thin)
            ## stringsAsFactors probably not needed below
            dfPi = as.data.frame(matrix(NA, sum(val), hM$nr),
                                 stringsAsFactors = TRUE)
            colnames(dfPi) = hM$rLNames
            for(r in seq_len(hM$nr)) {
                dfPi[,r] = factor(hM$dfPi[val,r])
            }
            if(is.null(partition.sp)) {
                pred1 = predict(hM1[[k]], post=postList, X=XVal, XRRR=XRRRVal,
                                studyDesign=dfPi, Yc=Yc[val,,drop=FALSE],
                                mcmcStep=mcmcStep, expected=expected)
                pred1Array = abind(pred1,along=3)
            } else {
                hM2 = duplicate(hM)
                if(is.null(hM2$rLPar)) {
                    hM2$rLPar = computeDataParameters(hM2)$rLPar
                }
                postList2 = duplicate(postList)
                for(r in seq_len(hM$nr)){
                    postEta = lapply(postList, function(c) c$Eta[[r]])
                    postAlpha = lapply(postList, function(c) c$Alpha[[r]])
                    predPostEta =
                        predictLatentFactor(unitsPred=levels(hM2$dfPi[,r]),
                                            units=levels(hM1$dfPi[,r]),
                                            postEta=postEta,
                                            postAlpha=postAlpha, rL=hM$rL[[r]])
                    for(i in seq_len(length(postList))) {
                        postList2[[i]]$Eta[[r]] = predPostEta[[i]]
                    }
                }

                pred1Array = array(dim=c(sum(val), hM$ns, postN))
                nfolds.sp = length(unique(partition.sp))
                for (i in 1:nfolds.sp){
                    if(verbose)
                        cat("   ... species partition", i, "\n")
                    train.sp = (partition.sp!=i)
                    val.sp = (partition.sp==i)
                    YcFull = hM$Y
                    YcFull[val,val.sp] = NA
                    pred2 = predict(hM2, post=postList2, X=hM2$X, XRRR=hM2$XRRR,
                                    studyDesign=hM2$studyDesign, Yc=YcFull,
                                    mcmcStep=mcmcStep, expected=expected)
                    pred2Array = abind(pred2,along=3)
                    pred1Array[,val.sp,] = pred2Array[val, val.sp,]
                }
            }
            predArray[val,,] = pred1Array
        }
    }
    return(predArray)
}
