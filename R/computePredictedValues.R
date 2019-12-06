#' @title computePredictedValues
#'
#' @description Computes predicted values from the fitted \code{Hmsc} model
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param partition partition vector for cross-validation created by \code{\link{createPartition}}
#' @param partition.sp species partitioning vector for conditional cross-validation
#' @param start index of first MCMC sample included
#' @param thin thinning interval of posterior distribution
#' @param Yc response matrix on which the predictions are to be conditioned
#' @param mcmcStep number of MCMC steps used to make conditional predictions
#' @param expected whether expected values (TRUE) or realizations (FALSE) are to be predicted
#' @param initPar a named list of parameter values used for initialization of MCMC states
#' @param nParallel number of parallel processes by which the chains are executed
#' @param nChains number of independent MCMC chains to be run
#' @param updater a named list, specifying which conditional updaters should be ommitted
#' @param verbose the interval between MCMC steps printed to the console
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
#' partition.sp = 1:m$ns, mcmcStep = 100)
#' }
#'
#' @importFrom stats predict
#' @importFrom abind abind
#' @export

computePredictedValues = function(hM, partition=NULL, partition.sp=NULL, start=1, thin=1,
                                  Yc=NULL, mcmcStep=1, expected=TRUE, initPar=NULL,
                                  nParallel=1, nChains = length(hM$postList), updater=list(), verbose = hM$verbose){
   if(is.null(partition)){
      postList = poolMcmcChains(hM$postList, start=start, thin=thin)
      pred = predict(hM, post=postList, Yc=Yc, mcmcStep=1, expected=expected)
      predArray = abind(pred, along=3)
   } else{
      if(length(partition) != hM$ny){
         stop("HMSC.computePredictedValues: partition parameter must be a vector of length ny")
      }
      nfolds = length(unique(partition))
      postN = Reduce(sum, lapply(hM$postList, length))
      predArray = array(NA,c(hM$ny,hM$ns,postN))
      for (k in 1:nfolds){
         print(sprintf("Cross-validation, fold %d out of %d", k, nfolds))
         train = (partition!=k)
         val = (partition==k)
         dfPi = as.data.frame(matrix(NA,sum(train),hM$nr))
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
         hM1 = Hmsc(Y=hM$Y[train,,drop=FALSE], X=XTrain, XRRR=XRRRTrain, ncRRR = hM$ncRRR, XSelect = hM$XSelect, distr=hM$distr, studyDesign=dfPi, Tr=hM$Tr, C=hM$C, ranLevels=hM$rL)
         setPriors(hM1, V0=hM$V0, f0=hM$f0, mGamma=hM$mGamma, UGamma=hM$UGamma, aSigma=hM$aSigma, bSigma=hM$bSigma,
                   nu=hM$nu, a1=hM$a1, b1=hM$b1, a2=hM$a2, b2=hM$b2, rhopw=hM$rhowp)
         hM1$YScalePar = hM$YScalePar
         hM1$YScaled = (hM1$Y - matrix(hM1$YScalePar[1,],hM1$ny,hM1$ns,byrow=TRUE)) / matrix(hM1$YScalePar[2,],hM1$ny,hM1$ns,byrow=TRUE)
         hM1$XInterceptInd = hM$XInterceptInd
         hM1$XScalePar = hM$XScalePar
         switch(class(hM$X)[1L],
                matrix = {
                   hM1$XScaled = (hM1$X - matrix(hM1$XScalePar[1,],hM1$ny,hM1$ncNRRR,byrow=TRUE)) / matrix(hM1$XScalePar[2,],hM1$ny,hM1$ncNRRR,byrow=TRUE)
                },
                list = {
                   hM1$XScaled = list()
                   for (zz in seq_len(length(hM1$X))){
                      hM1$XScaled[[zz]] = (hM1$X[[zz]] - matrix(hM1$XScalePar[1,],hM1$ny,hM1$ncNRRR,byrow=TRUE)) / matrix(hM1$XScalePar[2,],hM1$ny,hM1$ncNRRR,byrow=TRUE)
                   }
                }
         )
         if(hM1$ncRRR>0){
            hM1$XRRRScalePar = hM$XRRRScalePar
            hM1$XRRRScaled = (hM1$XRRR - matrix(hM1$XRRRScalePar[1,],hM1$ny,hM1$ncORRR,byrow=TRUE)) / matrix(hM1$XRRRScalePar[2,],hM1$ny,hM1$ncORRR,byrow=TRUE)
         }
         hM1$TrInterceptInd = hM$TrInterceptInd
         hM1$TrScalePar = hM$TrScalePar
         hM1$TrScaled = (hM1$Tr - matrix(hM1$TrScalePar[1,],hM1$ns,hM1$nt,byrow=TRUE)) / matrix(hM1$TrScalePar[2,],hM1$ns,hM1$nt,byrow=TRUE)
         hM1 = sampleMcmc(hM1, samples=hM$samples, thin=hM$thin, transient=hM$transient, adaptNf=hM$adaptNf,
                          initPar=initPar, nChains=nChains, nParallel=nParallel, updater = updater, verbose = verbose)
         postList = poolMcmcChains(hM1$postList, start=start)
         dfPi = as.data.frame(matrix(NA,sum(val),hM$nr))
         colnames(dfPi) = hM$rLNames
         for (r in seq_len(hM$nr)){
            dfPi[,r] = factor(hM$dfPi[val,r])
         }
         if(is.null(partition.sp)){
            pred1 = predict(hM1, post=postList, X=XVal, XRRR=XRRRVal, studyDesign=dfPi, Yc=Yc[val,,drop=FALSE], mcmcStep=mcmcStep, expected=expected)
            pred1Array = abind(pred1,along=3)
         } else {
            pred1Array =  array(dim=c(sum(val),hM$ns,postN))
            nfolds.sp = length(unique(partition.sp))
            for (i in 1:nfolds.sp){
               train.sp = (partition.sp!=i)
               val.sp = (partition.sp==i)
               Yc = matrix(NA,nrow=sum(val), ncol=hM$ns)
               Yc[,train.sp] = hM$Y[val,train.sp,drop=FALSE]
               pred2 = predict(hM1, post=postList, X=XVal, XRRR=XRRRVal,studyDesign=dfPi, Yc=Yc, mcmcStep=mcmcStep, expected=expected)
               pred2Array = abind(pred2,along=3)
               pred1Array[,val.sp,] = pred2Array[,val.sp,]
            }
         }
         predArray[val,,] = pred1Array
      }
   }
   return(predArray)
}
