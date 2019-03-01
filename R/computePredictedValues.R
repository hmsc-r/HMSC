#' @title computePredictedValues
#'
#' @description Computes predicted values from the trained model
#' @param partition vector of observations partitioning for cross-validation
#' @param start index of first MCMC sample included
#'
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export

computePredictedValues = function(hM, partition=NULL, start=1, Yc=NULL, mcmcStep=1, expected=TRUE, initPar=NULL, nCores=1){
   if(is.null(partition)){
      postList = poolMcmcChains(hM$postList, start=start)
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
         switch(class(hM$X),
            matrix = {
               XTrain = hM$X[train,,drop=FALSE]
               XVal = hM$X[val,,drop=FALSE]
            },
            list = {
               XTrain = lapply(hM$X, function(a) a[train,,drop=FALSE])
               XVal = lapply(hM$X, function(a) a[val,,drop=FALSE])
            }
         )
         hM1 = Hmsc(Y=hM$Y[train,,drop=FALSE], X=XTrain, dist="probit", studyDesign=dfPi, Tr=hM$Tr, C=hM$C, ranLevels=hM$rL)
         hM1$distr = hM$distr
         # HOW TO BETTER SET THE DISTRIBUTION?
         # NEED TO INHERIT PRIORS, SCALINGS ETC. FROM hM
         hM1 = sampleMcmc(hM1, samples=hM$samples, thin=hM$thin, transient=hM$transient, adaptNf=hM$adaptNf, initPar=initPar, nChains=length(hM$postList), nCores=nCores)
         postList = poolMcmcChains(hM1$postList, start=start)
         dfPi = as.data.frame(matrix(NA,sum(val),hM$nr))
         colnames(dfPi) = hM$rLNames
         for (r in seq_len(hM$nr)){
            dfPi[,r] = factor(hM$dfPi[val,r])
         }
         pred1 = predict(hM1, post=postList, X=XVal, studyDesign=dfPi, Yc=Yc[val,,drop=FALSE], mcmcStep=mcmcStep, expected=expected)
         predArray[val,,] = abind(pred1,along=3)
      }
   }
   return(predArray)
}
