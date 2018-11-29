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

computePredictedValues = function(hM, partition=NULL, start=1, Yc=NULL, mcmcStep=1, expected=TRUE){
   if(is.null(partition)){
      postList=poolMcmcChains(hM$postList, start = start)
      pred = predict(hM, post=postList, Yc=Yc, mcmcStep=1, expected=expected)
      predArray = abind(pred,along=3)
   } else{
      if(length(partition) != hM$ny){
         stop("HMSC.computePredictedValues: partition parameter must be a vector of length ny")
      }
      nfolds = length(unique(partition))
      mpred = matrix(NA,hM$ny,hM$ns)
      for (k in 1:nfolds){
         print(sprintf("Cross-validation, fold %d out of %d", k, nfolds))
         train = (partition!=k)
         val = (partition==k)
         dfPi = as.data.frame(matrix(NA,sum(train),hM$nr))
         colnames(dfPi) = hM$rLNames
         for(r in seq_len(hM$nr)){
            dfPi[,r] = factor(hM$dfPi[train,r])
         }
         hM1 = Hmsc(Y=as.matrix(hM$Y[train,]), X=as.matrix(hM$X[train,]), Xs=as.matrix(hM$Xs[train,]), Xv=as.matrix(hM$Xv[train,]), dist="probit", studyDesign=dfPi, Tr=hM$Tr, C=hM$C, ranLevels=hM$rL)
         hM1$distr=hM$distr
         # HOW TO BETTER SET THE DISTRIBUTION?
         # NEED TO INHERIT PRIORS, SCALINGS ETC. FROM hM
         hM1 = sampleMcmc(hM1, samples=hM$samples, thin=hM$thin, transient=hM$transient,adaptNf=hM$adaptNf, nChains=1)
         postList = hM1$postList[[1]]
         postList = postList[start:length(postList)]
         dfPi = as.data.frame(matrix(NA,sum(val),hM$nr))
         colnames(dfPi) = hM$rLNames
         for (r in seq_len(hM$nr)){
            dfPi[,r] = factor(hM$dfPi[val,r])
         }
         pred1 = predict(hM1, post=postList, X=as.matrix(hM$X[val,,drop=FALSE]), studyDesign=dfPi, Yc=Yc[val,,drop=FALSE], mcmcStep=mcmcStep, expected=expected)
         predArray = abind(pred1,along=3)
      }
   }
   colnames(mpred) = hM$spNames
   return(predArray)
}
