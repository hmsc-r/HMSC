#' @title computePredictedValues
#'
#' @description Computes predicted values from the trained model
#' @param nfolds number of cross-validation folds
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

computePredictedValues = function(hM, nfolds=NULL, column=NULL, partition=NULL, start=1, Yc=NULL, mcmcStep=1, expected=TRUE){
   if(identical(nfolds,1) || (is.null(nfolds) && is.null(partition))){
      if(!is.null(column)){
         stop("HMSC.computePredictedValues: no column can be specified when calculating with single fold")
      }
      if(!is.null(partition)){
         stop("HMSC.computePredictedValues: no partition can be specified when calculating with single fold")
      }
      postList=poolMcmcChains(hM$postList, start = start)
      pred = predict(hM, post=postList, Yc=Yc, mcmcStep=1, expected=expected)
      mpred = apply(abind(pred,along=3),c(1,2),mean)
      attr(mpred, "partition") = rep(1,hM$ny)
   } else{
      if(!is.null(nfolds) && !is.null(partition)){
         stop("HMSC.computePredictedValues: both nfolds and partition parameters cannot be specified similtaniously")
      }
      if(!is.null(nfolds)){
         if(hM$nr > 0){
            if(is.null(column)){
               indResidual = which(apply(hM$dfPi, 2, function(a) length(unique(a)) ) == self$ny)
               if(length(indResidual)==1){
                  column = indResidual
               } else{
                  if(length(indResidual)==0){
                     stop("HMSC.computePredictedValues: none of latent factors correspond to residual level, specify column parameter manually")
                  } else{
                     stop("HMSC.computePredictedValues: multiple latent factors correspond to residual level, specify column parameter manually")
                  }
               }
            }
            np = length(unique(hM$dfPi[,column]))
            if(np < nfolds){
               stop("HMSC.computePredictedValues: nfolds must be no bigger than numero of units in specified random level")
            }
            tmp1 = data.frame(col=unique(hM$dfPi[,column]), part=sample(rep(1:nfolds,ceiling(np/nfolds)),np))
            colnames(tmp1)[1] = colnames(hM$dfPi)[column]
            tmp2 = merge(hM$dfPi, tmp1, all.x=TRUE, sort=FALSE)
            part = tmp2[,ncol(tmp2)]
         }
         else{
            part=sample(rep(1:nfolds,ceiling(hM$ny/nfolds)),hM$ny)
         }
      }
      if(!is.null(partition)){
         if(length(partition) != hM$ny){
            stop("HMSC.computePredictedValues: partition parameter must be a vector of length ny")
         }
         nfolds = length(unique(partition))
         if(any(!(partition %in% 1:nfolds))){
            stop("HMSC.computePredictedValues: partition must be a vector with values in 1:N, each value represented")
         }
         part = partition
      }
      mpred = matrix(NA,hM$ny,hM$ns)
      for (k in 1:nfolds){
         print(sprintf("Cross-validation, fold %d out of %d", k, nfolds))
         train = (part!=k)
         val = (part==k)
         dfPi = as.data.frame(matrix(NA,sum(train),hM$nr))
         for(r in seq_len(hM$nr)){
            dfPi[,r] = factor(hM$dfPi[train,r])
         }
         hM1 = Hmsc(Y=as.matrix(hM$Y[train,]), X=as.matrix(hM$X[train,]), Xs=as.matrix(hM$Xs[train,]), Xv=as.matrix(hM$Xv[train,]), dist="probit", dfPi=dfPi, Tr=hM$Tr, C=hM$C, rL=hM$rL)
         hM1$distr=hM$distr
         # HOW TO BETTER SET THE DISTRIBUTION?
         # NEED TO INHERIT PRIORS, SCALINGS ETC. FROM m
         hM1 = sampleMcmc(hM1, samples=hM$samples, thin=hM$thin, transient=hM$transient,adaptNf=hM$adaptNf, nChains=1)
         postList = hM1$postList[[1]]
         postList = postList[start:length(postList)]
         dfPi = as.data.frame(matrix(NA,sum(val),hM$nr))
         for (r in seq_len(hM$nr)){
            dfPi[,r] = factor(hM$dfPi[val,r])
         }
         pred1 = predict(hM1, post=postList, X=as.matrix(hM$X[val,,drop=FALSE]), dfPiNew=dfPi, rL=hM$rL, Yc=Yc[val,,drop=FALSE], mcmcStep=mcmcStep, expected=expected)
         mpred1 = apply(abind(pred1,along=3),c(1,2),mean)
         mpred[val,] = mpred1
      }
      attr(mpred, "partition") = part
   }
   colnames(mpred) = hM$spNames
   return(mpred)
}
