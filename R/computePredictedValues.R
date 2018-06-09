#' @title Hmsc$computePredictedValues
#'
#' @description Computes ...
#' @param nfolds
#' @param start
#'
#' @examples
#'


computePredictedValues = function(nfolds=NULL, column=NULL, partition=NULL, start=1){
   m = self
   if(identical(nfolds,1) || (is.null(nfolds) && is.null(partition))){
      if(!is.null(column)){
         stop("HMSC.computePredictedValues: no column can be specified when calculating with single fold")
      }
      if(!is.null(partition)){
         stop("HMSC.computePredictedValues: no partition can be specified when calculating with single fold")
      }
      postList=poolMcmcChains(m$postList, start = start)
      pred = m$predict(post = postList, expected = TRUE)
      mpred = apply(abind(pred,along=3),c(1,2),mean)
      attr(mpred, "partition") = rep(1,m$ny)
   } else{
      if(!is.null(nfolds) && !is.null(partition)){
         stop("HMSC.computePredictedValues: both nfolds and partition parameters cannot be specified similtaniously")
      }
      if(!is.null(nfolds)){
         if(m$nr > 0){
            if(is.null(column)){
               indResidual = which(apply(m$dfPi, 2, function(a) length(unique(a)) ) == self$ny)
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
            np = length(unique(m$dfPi[,column]))
            if(np < nfolds){
               stop("HMSC.computePredictedValues: nfolds must be no bigger than numero of units in specified random level")
            }
            tmp1 = data.frame(col=unique(m$dfPi[,column]), part=sample(rep(1:nfolds,ceiling(np/nfolds)),np))
            colnames(tmp1)[1] = colnames(m$dfPi)[column]
            tmp2 = merge(m$dfPi, tmp1, all.x=TRUE, sort=FALSE)
            part = tmp2[,ncol(tmp2)]
         }
         else{
            part=sample(rep(1:nfolds,ceiling(m$ny/nfolds)),m$ny)
         }
      }
      if(!is.null(partition)){
         if(length(partition) != m$ny){
            stop("HMSC.computePredictedValues: partition parameter must be a vector of length ny")
         }
         nfolds = length(unique(partition))
         if(any(!(partition %in% 1:nfolds))){
            stop("HMSC.computePredictedValues: partition must be a vector with values in 1:N, each value represented")
         }
         part = partition
      }
      mpred = matrix(NA,m$ny,m$ns)
      for (k in 1:nfolds){
         print(sprintf("Cross-validation, fold %d out of %d", k, nfolds))
         train = (part!=k)
         val = (part==k)
         dfPi = as.data.frame(matrix(NA,sum(train),m$nr))
         for(r in seq_len(m$nr)){
            dfPi[,r] = factor(m$dfPi[train,r])
         }
         m1 = Hmsc$new(Y=as.matrix(m$Y[train,]), X=as.matrix(m$X[train,]), Xs=as.matrix(m$Xs[train,]), Xv=as.matrix(m$Xv[train,]), dist="probit", dfPi=dfPi, Tr=m$Tr, C=m$C, rL=m$rL)
         m1$distr=m$distr
         # HOW TO BETTER SET THE DISTRIBUTION?
         # NEED TO INHERIT PRIORS, SCALINGS ETC. FROM m
         m1$sampleMcmc(samples=m$samples, thin=m$thin, transient=m$transient,adaptNf=m$adaptNf, nChains = 1)
         postList = m1$postList[[1]]
         postList = postList[start:length(postList)]
         dfPi = as.data.frame(matrix(NA,sum(val),m$nr))
         for (r in seq_len(m$nr)){
            dfPi[,r] = factor(m$dfPi[val,r])
         }
         pred1 = m1$predict(post = postList, X=as.matrix(m$X[val,]), dfPiNew = dfPi, rL = m$rL, expected = TRUE)
         mpred1 = apply(abind(pred1,along=3),c(1,2),mean)
         mpred[val,] = mpred1
      }
      attr(mpred, "partition") = part
   }
   colnames(mpred) = m$spNames
   return(mpred)
}

Hmsc$set("public", "computePredictedValues", computePredictedValues, overwrite=TRUE)
