#' @title Hmsc$computePredictedValues
#'
#' @description Computes ...
#' @param nfolds
#' @param start
#'
#' @examples
#'


computePredictedValues = function(nfolds=1, start=1){
   m = self
   if (nfolds==1){
      postList=poolMcmcChains(m$postList, start = start)
      pred = m$predict(post = postList, expected = TRUE)
      mpred = apply(abind(pred,along=3),c(1,2),mean)}
   else{
      mpred = matrix(NA,m$ny,m$ns)
      part = sample(rep(1:nfolds,ceiling(m$ny/nfolds)),m$ny,replace=FALSE)
      for (k in 1:nfolds){
         print(sprintf("Cross-validation, fold %d out of %d", k, nfolds))
         train = (part!=k)
         val = (part==k)
         dfPi = as.data.frame(matrix(NA,sum(train),m$nr))
         for (r in 1:m$nr){
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
         for (r in 1:m$nr){
            dfPi[,r] = factor(m$dfPi[val,r])
         }
         pred1 = m1$predict(post = postList, X=as.matrix(m$X[val,]), dfPiNew = dfPi, rL = m$rL, expected = TRUE)
         mpred1 = apply(abind(pred1,along=3),c(1,2),mean)
         mpred[val,] = mpred1
      }
   }
   colnames(mpred) = m$spNames
   return(mpred)
}

Hmsc$set("public", "computePredictedValues", computePredictedValues, overwrite=TRUE)
