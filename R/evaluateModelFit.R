#' @title evaluateModelFit
#'
#' @description Computes measures of model fit for a \code{Hmsc} model
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param predY array of predictions, typically posterior sample
#'
#' @return a list of measures of model fit
#'
#' @details All measures of model fit are based on comparing the posterior predictive distribution (\code{predY)})
#' to the observed values (\code{hM$Y}). The predicted distribution is first summarized to
#' a single matrix of predicted values by taking the posterior mean (for normal and probit models)
#' or posterior median (for Poisson models). All measures of model fit are given as vectors with
#' one value for each species.
#'
#' The kinds of measures of model fit depend on the type of response variable.
#' For all types of response variables, root-mean-square error (RMSE) between predicted
#' and observed values is computed.
#' For normal models, R2 is computed as squared pearson
#' correlation between observed and predicted values, times the sign of the correlation.
#' For probit models, Tjur R2 and AUC are
#' computed. For Poisson models, a pseudo-R2 is computed as
#' squared spearman correlation between observed and predicted values, times the sign of the correlation (SR2).
#' For Poisson models, the observed and predicted data are also truncated to occurrences (presence-absences),
#' for which the same measures are given as for the probit models (O.RMSE, O.AUC and O.TjurR2).
#' For Poisson models, the observed and predicted data are also subsetted to conditional on presence,
#' for which the root-mean-square error and pseudo-R2 based on squared spearman correlation
#' are computed (C.RMSE, C.SR2).
#'
#' The measures O.RMSE, O.AUC, O.TjurR2, C.RMSE and C.SR2 can be computed only if the option
#' \code{expected=FALSE} has been used when making the predictions
#'
#' If the model includes a mixture of response variable types, the resulting measures of model fit
#' contain NA's for those response variables for which they cannot be computed.
#'
#' @examples
#' # Evaluate model fit
#' preds = computePredictedValues(TD$m)
#' MF = evaluateModelFit(hM=TD$m, predY=preds)
#'
#' # Evaluate model performance based on cross validation: this will be slow
#' \dontrun{
#' partition = createPartition(TD$m, nfolds = 2)
#' predsCV1 = computePredictedValues(TD$m, partition=partition)
#' MF = evaluateModelFit(hM=TD$m, predY=predsCV1)
#' }
#' @importFrom stats cor median
#' @importFrom pROC auc
#' @importFrom abind abind
#'
#' @export

evaluateModelFit = function(hM, predY){

   computeRMSE = function(Y, predY){
      ns = dim(Y)[2]
      RMSE = rep(NA,ns)
      for (i in 1:ns){
         RMSE[i] = sqrt(mean((Y[,i]-predY[,i])^2, na.rm=TRUE))
      }
      return(RMSE)
   }

   computeR2 = function(Y, predY, method="pearson"){
      ns = dim(Y)[2]
      R2 = rep(NA,ns)
      for (i in 1:ns){
         co = cor(Y[,i], predY[,i], method=method, use='pairwise')
         R2[i] = sign(co)*co^2
      }
      return(R2)
   }

   computeAUC = function(Y, predY){
      ns = dim(Y)[2]
      AUC = rep(NA,ns)
      ## take care that Y has only levels {0,1} as specified in auc() below
      Y <- ifelse(Y > 0, 1, 0)
      for (i in 1:ns){
         sel = !is.na(Y[,i])
         if(length(unique(Y[sel,i]))==2)
             AUC[i] = auc(Y[sel,i],predY[sel,i], levels=c(0,1),direction="<")
      }
      return(AUC)
   }

   computeTjurR2 = function(Y, predY) {
      ns = dim(Y)[2]
      R2 = rep(NA, ns)
      for (i in 1:ns) {
         R2[i] = mean(predY[which(Y[, i] == 1), i]) - mean(predY[which(Y[,i] == 0), i])
      }
      return(R2)
   }

   median2 = function(x){return (median(x,na.rm=TRUE))}
   mean2 = function(x){return (mean(x,na.rm=TRUE))}

   mPredY = matrix(NA, nrow=hM$ny, ncol=hM$ns)
   sel = hM$distr[,1]==3
   if (sum(sel)>0){
       mPredY[,sel] = as.matrix(apply(abind(predY[,sel,,drop=FALSE], along=3),
                                      c(1,2), median2))
   }
   sel = !hM$distr[,1]==3
   if (sum(sel)>0){
       mPredY[,sel] = as.matrix(apply(abind(predY[,sel,,drop=FALSE], along=3),
                                      c(1,2), mean2))
   }

   RMSE = computeRMSE(hM$Y, mPredY)
   R2 = NULL
   AUC = NULL
   TjurR2 = NULL
   SR2 = NULL
   O.AUC = NULL
   O.TjurR2 = NULL
   O.RMSE = NULL
   C.SR2 = NULL
   C.RMSE = NULL

   sel = hM$distr[,1]==1
   if (sum(sel)>0){
      R2 = rep(NA,hM$ns)
      R2[sel] = computeR2(hM$Y[,sel,drop=FALSE],mPredY[,sel,drop=FALSE])
   }

   sel = hM$distr[,1]==2
   if (sum(sel)>0){
      AUC = rep(NA,hM$ns)
      TjurR2 = rep(NA,hM$ns)
      AUC[sel] = computeAUC(hM$Y[,sel,drop=FALSE], mPredY[,sel,drop=FALSE])
      TjurR2[sel] = computeTjurR2(hM$Y[,sel,drop=FALSE], mPredY[,sel,drop=FALSE])
   }

   sel = hM$distr[,1]==3
   if (sum(sel)>0){
      SR2 = rep(NA,hM$ns)
      O.AUC = rep(NA,hM$ns)
      O.TjurR2 = rep(NA,hM$ns)
      O.RMSE = rep(NA,hM$ns)
      C.SR2 = rep(NA,hM$ns)
      C.RMSE = rep(NA,hM$ns)
      SR2[sel] = computeR2(hM$Y[,sel,drop=FALSE],mPredY[,sel,drop=FALSE], method="spearman")
      predO = 1*(predY[,sel,,drop=FALSE]>0)
      mPredO = as.matrix(apply(abind(predO,along=3),c(1,2),mean2))
      O.AUC[sel] = computeAUC(1*(hM$Y[,sel,drop=FALSE]>0),mPredO)
      O.TjurR2[sel] = computeTjurR2(1*(hM$Y[,sel,drop=FALSE]>0), mPredO)
      O.RMSE[sel] = computeRMSE(1*(hM$Y[,sel,drop=FALSE]>0), mPredO)
      mPredCY = mPredY[,sel,drop=FALSE]/mPredO
      CY=hM$Y[,sel,drop=FALSE]
      CY[CY==0]=NA
      C.SR2[sel] = computeR2(CY, mPredCY, method="spearman")
      C.RMSE[sel] = computeRMSE(CY, mPredCY)
   }

   MF = list(RMSE=RMSE)
   if (!is.null(R2)){MF$R2 = R2}
   if (!is.null(AUC)){MF$AUC = AUC}
   if (!is.null(TjurR2)){MF$TjurR2 = TjurR2}
   if (!is.null(SR2)){MF$SR2 = SR2}
   if (!is.null(O.AUC)){MF$O.AUC = O.AUC}
   if (!is.null(O.TjurR2)){MF$O.TjurR2 = O.TjurR2}
   if (!is.null(O.RMSE)){MF$O.RMSE = O.RMSE}
   if (!is.null(C.SR2)){MF$C.SR2 = C.SR2}
   if (!is.null(C.RMSE)){MF$C.RMSE = C.RMSE}

   return(MF)
}
