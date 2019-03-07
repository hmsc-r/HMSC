#' @title Hmsc$evaluateModelFit
#'
#' @description Computes several measures of model fit
#' @param predY array of predictions, typically posterior sample
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

evaluateModelFit = function(hM, predY){

  computeRMSE = function(Y, predY){
    ny=dim(Y)[1]
    ns = dim(Y)[2]
    RMSES = rep(NA,ns)
    RMSEY = rep(NA,ny)
    RMSEA = NA
    for (i in 1:ns){
      RMSES[i] = sqrt(mean((Y[,i]-predY[,i])^2, na.rm=TRUE))
    }
    for (i in 1:ny){
      RMSEY[i] = sqrt(mean((Y[i,]-predY[i,])^2, na.rm=TRUE))
    }
    RMSEA=sqrt(mean((as.vector(Y)-as.vector(predY))^2, na.rm=TRUE))
    RMSE = list(S=RMSES, Y=RMSEY, A=RMSEA)
    return(RMSE)
  }

  computeR2 = function(Y, predY, method="pearson"){
    ny=dim(Y)[1]
    ns = dim(Y)[2]
    R2S = rep(NA,ns)
    R2Y = rep(NA,ny)
    R2A = NA
    for (i in 1:ns){
      R2S[i] = cor(Y[,i], predY[,i], method=method, use='pairwise')^2
    }
    for (i in 1:ny){
      R2Y[i] = cor(Y[i,], predY[i,], method=method, use='pairwise')^2
    }
    R2A=cor(as.vector(Y), as.vector(predY), method=method, use='pairwise')^2
    R2 = list(S=R2S, Y=R2Y, A=R2A)
    return(R2)
  }

  computeAUC = function(Y, predY){
    ny=dim(Y)[1]
    ns = dim(Y)[2]
    AUCS = rep(NA,ns)
    AUCY = rep(NA,ny)
    AUCA = NA
    for (i in 1:ns){
      if(length(unique(Y[,i]))==2) {AUCS[i] = auc(Y[,i],predY[,i])}
    }
    for (i in 1:ny){
      if(length(unique(Y[i,]))==2) {AUCY[i] = auc(Y[i,],predY[i,])}
    }
    resp = as.vector(Y)
    if (length(unique(resp))==2) {AUCA = auc(resp,as.vector(predY))}
    AUC = list(S=AUCS, Y=AUCY, A=AUCA)
    return(AUC)
  }

  computeTjurR2 = function(Y, predY){
    ny=dim(Y)[1]
    ns = dim(Y)[2]
    R2S = rep(NA,ns)
    R2Y = rep(NA,ny)
    R2A = NA
    for (i in 1:ns){
      R2S[i] = mean(predY[Y[,i]==1,i]) - mean(predY[Y[,i]==0,i])
    }
    for (i in 1:ny){
      R2Y[i] = mean(predY[i,Y[i,]==1]) - mean(predY[i,Y[i,]==0])
    }
    R2A = mean(as.vector(predY[as.vector(Y==1)])) - mean(as.vector(predY[as.vector(Y==0)]))
    R2 = list(S=R2S, Y=R2Y, A=R2A)
    return(R2)
  }

  median2 = function(x){return (median(x,na.rm=TRUE))}
  mean2 = function(x){return (mean(x,na.rm=TRUE))}

  MF = NA
  if(sum(hM$distr[,1]==3)==hM$ns){
    mPredY = apply(abind(predY,along=3),c(1,2),median2)
  } else {
    mPredY = apply(abind(predY,along=3),c(1,2),mean2)
  }

  if(sum(hM$distr[,1]==1)==hM$ns){
    R2 = computeR2(hM$Y, mPredY)
    RMSE = computeRMSE(hM$Y, mPredY)
    MF = list(R2=R2, RMSE=RMSE)
  }

  if(sum(hM$distr[,1]==2)==hM$ns){
    AUC = computeAUC(hM$Y, mPredY)
    TjurR2 = computeTjurR2(hM$Y, mPredY)
    RMSE = computeRMSE(hM$Y, mPredY)
    MF = list(AUC=AUC, TjurR2=TjurR2, RMSE=RMSE)
  }

  if(sum(hM$distr[,1]==3)==hM$ns){
    predO = 1*(predY>0)
    mPredO = apply(abind(predO,along=3),c(1,2),mean2)
    mPredCY = mPredY/mPredO
    SR2 = computeR2(hM$Y, mPredY, method="spearman")
    RMSE = computeRMSE(hM$Y, mPredY)
    O.AUC = computeAUC(1*(hM$Y>0), mPredO)
    O.TjurR2 = computeTjurR2(1*(hM$Y>0), mPredO)
    O.RMSE = computeRMSE(1*(hM$Y>0), mPredO)
    CY=hM$Y
    CY[CY==0]=NA
    C.SR2 = computeR2(CY, mPredCY, method="spearman")
    C.RMSE = computeRMSE(CY, mPredCY)
    MF = list(SR2=SR2, RMSE=RMSE, O.AUC=O.AUC, O.TjurR2=O.TjurR2, O.RMSE=O.RMSE, C.SR2=C.SR2, C.RMSE=C.RMSE)
  }

  return(MF)
}
