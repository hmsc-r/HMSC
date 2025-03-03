#' @title coralPredict
#'
#' @description Makes predictions from CORAL models
#'
#' @param XData covariates data.frame
#' @param XFormula arg2
#' @param mu arg3
#' @param V arg4
#' @param prob arg5
#'
#' @return
#' Matrix of size ny x ns with CORAL predictions
#'
#' @export

coralPredict = function(XData, XFormula, mu, V, prob=TRUE){
   # predicts from CORAL model using means and variances of species-specific estimated posterior MVN approximation
   ns = nrow(mu)
   X = model.matrix(XFormula, XData)
   L = X %*% t(mu)
   V = array(V, c(nrow(mu),ncol(mu),ncol(mu)))
   if(prob == TRUE){
      LHat = matrix(NA, nrow(XData), ns)
      for(j in 1:ns){
         LHat[,j] = L[,j] / sqrt(1 + rowSums((X%*%V[j,,])*X))
      }
      P = pnorm(LHat)
   } else{
      P = L
   }
   rownames(P) = rownames(XData)
   colnames(P) = rownames(mu)
   return(P)
}
