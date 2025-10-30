#' @title coralPredict
#'
#' @description Makes predictions from CORAL models
#'
#' @param XData dataframe of covariates for the sampling units to be predicted
#' @param XFormula a \code{\link{formula}}-class object for fixed effects
#' @param mu matrix of CORAL prior or posterior means
#' @param V matrix of CORAL prior or posterior flattened variances
#' @param prob whether to return predictions as probabilities (TRUE) or latent linear predictor (FALSE)
#'
#' @return
#' Matrix of CORAL predictions with predicted sampling units in the rows and species in the columns
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
