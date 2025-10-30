#' @title coralTrain
#'
#' @description Trains CORAL models
#'
#' @param Y covariates data.frame
#' @param XData a data frame of measured covariates for fixed effects, same as in the backbone \code{Hmsc(...)} constructor
#' @param XFormula a \code{\link{formula}}-class object for fixed effects
#' @param prior.mu matrix of CORAL prior means
#' @param prior.V matrix of CORAL prior variances
#' @param transient number of transient MCMC steps
#' @param samples number of posterior MCMC samples
#' @param thin thinning interval for MCMC samples
#'
#' @return
#' list with means and covariance matrices of CORAL posteriors
#'
#' @details
#' Arguments \code{transient}, \code{samples} and \code{thin} are passed to the \code{MCMCprobit(...)} call
#' sequentially made for for each species.
#'
#'
#' @importFrom stats as.formula
#' @importFrom MCMCpack MCMCprobit
#'
#' @export

coralTrain = function(Y, XData, XFormula, prior.mu, prior.V, transient, samples, thin){
   formula.coral = as.formula(paste("y", deparse1(as.formula(XFormula)), sep=""))
   beta.mean = matrix(NA, nrow(prior.mu), ncol(prior.mu))
   beta.var = matrix(NA, nrow(prior.V), ncol(prior.V))
   for(j in 1:nrow(prior.mu)){
      data.coral = cbind(data.frame(y=Y[,j]), XData)
      V = matrix(prior.V[j,], ncol(prior.mu), ncol(prior.mu))
      Q = chol2inv(chol(V))
      probit.coral = MCMCprobit(formula=formula.coral, data=data.coral,
                                burnin=transient, mcmc=samples, thin=thin, beta.start=0,
                                b0=prior.mu[j,], B0=Q)
      beta.mean[j,] = colMeans(probit.coral)
      beta.var[j,] = as.vector(var(probit.coral))
   }
   rownames(beta.var) = rownames(beta.mean) = colnames(Y)
   return(list(beta.mean=beta.mean, beta.var=beta.var))
}
