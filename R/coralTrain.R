#' @title coralTrain
#'
#' @description Trains CORAL models
#'
#' @param Y covariates data.frame
#' @param XData arg2
#' @param XFormula arg3
#' @param prior.mu arg4
#' @param prior.V arg5
#' @param transient arg6
#' @param samples arg7
#' @param thin arg8
#'
#' @return
#' list with means and covariance matrices of CORAL posteriors
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
