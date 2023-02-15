#' @title getPostEstimate
#'
#' @description Calculates mean, support and other posterior quantities for a specified model parameter
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param parName name of the parameter to be summarized. Can take value of
#'  model's baseline parameters, "Omega" or "OmegaCor".
#' @param r the random level for which to calculate the parameter. Has effect only for Eta, Lambda, Omega and OmegaCor.
#' @param x values of covariates for covariate dependent omega
#' @param q vector of quantiles to calculate.
#' @param chainIndex which posterior chains to use for summarization (defaults to all)
#' @param start index of first MCMC sample included
#' @param thin thinning interval of posterior distribution
#'
#'
#' @return
#' A named list of posterior quantities.
#'
#' @examples
#' # Get posterior mean and support for species' responses to environmental covariates
#' postBeta = getPostEstimate(TD$m, parName='Beta')
#'
#' # Get posterior mean and support for species' responses to latent factors for the first random level
#' postLambda = getPostEstimate(TD$m, parName='Lambda', r=1)
#'
#' @importFrom stats cov2cor
#' @importFrom abind abind
#'
#' @export


getPostEstimate = function(hM, parName, r=1, x=NULL, q=c(), chainIndex=1:length(hM$postList), start=1, thin=1){
   bind0 = function(...){
      abind(...,along=0)
   }
   postList = poolMcmcChains(hM$postList, chainIndex=chainIndex, start=start, thin=thin)

   if(parName %in% c("Beta", "Gamma", "V", "rho", "sigma","wRRR")){
      valList = lapply(postList, function(a) a[[parName]])
   }
   if(parName %in% c("Eta", "Lambda", "Psi", "Delta")){
      valList = lapply(postList, function(a) a[[parName]][[r]])
   }
   if(parName %in% c("Alpha")){
      valList = lapply(postList, function(a) hM$rL[[r]]$alphapw[a[[parName]][[r]],1])
   }
   if(parName %in% c("Omega", "OmegaCor")){
      if(hM$ranLevels[[r]]$xDim == 0){
         valList = lapply(postList, function(a) crossprod(a[["Lambda"]][[r]]))
      } else{
         if(is.null(x))
            x = c(1, rep(0, hM$ranLevels[[r]]$xDim-1))
         getOmega = function(a){
            lambda = a[["Lambda"]][[r]]
            dim = dim(lambda)
            return( crossprod(rowSums(lambda * array(rep(x, each=prod(dim[1:2])),dim), dims=2)) )
         }
         valList = lapply(postList, getOmega)
      }
      if(parName %in% c("OmegaCor")){
         valList = lapply(valList, function(Om) cov2cor(Om))
      }
   }
   val = do.call(bind0, valList)
   parDim = dim(valList[[1]])
   parDimLength = length(parDim)
   valMean = colMeans(val)
   valSupport = colMeans(val>0)
   valSupportNeg = colMeans(val<0)
   if(parName=="Beta"){
      colnames(valMean) = hM$spNames
      colnames(valSupport) = hM$spNames
      colnames(valSupportNeg) = hM$spNames
   }
   res = list(mean=valMean, support=valSupport, supportNeg = valSupportNeg)
   if(length(q) > 0)
      res$q = apply(val, 1+(1:parDimLength), quantile, q)
   return(res)
}

