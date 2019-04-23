#' @title Hmsc$getPostEstimate
#'
#' @description Calculates mean, support and other posterior quantities for model parameter
#'
#' @param parName The name of the parameter to be summarized. Can take value of
#'  model's baseline parameters, "Omega" or "OmegaCor".
#' @param r Which random level calculate the parameter for. Has effect only for Eta, Lambda, Omega and OmegaCor.
#' @param q Which quantiles to calculate.
#' @param chainIndex Which posterir chains to use for summarization (defaults to all)
#'
#'
#' @return
#' A named list of posterior quantities. Element mean
#'
#' @examples
#'
#' @export


getPostEstimate = function(hM, parName, r=1, x=NULL, q=c(), chainIndex=1:length(hM$postList)){
   bind0 = function(...){
      abind(...,along=0)
   }
   postList = poolMcmcChains(hM$postList, chainIndex=chainIndex)

   if(parName %in% c("Beta", "Gamma", "V", "sigma")){
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
   res = list(mean=valMean, support=valSupport, supportNeg = valSupportNeg)
   if(length(q) > 0)
      res$q = apply(val, 1+(1:parDimLength), quantile, q)
   return(res)
}

