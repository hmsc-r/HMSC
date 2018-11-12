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
#'
#'
#' @seealso
#'
#' 
#' @examples
#'


getPostEstimate = function(parName, r=1, q=c(), chainIndex=1:length(self$postList)){
   m = self
   bind0 = function(...){
      abind(...,along=0)
   }
   postList = poolMcmcChains(m$postList, chainIndex=chainIndex)

   if(parName %in% c("Beta", "Gamma", "V", "sigma")){
      valList = lapply(postList, function(a) a[[parName]])
   }
   if(parName %in% c("Eta", "Lambda", "Psi", "Delta")){
      valList = lapply(postList, function(a) a[[parName]][[r]])
   }
   if(parName %in% c("Alpha")){
      valList = lapply(postList, function(a) m$rL[[r]]$alphapw[a[[parName]][[r]],1])
   }
   if(parName %in% c("Omega", "OmegaCor")){
      valList = lapply(postList, function(a) crossprod(a[["Lambda"]][[r]]))
      if(parName %in% c("OmegaCor")){
         valList = lapply(valList, function(Om) cov2cor(Om))
      }
   }
   val = do.call(bind0, valList)
   parDim = dim(valList[[1]])
   parDimLength = length(parDim)
   valMean = apply(val, 1+(1:parDimLength), mean)
   valSupport = apply(val>0, 1+(1:parDimLength), mean)
   res = list(mean=valMean, support=valSupport)
   if(length(q) > 0)
      res$q = apply(val, 1+(1:parDimLength), quantile, q)
   return(res)
}

Hmsc$set("public", "getPostEstimate", getPostEstimate, overwrite=TRUE)


