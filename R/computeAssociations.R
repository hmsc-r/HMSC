#' @title computeAssociations
#'
#' @description Computes the species association matrices associated with each random level
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param start index of first MCMC sample included
#' @param thin thinning interval of posterior distribution
#'
#' @return list of association matrices (\eqn{\omega}) corresponding to each random level in the model
#'
#' @examples
#' # Compute the associations (residual correlations) between species from a HMSC model
#' assoc = computeAssociations(TD$m)
#'
#' @importFrom abind abind
#'
#' @export

computeAssociations = function(hM, start=1, thin=1){
   OmegaCor = vector("list", hM$nr)
   postList=poolMcmcChains(hM$postList, start = start, thin = thin)
   getOmegaCor = function(a,r=r)
      return(cov2cor(crossprod(a$Lambda[[r]])))
   for (r in seq_len(hM$nr)){
      OmegaCor1 = lapply(postList, getOmegaCor, r=r)
      mOmegaCor1 = apply(abind(OmegaCor1,along=3),c(1,2),mean)
      OmegaCor2 = lapply(OmegaCor1, function(a) return(a>0))
      support1 = apply(abind(OmegaCor2,along=3),c(1,2),mean)
      colnames(mOmegaCor1) = hM$spNames
      rownames(mOmegaCor1) = hM$spNames
      colnames(support1) = hM$spNames
      rownames(support1) = hM$spNames
      tmp = list()
      tmp$mean = mOmegaCor1
      tmp$support = support1
      OmegaCor[[r]] = tmp
   }
   return(OmegaCor)
}
