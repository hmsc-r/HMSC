#' @title setPriors.HmscKroneckerRandomLevel
#'
#' @description Sets or resets priors to the Hmsc object
#' @param rL a fitted \code{HmscSpatialRandomLevel} model object
#' @param nu,a1,b1,a2,b2 parameters of the multiplicative gamma process shrinking prior
#' @param nfMax maximum number of latent factors to be sampled
#' @param nfMin minimum number of latent factors to be sampled
#' @param alphaPrior #TODO fill-in
#' @param setDefault logical indicating whether default priors should be used
#' @param ... other arguments (ignored)
#'
#'
#' @return Modified HmscKroneckerRandomLevel object
#'
#' @importFrom matrixStats rowProds
#'
#' @export

setPriors.HmscKroneckerRandomLevel = function(rL, nu=NULL, a1=NULL, a2=NULL, b1=NULL, b2=NULL, nfMax=NULL, nfMin=NULL, alphaPrior=NULL, setDefault=FALSE)
{
   stopifnot(inherits(rL, "HmscKroneckerRandomLevel"))
   rL = setPriors.HmscRandomLevel(rL, nu=nu, a1=a1, a2=a2, b1=b1, b2=b2, nfMax=nfMax, nfMin=nfMin, setDefault)
   if(!is.null(alphaPrior)){
      rLSpatialVec = unlist(lapply(rL$rLList, function(rL) rL$sDim>0))
      if(length(alphaPrior)!=2 || length(alphaPrior)!=rL$kN)
         stop("HmscKroneckerRandomLevel.setPriors: alphaPrior must be either a list of length 2 or a list of length equal to number of kronecker elements")
      if(all(names(alphaPrior) == names(rL$rLList))){
         alphaPriorIsNullVec = unlist(lapply(alphaPrior, is.null))
         if(any(alphaPriorIsNullVec[rLSpatialVec]==TRUE))
            stop("HmscKroneckerRandomLevel.setPriors: alphaPrior elements for non-spatial levels must be alphapw-like 2-column matrices")
         if(any(alphaPriorIsNullVec[!rLSpatialVec]==FALSE))
            stop("HmscKroneckerRandomLevel.setPriors: alphaPrior elements for non-spatial levels must be NULL")
         alphaGridList = lapply(alphaPrior, function(a) a[,1])
         alphaProbList = lapply(alphaPrior, function(a) a[,2])
         alphaProbMatrix = as.matrix(expand.grid(alphaProbList[rLSpatialVec]))
         alphaProbVec = rowProds(alphaProbMatrix)
         dimVec = unlist(lapply(alphaPrior[rLSpatialVec], function(a) nrow(a)))
         dimNames = vector("list", sum(rLSpatialVec))
         names(dimNames) = names(alphaPrior)[rLSpatialVec]
         alphaProb = array(alphaProbVec, dimVec, dimNames)
      } else if(names(alphaPrior)==c("alphaGridList","alphaProb")){
         alphaGridList = alphaPrior$alphaGridList
         alphaProb = alphaPrior$alphaProb
      } else{
         stop("HmscKroneckerRandomLevel.setPriors: alphaPrior must be either a list of length 2 or a list of length equal to number of kronecker elements")
      }
      rL$alphaPrior = list(alphaGridList=alphaGridList, alphaProb=alphaProb)
   } else if(setDefault && rL$sDim>0){
      warning("HmscKroneckerRandomLevel.setPriors: default setting of alphaPrior value for HmscKroneckerRandomLevel is not implemented. No changes are made to existing prior.")
   }
   return(rL)
}
