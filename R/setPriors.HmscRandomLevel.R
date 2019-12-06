#' @title setPriors.HmscRandomLevel
#'
#' @description Sets or resets priors to the Hmsc object
#' @param rL a fitted \code{HmscRandomLevel} model object
#' @param nu,a1,b1,a2,b2 parameters of the multiplicative gamma process shrinking prior
#' @param alphapw discrete grid prior for spatial scale parameter
#' @param nfMax maximum number of latent factors to be sampled
#' @param nfMin minimum number of latent factors to be sampled
#' @param setDefault logical indicating whether default priors should be used
#' @param ... other arguments (ignored)
#'
#'
#' @return Modified HmscRandomLevel object
#'
#'
#' @export

setPriors.HmscRandomLevel = function(rL, nu=NULL, a1=NULL, a2=NULL, b1=NULL, b2=NULL, alphapw=NULL, nfMax=NULL, nfMin=NULL, setDefault=FALSE, ...)
{
   stopifnot(inherits(rL, "HmscRandomLevel"))
   xDim = max(rL$xDim, 1)
   if(!is.null(nu)){
      if(length(nu) == 1){
         rL$nu = rep(nu, xDim)
      } else{
         if(length(nu) == xDim){
            rL$nu = nu
         } else
            stop("HmscRandomLevel.setPriors: length of nu argument must be either 1 or rL$xDim")
      }
   } else if(setDefault){
      rL$nu = rep(3, xDim)
   }
   if(!is.null(a1)){
      if(length(a1) == 1){
         rL$a1 = rep(a1, xDim)
      } else{
         if(length(a1) == xDim){
            rL$a1 = a1
         } else
            stop("HmscRandomLevel.setPriors: length of a1 argument must be either 1 or rL$xDim")
      }
   } else if(setDefault){
      rL$a1 = rep(50, xDim)
   }
   if(!is.null(b1)){
      if(length(b1) == 1){
         rL$b1 = rep(b1, xDim)
      } else{
         if(length(b1) == xDim){
            rL$b1 = b1
         } else
            stop("HmscRandomLevel.setPriors: length of b1 argument must be either 1 or rL$xDim")
      }
   } else if(setDefault){
      rL$b1 = rep(1, xDim)
   }
   if(!is.null(a2)){
      if(length(a2) == 1){
         rL$a2 = rep(a2, xDim)
      } else{
         if(length(a2) == xDim){
            rL$a2 = a2
         } else
            stop("HmscRandomLevel.setPriors: length of a2 argument must be either 1 or rL$xDim")
      }
   } else if(setDefault){
      rL$a2 = rep(50, xDim)
   }
   if(!is.null(b2)){
      if(length(b2) == 1){
         rL$b2 = rep(b2, xDim)
      } else{
         if(length(b2) == xDim){
            rL$b2 = b2
         } else
            stop("HmscRandomLevel.setPriors: length of b2 argument must be either 1 or rL$xDim")
      }
   } else if(setDefault){
      rL$b2 = rep(1, xDim)
   }
   if(!is.null(alphapw)){
      if(rL$sDim == 0)
         stop("HmscRandomLevel.setPriors: prior for spatial scale was given, but not spatial coordinates were specified")
      if(ncol(alphapw)!=2)
         stop("HmscRandomLevel.setPriors: alphapw must be a matrix with two columns")
      rL$alphapw = alphapw
   } else if(setDefault && rL$sDim>0){
      alphaN = 100
      if(is.null(rL$distMat)){
         enclosingRectDiag = sqrt(sum(apply(rL$s, 2, function(c) diff(range(c)))^2))
      } else {
         enclosingRectDiag = max(rL$distMat)
      }
      rL$alphapw = cbind(enclosingRectDiag*c(0:alphaN)/alphaN, c(0.5,rep(0.5/alphaN,alphaN)))
   }
   if(!is.null(nfMax)){
      rL$nfMax = nfMax
   } else if(setDefault){
      rL$nfMax = Inf
   }
   if(!is.null(nfMin)){
      if(nfMin > nfMax)
         stop("HmscRandomLevel.setPriors: nfMin must be not greater than nfMax")
      rL$nfMin = nfMin
   } else if(setDefault){
      rL$nfMin = 2
   }
   return(rL)
}
