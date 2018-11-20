#' @title setPriors.HmscRandomLevel
#'
#' @description Sets or resets priors to the Hmsc object
#'
#' @param priors
#'
#'
#' @return Modified HmscRandomLevel object
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export

setPriors.HmscRandomLevel = function(rL, priors=NULL, nu=NULL, a1=NULL, a2=NULL, b1=NULL, b2=NULL, alphapw=NULL, nfMax=NULL, nfMin=NULL, setDefault=FALSE){
   stopifnot(class(rL) == "HmscRandomLevel")

   if(!is.null(nu)){
      rL$nu = nu
   } else if(setDefault){
      rL$nu = 3
   }
   if(!is.null(a1)){
      rL$a1 = a1
   } else if(setDefault){
      rL$a1 = 50
   }
   if(!is.null(b1)){
      rL$b1 = b1
   } else if(setDefault){
      rL$b1 = 1
   }
   if(!is.null(a2)){
      rL$a2 = a2
   } else if(setDefault){
      rL$a2 = 50
   }
   if(!is.null(b2)){
      rL$b2 = b2
   } else if(setDefault){
      rL$b2 = 1
   }
   if(!is.null(alphapw)){
      if(rL$sDim == 0)
         stop("HmscRandomLevel.setPriors: prior for spatial scale was given, but not spatial coordinates were specified")
      if(ncol(alphapw)!=2)
         stop("HmscRandomLevel.setPriors: alphapw must be a matrix with two columns")
      rL$alphapw = alphapw
   } else if(setDefault && rL$sDim>0){
      alphaN = 100
      enclosingRectDiag = sqrt(sum(apply(rL$s, 2, function(c) diff(range(c)))^2))
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
      rL$nfMin = 1
   }
   return(rL)
}
