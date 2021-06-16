#' @title setPriors.HmscSpatialRandomLevel
#'
#' @description Sets or resets priors to the Hmsc object
#' @param rL a fitted \code{HmscSpatialRandomLevel} model object
#' @param nu,a1,b1,a2,b2 parameters of the multiplicative gamma process shrinking prior
#' @param alphapw discrete grid prior for spatial scale parameter
#' @param nfMax maximum number of latent factors to be sampled
#' @param nfMin minimum number of latent factors to be sampled
#' @param setDefault logical indicating whether default priors should be used
#' @param ... other arguments (ignored)
#'
#'
#' @return Modified HmscSpatialRandomLevel object
#'
#' @importFrom methods is
#' @importFrom sp bbox coordinates `coordinates<-` proj4string
#'     `proj4string<-` spDists
#'
#' @export

setPriors.HmscSpatialRandomLevel = function(rL, nu=NULL, a1=NULL, a2=NULL, b1=NULL, b2=NULL, nfMax=NULL, nfMin=NULL, alphapw=NULL, setDefault=FALSE)
{
   stopifnot(inherits(rL, "HmscSpatialRandomLevel"))
   rL = setPriors.HmscRandomLevel(rL, nu=nu, a1=a1, a2=a2, b1=b1, b2=b2, nfMax=nfMax, nfMin=nfMin, setDefault)
   if(!is.null(alphapw)){
      if(rL$sDim == 0)
         stop("HmscSpatialRandomLevel.setPriors: prior for spatial scale was given, but not spatial coordinates were specified")
      if(ncol(alphapw)!=2)
         stop("HmscSpatialRandomLevel.setPriors: alphapw must be a matrix with two columns")
      rL$alphapw = alphapw
   } else if(setDefault && rL$sDim>0){
      alphaN = 100
      if(is.null(rL$distMat)){
         if (is(rL$s, "Spatial")) {
            ## find diagonal from the bounding box instead of
            ## evaluating all spatial distances (that can be a huge
            ## task) similarly as with non-spatial points
            enclosingRect <- as.data.frame(t(bbox(rL$s)))
            coordinates(enclosingRect) <- colnames(enclosingRect)
            proj4string(enclosingRect) <- proj4string(rL$s)
            enclosingRectDiag <- max(spDists(enclosingRect))
         } else {
            enclosingRectDiag = sqrt(sum(apply(rL$s, 2, function(c) diff(range(c)))^2))
         }
      } else {
         enclosingRectDiag = max(rL$distMat)
      }
      rL$alphapw = cbind(enclosingRectDiag*c(0:alphaN)/alphaN, c(0.5,rep(0.5/alphaN,alphaN)))
   }
   return(rL)
}
