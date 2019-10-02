#' @title Create an \code{Hmsc} random level
#'
#' @description Specifies the structure of a random factor, including whether the random factor is assumed to
#'   be spatially explicit or not, the spatial coordinates and the potential structure of covariate-dependent random factors.
#'
#' @param sData a dataframe containing spatial or temporal coordinates of units of the random level
#' @param sMethod a string specifying which spatial method to be used. Possible values are \code{Full}, \code{GPP} and \code{NNGP}
#' @param distMat a distance matrix containing the distances between units of the random level
#' @param xData a dataframe containing the covariates measured at the units of the random level for covariate-dependent
#'   associations
#' @param units a vector, specifying the names of the units of a non-structured level
#' @param N number of unique units on this level
#' @param nNeighbours a scalar specifying the number of neighbours to be used in case the spatial method is set to \code{NNGP}. Only positive values smaller than the total number of plots are allowed.
#' @param sKnot a dataframe containing the knot locations to be used for the gaussian predictive process if sMethod is set to \code{GPP}
#' @return a \code{HmscRandomLevel}-class object that can be used for \code{Hmsc}-class object construction
#'
#' @details Only one of \code{sData}, \code{distMat}, \code{xData}, \code{units} and \code{N} arguments can be
#'   provided.
#'
#'   As a good practice, we recommend to specify all available units for a random level, even if some of those are not
#'   used for training the model.
#'
#'
#' @seealso [setPriors.Hmsc()]
#'
#' @examples
#' # Setting a random level with 50 units
#' rL = HmscRandomLevel(units=TD$studyDesign$sample)
#'
#' # Setting a spatial random level
#' rL = HmscRandomLevel(sData=TD$xycoords)
#'
#' # Setting a covariate-dependent random level.
#' rL = HmscRandomLevel(xData=data.frame(x1=rep(1,length(TD$X$x1)),x2=TD$X$x2))
#'
#' @export

HmscRandomLevel = function(sData=NULL, sMethod = "Full", distMat=NULL, xData=NULL, units=NULL, N=NULL, nNeighbours=NULL, sKnot=NULL){
   rL = structure(list(pi=NULL, s=NULL, sDim=NULL, spatialMethod=NULL, x=NULL, xDim=NULL, N=NULL, distMat=NULL, #
      nfMax=NULL, nfMin=NULL, nNeighbours=NULL, nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL, alphapw=NULL), class="HmscRandomLevel")
   if(nargs()==0)
      stop("HmscRandomLevel: At least one argument must be specified")
   if(!is.null(distMat) && !is.null(sData)){
      stop("HmscRandomLevel: sData and distMat cannot both be specified")
   }
   if(!is.null(sData)){
      rL$s = sData
      rL$N = nrow(sData)
      rL$pi = sort(rownames(sData))
      rL$sDim = ncol(sData)
      rL$spatialMethod = sMethod
      rL$nNeighbours = nNeighbours
      rL$sKnot = sKnot
   } else
      rL$sDim = 0
   if(!is.null(distMat)){
      rL$distMat = distMat
      rL$N = nrow(distMat)
      rL$pi = sort(rownames(distMat))
      rL$spatialMethod = sMethod
      rL$sDim = Inf
   }
   if(!is.null(xData)){
      if(!is.null(rL$pi)){
         if(any(!(rownames(xData)%in%rL$pi)))
            stop("HmscRandomLevel: duplicated specification of unit names")
      } else{
         rL$pi = sort(rownames(xData))
         rL$N = nrow(xData)
      }
      rL$xDim = ncol(xData)
      rL$x = xData
   } else
      rL$xDim = 0

   if(!is.null(units)){
      if(!is.null(rL$pi))
         stop("HmscRandomLevel: duplicated specification of unit names")
      rL$pi = as.factor(units)
      rL$N = length(units)
      rL$sDim = 0
   }

   if(!is.null(N)){
      if(!is.null(rL$pi))
         stop("HmscRandomLevel: duplicated specification of the number of units")
      rL$N = N
      rL$pi = as.factor(1:N)
      rL$sDim = 0
   }

   rL = setPriors(rL, setDefault=TRUE)
   return(rL)
}
