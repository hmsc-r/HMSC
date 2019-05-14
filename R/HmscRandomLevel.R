#' @title Create single level of random factors in HMSC
#'
#' @description Specifies the structure for a single level of random factors, including whether the level is assumed to
#'   be spatial or not, the spatial coordinates and the potential structure of covariate-dependent nonstationarity.
#'
#'
#' @param sData a dataframe containing spatial or temporal coordinates of units of the random level
#' @param SMethod a string specifying which spatial method to be used. Possible values are \code{Full}, \code{GPP} and \code{NNGP}
#' @param distMat a distance matrix containing the distances between units of the random level
#' @param xData a dataframe containing the covariates measured at the units of the random level for covariate-dependent
#'   associations
#' @param units a vector, specifying the names of the units of a non-structured level
#' @param N number of unique units on this level
#' @param nNeighbours a scalar specifying the number of neighbours to be used in case the spatial method is set to \code{NNGP}. Only positive values smaller than the total number of plots are allowed.
#'
#' @return a \code{HmscRandomLevel}-class object that can be used for \code{Hmsc}-class object construction
#'
#' @details Only one of \code{sData}, \code{distMat}, \code{xData}, \code{units} and \code{N} arguments shall be
#'   provided (implmentation for \code{sData} and \code{xData} is coming later). Implementation for \code{GPP} is coming later
#'
#'   As a good practice we recommend to specify all available units for a random level, even if some of those are not
#'   used for training the model.
#'
#'
#' @seealso [setPriors.Hmsc()]
#'
#' @examples
#' rL = HmscRandomLevel(sData=data.frame(s1=c(1:10),s2=c(10:1)))
#' rL = HmscRandomLevel(units=as.factor(1:10))
#'
#' @export

HmscRandomLevel = function(sData=NULL, sMethod = "Full", distMat=NULL, xData=NULL, units=NULL, N=NULL, nNeighbours=NULL){
   rL = structure(list(pi=NULL, s=NULL, sDim=NULL, spatialMethod=NULL, x=NULL, xDim=NULL, N=NULL, distMat=NULL, #
      nfMax=NULL, nfMin=NULL, nNeighbours=NULL, nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL, alphapw=NULL), class="HmscRandomLevel")
   if(nargs()==0)
      stop("HmscRandomLevel: At least one argument should be specified")
   if(!is.null(distMat) && !is.null(sData)){
      stop("HmscRandomLevel: both sData and distMat arguments cannot be specified")
   }
   if(!is.null(sData)){
      rL$s = sData
      rL$N = nrow(sData)
      rL$pi = sort(rownames(sData))
      rL$sDim = ncol(sData)
      rL$spatialMethod = sMethod
      rL$nNeighbours = nNeighbours
   } else
      rL$sDim = 0
   if(!is.null(distMat)){
      rL$distMat = distMat
      rL$N = nrow(distMat)
      rL$pi = sort(rownames(distMat))
      rL$sDim = Inf
   }
   if(!is.null(xData)){
      if(!is.null(rL$pi)){
         if(any(!(rownames(xData)%in%rL$pi)))
            stop("HmscRandomLevel: duplicated specification of units names")
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
         stop("HmscRandomLevel: duplicated specification of units names")
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
