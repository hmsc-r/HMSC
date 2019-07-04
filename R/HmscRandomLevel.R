#' @title Create single level of random factors in HMSC
#'
#' @description Specifies the structure for a single level of random factors, including whether the level is assumed to
#'   be spatial or not, the spatial coordinates and the potential structure of covariate-dependent nonstationarity.
#'
#'
#' @param sData a dataframe containing spatial or temporal coordinates of units of the random level
#' @param distMat a distance matrix containing the distances between units of the random level
#' @param xData a dataframe containing the covariates measured at the units of the random level for covariate-dependent
#'   associations
#' @param units a vector, specifying the names of the units of a non-structured level
#' @param N number of unique units on this level
#'
#' @return a \code{HmscRandomLevel}-class object that can be used for \code{Hmsc}-class object construction
#'
#' @details Only one of \code{sData}, \code{distMat}, \code{xData}, \code{units} and \code{N} arguments shall be
#'   provided (implmentation for \code{sData} and \code{xData} is coming later).
#'
#'   As a good practice we recommend to specify all available units for a random level, even if some of those are not
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
#' # Setting a covariate dependent random level.
#' rL = HmscRandomLevel(xData=data.frame(x1=rep(1,length(TD$X$x1)),x2=TD$X$x2))
#'
#' @export

HmscRandomLevel = function(sData=NULL, distMat=NULL, xData=NULL, units=NULL, N=NULL){
   rL = structure(list(pi=NULL, s=NULL, sDim=NULL, x=NULL, xDim=NULL, N=NULL, distMat=NULL, #
      nfMax=NULL, nfMin=NULL, nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL, alphapw=NULL), class="HmscRandomLevel")
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
