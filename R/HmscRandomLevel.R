#' @title Create an \code{Hmsc} random level
#'
#' @description Specifies the structure of a random factor, including
#'     whether the random factor is assumed to be spatially explicit
#'     or not, the spatial coordinates and the potential structure of
#'     covariate-dependent random factors.
#' @param sData a matrix or a dataframe containing spatial or temporal
#'     coordinates of units of the random level, or a similar
#'     \code{SpatialPoints} structure of the \CRANpkg{sp} package. If
#'     spatial coordinates are unprojected longitude and latitude,
#'     great circle distances will be calculated internally. All
#'     spatial locations should be unique. If you have several
#'     observations in the same point, they should be identified by
#'     the random levels.
#' @param sMethod a string specifying which spatial method to be
#'     used. Possible values are \code{"Full"}, \code{"GPP"} and
#'     \code{"NNGP"}
#' @param distMat a distance matrix containing the distances between
#'     units of the random level, with unit names as rownames, or a
#'     \code{\link{dist}} structure with location
#'     Labels. \code{distMat} cannot be used with \code{"GPP"} spatial
#'     model.
#' @param xData a dataframe containing the covariates measured at the
#'     units of the random level for covariate-dependent associations
#' @param units a vector, specifying the names of the units of a
#'     non-structured level
#' @param N number of unique units on this level
#' @param nNeighbours a scalar specifying the number of neighbours to
#'     be used in case the spatial method is set to \code{NNGP}. Only
#'     positive values smaller than the total number of plots are
#'     allowed.
#' @param sKnot a dataframe containing the knot locations to be used
#'     for the Gaussian predictive process if \code{sMethod} is set to
#'     \code{"GPP"}. Suitable data can be produced with
#'     \code{\link{constructKnots}}. The knot locations shall not
#'     duplicate \code{sData}.
#'
#' @param longlat Interpret coordinate data \code{sData} as longitude
#'     and latitude in decimal degrees. If this is \code{TRUE}, great
#'     circle distances will be used instead of Euclidean distances.
#'
#' @return a \code{HmscRandomLevel}-class object that can be used for \code{Hmsc}-class object construction
#'
#' @details Only one of \code{sData}, \code{distMat}, \code{xData}, \code{units} and \code{N} arguments can be
#'   provided.
#'
#'   As a good practice, we recommend to specify all available units for a random level, even if some of those are not
#'   used for training the model.
#'
#'
#' @seealso \code{\link{setPriors.Hmsc}} to change the default priors
#'     of an existing \code{HmscRandomLevel} object.
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
#' @importFrom methods is
#' @importFrom sp coordinates `coordinates<-` CRS is.projected
#'
#' @export

HmscRandomLevel =
    function(sData=NULL, sMethod = "Full", distMat=NULL, xData=NULL, units=NULL,
             N=NULL, nNeighbours=10, sKnot=NULL, longlat = FALSE)
{
   rL = structure(list(pi=NULL, s=NULL, sDim=NULL, spatialMethod=NULL, x=NULL, xDim=NULL, N=NULL, distMat=NULL, #
      nfMax=NULL, nfMin=NULL, nNeighbours=NULL, nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL, alphapw=NULL), class="HmscRandomLevel")
   if(nargs()==0)
      stop("at least one argument must be specified")
   if(!is.null(distMat) && !is.null(sData)){
      stop("sData and distMat cannot both be specified")
   }
   if(!is.null(sData)){
      ## longitude & latitude data?
      if (longlat) {
         sData <- as.data.frame(sData)
         coordinates(sData) <- colnames(sData)
         proj4string(sData) <- CRS("+proj=longlat")
      }
      ## Retain Spatial data if they are non-projected (longlat),
      ## otherwise extract only coordinates
      if (is(sData, "Spatial") && is.projected(sData)) {
         sData <- coordinates(sData) # no longer "Spatial"
      }
      rL$s = sData
      ## several standard functions do not work with Spatial (sp) data
      if (is(sData, "Spatial")) {
         if (any(duplicated(coordinates(sData))))
            stop("sData locations must be unique")
         rL$N <- nrow(coordinates(sData))
         rL$pi <- as.factor(sort(row.names(sData)))
         rL$sDim <- ncol(coordinates(sData))
      } else {
         if (any(duplicated(sData)))
            stop("sData locations must be unique")
         rL$N = nrow(sData)
         rL$pi = as.factor(sort(rownames(sData)))
         rL$sDim = ncol(sData)
      }
      rL$spatialMethod = sMethod
      rL$nNeighbours = nNeighbours
      rL$sKnot = sKnot
   } else
       rL$sDim = 0
   ## we test against duplicated location in sData, but not here: zero
   ## distances between locations give a rank-deficit distance
   ## matrix with error in computeDataParameters
   if(!is.null(distMat)) {
      if (inherits(distMat, "dist"))
         distMat <- as.matrix(distMat)
      rL$distMat = distMat
      rL$N = nrow(distMat)
      rL$pi = as.factor(sort(rownames(distMat)))
      if (length(rL$pi) == 0) # rownames missing
         stop("'distMat' should have rownames for random levels")
      rL$spatialMethod = sMethod
      if (sMethod == "NNGP")
          rL$nNeighbours = nNeighbours
      rL$sDim = Inf
   }
   ## check that data are adequate for sMethod
   if (sMethod == "GPP" && is.null(rL$s))
       stop("sMethod GPP needs sData of coordinates")
   if (sMethod == "GPP" && is.null(rL$sKnot))
       stop("sMethod GPP needs sKnot of coordinates")

   if(!is.null(xData)){
      if(!is.null(rL$pi)){
         if(any(!(rownames(xData)%in%rL$pi)))
            stop("duplicated specification of unit names")
      } else{
         rL$pi = sort(as.factor(rownames(xData)))
         rL$N = nrow(xData)
      }
      rL$xDim = ncol(xData)
      rL$x = xData
   } else
      rL$xDim = 0

   if(!is.null(units)){
      if(!is.null(rL$pi))
         stop("duplicated specification of unit names")
      rL$pi = as.factor(units)
      rL$N = length(units)
      rL$sDim = 0
   }

   if(!is.null(N)){
      if(!is.null(rL$pi))
         stop("duplicated specification of the number of units")
      rL$N = N
      rL$pi = as.factor(1:N)
      rL$sDim = 0
   }

   rL = setPriors(rL, setDefault=TRUE)
   rL$call <- match.call()
   rL
}
