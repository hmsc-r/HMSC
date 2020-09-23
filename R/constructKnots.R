#' @title constructKnots
#'
#' @description Construct a Regular Grid of Knot Locations for Spatial GPP Model
#'
#' @details This is a helper function for spatial Hmsc models with the
#'     spatial method set to GPP where user must provide knot
#'     locations. Knot locations with a distance greater than
#'     \code{minKnotDist} to the nearest data point are dropped from
#'     the grid. If the input locations are
#'     \code{\link[sp]{SpatialPoints}} data, these are treated like
#'     Euclidean coordinates, and if the points are not projected, a
#'     warning is issued.
#'
#' @param sData a dataframe containing spatial or temporal coordinates of units of the random level
#' @param nKnots the number of knots wanted on the spatial dimension with the shortest range
#' @param knotDist the distance between the wanted knots
#' @param minKnotDist the minimum distance of a knot to the nearest data point
#'
#' @return a data frame with knot locations
#'
#' @details Only one of \code{nKnots} and \code{minKnotDist} arguments can be provided.
#'
#'
#' @examples
#' #Creating knots for some 2 dimensional spatial data
#' n = 100
#' xycoords = matrix(runif(2*n),ncol=2)
#' xyKnots = constructKnots(xycoords,knotDist = 0.2, minKnotDist = 0.5)
#'
#' @importFrom methods is
#' @importFrom FNN knnx.dist
#' @importFrom sp coordinates `coordinates<-` is.projected
#'     proj4string `proj4string<-`
#' @export

constructKnots =
    function(sData, nKnots = NULL, knotDist = NULL, minKnotDist = NULL)
{
   if(!is.null(nKnots) && !is.null(knotDist)){
      stop("nKnots and knotDist cannot both be specified")
   }
   ## get coordinates of spatial points, but warn if these are not projected
   if (is(sData, "Spatial")) {
      if (!is.projected(sData)) {
         warning("producing regular grid, but spatial points are not projected")
         proj4 <- proj4string(sData)
      }
      sData <- coordinates(sData)
   }
   mins = apply(sData,2,min)
   maxs = apply(sData,2,max)
   if(is.null(knotDist)){
      if(is.null(nKnots)){
         nKnots = 10
      }
      knotDist = min(maxs-mins)/nKnots
   }
   knotindices = list()
   for(d in 1:ncol(sData)){
      knotindices[[d]] = seq(mins[d],maxs[d],by=knotDist)
   }
   sKnot = expand.grid(knotindices)  #Full grid
   Dist = knnx.dist(sData,sKnot,k=1)
   if(is.null(minKnotDist)){
      minKnotDist = 2*knotDist
   }
   sKnot = sKnot[Dist < minKnotDist,]
   ## set original proj4string for non-projected Spatial data
   if (exists("proj4", inherits = FALSE)) {
      coordinates(sKnot) <- colnames(sKnot)
      proj4string(sKnot) <- proj4
   }
   sKnot
}
