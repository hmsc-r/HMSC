#' @title Create single level of random factors in HMSC
#'
#' @description Specifies the structure for the single level of random factors, whether the level
#'  is assumed to be spatial or not, the spatial coordinates and the covariate-dependent nonstationarity.
#'
#'
#' @param pi determines the unique IDs for the distinct units on this level of random factors
#' @param s matrix of coordinates in the
#' @param sDim number of spatial dimensions
#' @param N number of unique units on this level
#'
#' @return
#'
#' @examples
#' HmscRandomLevel$new(data.frame(s1=c(1:10),s2=c(10:1)))
#' HmscRandomLevel$new(pi=as.factor(1:10))
#'
#' @export

HmscRandomLevel = function(sData=NULL, distMat=NULL, XData=NULL, XFormula=NULL, xData=NULL, units=NULL, N=NULL, priors=NULL){
   rL = structure(list(pi=NULL, s=NULL, sDim=NULL, x=NULL, xDim=NULL, N=NULL, distMat=NULL, #
      nfMax=NULL, nfMin=NULL, nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL, alphapw=NULL), class="HmscRandomLevel")
   if(nargs()==0)
      stop("HmscRandomLevel: At least one argumnet should be specified")
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

   if(!is.null(priors)){
      rL = setPriors(rL, priors=priors, setDefault=TRUE)
   } else{
      rL = setPriors(rL, setDefault=TRUE)
   }
   return(rL)
}
