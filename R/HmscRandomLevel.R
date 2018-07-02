#' @docType class
#' @title Single level of random factors in HMSC
#'
#' @description Specifies the structure for the single level of random factors, whether the level
#'  is assumed to be spatial or not, the spatial coordinates,
#'
#'
#' @param pi determines the unique IDs for the distinct units on this level of random factors
#' @param s matrix of coordinates in the
#' @param sDim number of spatial dimensions
#' @param N number of unique units on this level
#'
#' @examples
#' HmscRandomLevel$new(data.frame(s1=c(1:10),s2=c(10:1)))
#' HmscRandomLevel$new(pi=as.factor(1:10))
#'
#' @export

HmscRandomLevel <- R6::R6Class("HmscRandomLevel",
   public = list(
      pi = NULL,
      s = NULL,
      sDim = NULL,
      N = NULL,

      nfMax = NULL,
      nfMin = NULL,

      nu = NULL,
      a1 = NULL,
      b1 = NULL,
      a2 = NULL,
      b2 = NULL,
      alphapw = NULL,

      initialize = function(data=NULL, N=NULL, pi=NULL, priors=NULL){
         if(nargs()==0)
            stop("HmscRandomLevel: At least one argumnet should be specified")
         if(!is.null(data)){
            self$s = data
            self$N = nrow(data)
            self$pi = rownames(data)
            self$sDim = ncol(data)
         }
         if(!is.null(pi)){
            if(!is.null(self$pi) && self$pi != pi)
               stop("!!!Write some mistake output!!!. Specified names for units at latent factors' level must
                  conside with rows of data")
            self$pi = as.factor(pi)
            self$N = length(pi)
            self$sDim = 0
         }
         if(!is.null(N)){
            self$N = N
            self$pi = as.factor(1:N)
            self$sDim = 0
         }
         if(!is.null(priors)){
            self$setPriors(priors=priors, setDefault=TRUE)
         } else{
            self$setPriors(setDefault=TRUE)
         }
      },
      setPriors = function(priors=NULL, nu=NULL, a1=NULL, a2=NULL, b1=NULL, b2=NULL, alphapw=NULL, nfMax=NULL, nfMin=NULL, setDefault=FALSE){
         if(!is.null(nu)){
            self$nu = nu
         } else if(setDefault){
            self$nu = 3
         }
         if(!is.null(a1)){
            self$a1 = a1
         } else if(setDefault){
            self$a1 = 50
         }
         if(!is.null(b1)){
            self$b1 = b1
         } else if(setDefault){
            self$b1 = 1
         }
         if(!is.null(a2)){
            self$a2 = a2
         } else if(setDefault){
            self$a2 = 50
         }
         if(!is.null(b2)){
            self$b2 = b2
         } else if(setDefault){
            self$b2 = 1
         }
         if(!is.null(alphapw)){
            if(self$sDim == 0)
               stop("HmscRandomLevel.setPriors: prior for spatial scale was given, but not spatial coordinates were specified")
            if(ncol(alphapw)!=2)
               stop("HmscRandomLevel.setPriors: alphapw must be a matrix with two columns")
            self$alphapw = alphapw
         } else if(setDefault && self$sDim>0){
            alphaN = 100
            enclosingRectDiag = sqrt(sum(apply(self$s, 2, function(c) diff(range(c)))^2))
            self$alphapw = cbind(enclosingRectDiag*c(0:alphaN)/alphaN, c(0.5,rep(0.5/alphaN,alphaN)))
         }
         if(!is.null(nfMax)){
            self$nfMax = nfMax
         } else if(setDefault){
            self$nfMax = Inf
         }
         if(!is.null(nfMin)){
            if(nfMin > nfMax)
               stop("HmscRandomLevel.setPriors: nfMin must be not greater than nfMax")
            self$nfMin = nfMin
         } else if(setDefault){
            self$nfMin = 1
         }

      }
   )
)
