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
      priors = NULL,

      initialize = function(data=NULL, N=NULL, pi=NULL, priors=NULL){
         if(nargs()==0)
            stop("At least one argumnet should be specified")
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
      },
      setPriors = function(prior=NULL, mu=NULL, a1=NULL, a2=NULL, b1=NULL, b2=NULL, alphapw=NULL, nfFix=NULL, nfMax=NULL, nfMin=NULL){

      }
   )
)
