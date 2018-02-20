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

Hmsc <- R6::R6Class("Hmsc",
   public = list(
      # data
      Y = NULL,
      X = NULL,
      rL = NULL,
      Xs = NULL,
      Xv = NULL,
      Tr = NULL,
      C = NULL,
      Pi = NULL,

      # dimensions
      ny = NULL,
      ns = NULL,
      nc = NULL,
      nr = NULL,
      nt = NULL,
      nf = NULL,
      ncs = NULL,
      ncv = NULL,

      # names
      spNames = NULL,
      covNames = NULL,
      trNames = NULL,
      levelNames = NULL,

      dist = NULL,

      # priors
      V0=NULL, f0=NULL,
      mGamma=NULL, UGamma=NULL,
      aSigma=NULL, bSigma=NULL,

      # sampling parameters
      samples = NULL, thin = NULL, adaptLf = NULL, saveToDisk=NULL,
      initPar = NULL, repN=NULL,
      randSeed = NULL,

      # posterior
      postList=NULL, repList=NULL,

      initialize = function(Y=NULL, X=NULL, Pi=NULL, rL=NULL, Xs=NULL, Xv=NULL, Tr=NULL, C=NULL, dist="normal", priors=NULL){
         # combine Hmsc and set data functions from Matlab
         self$setData(Y=Y, X=X, Pi=Pi, rL=rL)
      },
      setData = function(Y=NULL, X=NULL, Pi=NULL, rL=NULL, Xs=NULL, Xv=NULL, dist="normal", spNames=NULL,
         trNames=NULL, covNames=NULL, ...){},
      setPriors = function(priors=NULL){},
      setMcmcParameters = function(){},
      sampleMcmc = function(){},

      getPosterior = function(){
         # combines setPostThinning, saves a copy of that in the class and returns the posterior as mcmc object
      }
   ),
   private = list(
      computeInitialParameters = function(initPar){}
   )
)
