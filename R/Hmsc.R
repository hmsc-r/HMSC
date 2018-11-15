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
#' @return
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
      XData=NULL, XFormula=NULL, X=NULL, XScaled=NULL, XInterceptInd=NULL,
      rL = NULL,
      Xs = NULL,
      Xv = NULL,
      TrData=NULL, TrFormula=NULL, Tr=NULL, TrScaled=NULL, TrInterceptInd=NULL,
      C=NULL, phyloTree=NULL,
      dfPi = NULL,
      Pi = NULL,
      distr = NULL,

      # dimensions
      ny = NULL,
      ns = NULL,
      nc = NULL,
      nr = NULL,
      nt = NULL,
      nf = NULL,
      ncs = NULL,
      ncv = NULL,
      np = NULL,

      # names
      spNames = NULL,
      covNames = NULL,
      trNames = NULL,
      levelNames = NULL,

      # scaling
      XScalePar=NULL, TrScalePar=NULL,

      # priors
      V0=NULL, f0=NULL,
      mGamma=NULL, UGamma=NULL,
      aSigma=NULL, bSigma=NULL,
      nu=NULL, a1=NULL, b1=NULL, a2=NULL, b2=NULL,
      rhopw=NULL,

      # sampling parameters
      samples=NULL, transient=NULL, thin=NULL, adaptNf=NULL, saveToDisk=NULL,
      initPar=NULL, repN=NULL,
      randSeed=NULL,

      # posterior
      postList=NULL, repList=NULL,

      initialize = function(Y=NULL, XFormula=~., XData=NULL, X=NULL, XScale=TRUE,
         dfPi=NULL, rL=NULL, Xs=NULL, Xv=NULL,
         TrFormula=NULL, TrData=NULL, Tr=NULL, TrScale=TRUE,
         phyloTree=NULL, C=NULL,
         distr="normal", priors=NULL){
         self$setData(Y=Y, XFormula=XFormula, XData=XData, X=X, XScale=XScale,
            dfPi=dfPi, rL=rL,
            TrFormula=TrFormula, TrData=TrData, Tr=Tr,
            phyloTree=phyloTree, C=C,
            distr=distr)
      },
      setData = function(Y=NULL, XFormula=~., XData=NULL, X=NULL, XScale=TRUE,
         dfPi=NULL, rL=NULL, Xs=NULL, Xv=NULL,
         TrFormula=NULL, TrData=NULL, Tr=NULL, TrScale=TRUE,
         phyloTree=NULL, C=NULL,
         distr="normal", spNames=NULL,
         trNames=NULL, covNames=NULL, levelNames=NULL){},
      setPriors = function(priors=NULL){},
      sampleMcmc = function(){},

      computeAssociations = function(start=1){},
      computeMarginalEffects = function(ngrid=20, prob=c(0.025,0.5,0.975), clist=1:self$nc){},
      computePredictedValues = function(nfolds=1, start=1){},
      computeR2 = function(predY){},
      computeVariancePartitioning = function(group, groupnames, start=1){},
      convertToCodaObject = function(start=1, spNamesNumbers=c(TRUE,TRUE), covNamesNumbers=c(TRUE,TRUE), trNamesNumbers=c(TRUE,TRUE),
         Beta=TRUE, Gamma=TRUE, V=TRUE, Sigma=TRUE, Rho=TRUE, Eta=TRUE, Lambda=TRUE, Alpha=TRUE,
         Omega=TRUE, Psi=TRUE, Delta=TRUE){},
      plotMarginalEffects = function(pred, covariate, measure, index=1){},
      plotVariancePartitioning = function(VP){}
   ),
   private = list(
      computeInitialParameters = function(initPar){},
      computeDataParameters = function(){}
   )
)
