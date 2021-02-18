#' @title setPriors.Hmsc
#'
#' @description Sets or resets priors to the \code{Hmsc} object
#' @param hM a fitted \code{Hmsc} model object
#' @param V0 scale matrix in the Wishart prior distribution for the V matrix
#' @param f0 number of degrees of freedom in the Wishart prior distribution for the V matrix
#' @param mGamma mean for the prior multivariate Gaussian distribution for Gamma parameters
#' @param UGamma covariance matrix for the prior multivariate Gaussian distribution for Gamma parameters
#' @param aSigma shape parameter for the prior gamma distribution for the variance parameter, only for normal & lognormal Poisson models
#' @param bSigma rate parameter for the prior gamma distribution for the variance parameter, only for normal & lognormal Poisson models
#' @param nuRRR,a1RRR,b1RRR,a2RRR,b2RRR parameters of the multiplicative gamma process shrinking prior for reduced rank regression
#' @param rhopw discrete grid prior for phylogenetic signal, should be a matrix of 2 columns
#' @param setDefault logical indicating whether default priors should be used
#' @param \dots other parameters passed to the function.
#'
#' @return Modified \code{Hmsc} object
#'
#' @export

setPriors.Hmsc = function(hM, V0=NULL, f0=NULL, mGamma=NULL,
   UGamma=NULL, aSigma=NULL, bSigma=NULL, nuRRR=NULL, a1RRR=NULL,
   b1RRR=NULL, a2RRR=NULL, b2RRR=NULL, rhopw=NULL, setDefault=FALSE, ...){

   if(!is.null(V0)){
      if(!isSymmetric(V0) || nrow(V0) != hM$nc || ncol(V0) != hM$nc)
         stop("V0 must be a positive definite matrix of size equal to number of covariates nc")
      hM$V0 = V0
   } else if(setDefault){
      hM$V0 = diag(hM$nc)
   }

   if(!is.null(f0)){
      if(f0 < hM$nc)
         stop("f0 must be greater than number of covariates in the model nc")
      hM$f0 = f0
   } else if(setDefault){
      hM$f0 = hM$nc+1
   }

   if(!is.null(mGamma)){
      if(length(mGamma) != hM$nc*hM$nt)
         stop("mGamma must be a vector of length equal to number of covariates times traits: nc x nt")
      hM$mGamma = mGamma
   } else if(setDefault){
      hM$mGamma = rep(0, hM$nc*hM$nt)
   }

   if(!is.null(UGamma)){
      if(!isSymmetric(UGamma) || nrow(UGamma) != (hM$nc*hM$nt) || ncol(UGamma) != (hM$nc*hM$nt))
         stop("UGamma must be a positive definite matrix of size equal to nc x nt")
      hM$UGamma = UGamma
   } else if(setDefault){
      hM$UGamma = diag(hM$nc * hM$nt)
   }

   ## Default values for aSigma & bSigma: a vector of length 3 indexed
   ## by distr[,1] corresponding to "normal", never used ("probit"
   ## where sigma is not estimated), and "lognormal poisson"
   ## ("poisson" has the same distr[,1] index, but sigma is not
   ## estimated if distr[,2]==0)
   if(setDefault) {
      aSigmaDef <- c(1, 1, 0.5)
      bSigmaDef <- c(0.01, 0.01, 0.0001)
   }
   if(!is.null(aSigma)){
      hM$aSigma = aSigma
   } else if(setDefault){
      hM$aSigma = aSigmaDef[hM$distr[,1]]
   }

   if(!is.null(bSigma)){
      hM$bSigma = bSigma
   } else if(setDefault){
      hM$bSigma = bSigmaDef[hM$distr[,1]]
   }

   if(!is.null(rhopw)){
      if(is.null(hM$C))
         stop("prior for phylogeny given, but no phylogenic relationship matrix was specified")
      if(ncol(rhopw)!=2)
         stop("rhopw must be a matrix with two columns")
      hM$rhopw = rhopw
   } else if(setDefault){
      rhoN = 100
      hM$rhopw = cbind(c(0:rhoN)/rhoN, c(0.5,rep(0.5/rhoN,rhoN)))
   }
   if(!is.null(nuRRR)){
      hM$nuRRR = nuRRR
   } else if(setDefault){
      hM$nuRRR = 3
   }
   if(!is.null(a1RRR)){
      hM$a1RRR = a1RRR
   } else if(setDefault){
      hM$a1RRR = 1
   }
   if(!is.null(b1RRR)){
      hM$b1RRR = b1RRR
   } else if(setDefault){
      hM$b1RRR = 1
   }
   if(!is.null(a2RRR)){
      hM$a2RRR = a2RRR
   } else if(setDefault){
      hM$a2RRR = 50
   }
   if(!is.null(b2RRR)){
      hM$b2RRR = b2RRR
   } else if(setDefault){
      hM$b2RRR = 1
   }
   return(hM)
}

