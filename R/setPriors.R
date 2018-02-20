#' @title setPriors
#'
#' @description Sets or resets priors to the Hmsc object
#'
#' @param V0 scale matrix in the Wishart prior distribution for the Beta parameters
#' @param f0 number of degreees of freedom in the Wishart prior distribution for the Beta parameters
#' @param mGamma mean for the prior multivariate Gaussian distribution for Gamma parameters
#' @param UGamma covariance matrix for the prior multivariate Gaussian distribution for Gamma parameters
#'
#' @examples
#'
#' @export

setPriors = function(V0=NULL, f0=NULL, mGamma=NULL, UGamma=NULL, aSigma=NULL, bSigma=NULL, setDefault=FALSE){
   if(!is.null(V0)){
      if(!isSymmetric(V0) || nrow(V0) != self$nc || ncol(V0) != self$nc)
         stop("HMSC.setPriors: V0 must be a positive definite matrix of size equal to number of covariates nc")
      self$V0 = V0
   } else if(setDefault){
      self$V0 = diag(self$nc)
   }

   if(!is.null(f0)){
      if(f0 < self$nc)
         stop("HMSC.setPriors: f0 must be greater than number of covariates in the model nc")
      self$f0 = f0
   } else if(setDefault){
      self$f0 = self$nc+1
   }

   if(!is.null(mGamma)){
      if(length(mGamma) != self$nc)
         stop("HMSC.setPriors: mGamma must be a vector of length equal to number of covariates nc")
      self$mGamma = mGamma
   } else if(setDefault){
      self$mGamma = rep(0, self$nc*self$nt)
   }

   if(!is.null(mGamma)){
      if(!isSymmetric(UGamma) || nrow(UGamma) != self$nc || ncol(UGamma) != self$nc)
         stop("HMSC.setPriors: UGamma must be a positive definite matrix of size equal to number of covariates nc")
      self$UGamma = UGamma
   } else if(setDefault){
      self$UGamma = diag(self$nc)
   }

   if(!is.null(aSigma)){
      self$aSigma = aSigma
   } else if(setDefault){
      self$aSigma = rep(1, self$ns)
   }

   if(!is.null(bSigma)){
      self$bSigma = bSigma
   } else if(setDefault){
      self$bSigma = rep(0.3, self$ns)
   }
   # print("Hi, it is setPriors!")
}

Hmsc$set("public", "setPriors", setPriors, overwrite=TRUE)

