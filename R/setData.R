#' @title setData
#'
#' @description Sets or resets various data elements to an Hmsc object
#'
#' @param Y matrix of species occurences or abundances
#' @param X matrix of measured covariates
#' @param Pi matrix of correspondence between sampling units and units on different levels of latent factors
#' @param rL list of HmscRandomLevel objects, specifying the structure and data for random levels
#'
#' @examples
#'
#' @export

setData = function(Y=NULL, X=NULL, Pi=NULL, rL=NULL, Xs=NULL, Xv=NULL, dist="normal", spNames=NULL,
   trNames=NULL, covNames=NULL, levelNames=NULL){
   self$Y = as.matrix(Y)
   self$ny = nrow(Y)
   self$ns = ncol(Y)

   self$X = as.matrix(X)
   if(nrow(X) != ny){
      stop("Hmsc.setData: the number of rows in X should be equal to number of rows in Y")
   }
   self$nc = ncol(X)

   self$Pi = as.matrix(Pi)
   if(nrow(X) != ny){
      stop("Hmsc.setData: the number of rows in Pi should be equal to number of rows in Y")
   }
   self$nr = ncol(Pi)

   self$rL = rL

   self$Tr = matrix(1,ns,1)
   self$nt = 1

   if(is.null(spNames)){
      self$spNames = colnames(Y)
   } else{
      self$spNames = spNames
      colnames(self$Y) = spNames
   }

   if(is.null(covNames)){
      self$covNames = colnames(X)
   } else{
      if(length(covNames))
      self$covNames = covNames
      colnames(self$X) = covNames
   }

   if(is.null(levelNames)){
      self$levelNames = colnames(Pi)
   } else{
      self$levelNames = levelNames
      colnames(self$Pi) = levelNames
   }

   self$setPriors(setDefault=TRUE)
}

Hmsc$set("public", "setData", setData, overwrite=TRUE)

