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

setData = function(Y=NULL, X=NULL, Pi=NULL, rL=NULL, Xs=NULL, Xv=NULL, distr="normal", spNames=NULL,
   trNames=NULL, covNames=NULL, levelNames=NULL){
   self$Y = as.matrix(Y)
   self$ny = nrow(Y)
   self$ns = ncol(Y)

   if(nrow(X) != ny){
      stop("Hmsc.setData: the number of rows in X should be equal to number of rows in Y")
   }
   self$X = as.matrix(X)
   self$nc = ncol(X)

   if(nrow(Pi) != ny){
      stop("Hmsc.setData: the number of rows in Pi should be equal to number of rows in Y")
   }
   self$Pi = matrix(as.numeric(as.matrix(Pi)),nrow(Pi),ncol(Pi)) # This should be fixed to enable arbitrary levels
   self$np = apply(self$Pi, 2, function(a) return(length(unique(a))))
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

   switch (distr,
      "normal" = {
         distr = matrix(0,ns,4)
         distr[,1] = 1
         distr[,2] = 1
         colnames(distr) = c("family","variance","link","something")
      },
      "probit" = {
         distr = matrix(0,ns,4)
         distr[,1] = 2
         distr[,2] = 0
      }
   )
   self$distr = distr

   self$setPriors(setDefault=TRUE)
}

Hmsc$set("public", "setData", setData, overwrite=TRUE)

