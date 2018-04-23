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

setData = function(Y=NULL, X=NULL, dfPi=NULL, rL=NULL, Xs=NULL, Xv=NULL, Tr=NULL, C=NULL, distr="normal", spNames=NULL,
   trNames=NULL, covNames=NULL, levelNames=NULL){
   self$Y = as.matrix(Y)
   self$ny = nrow(Y)
   self$ns = ncol(Y)

   if(nrow(X) != ny){
      stop("Hmsc.setData: the number of rows in X should be equal to number of rows in Y")
   }
   self$X = as.matrix(X)
   self$nc = ncol(X)

   if(is.null(dfPi)){
      self$dfPi = NULL
      self$Pi = NULL
      self$np = NULL
      self$nr = 0
   } else{
      if(nrow(dfPi) != ny){
         stop("Hmsc.setData: the number of rows in Pi should be equal to number of rows in Y")
      }
      self$dfPi = dfPi
      # tmp = as.data.frame(dfPi)
      self$Pi = matrix(NA,nrow(dfPi),ncol(dfPi)) # This should be fixed to enable arbitrary levels
      for(r in 1:ncol(dfPi))
         self$Pi[,r] = as.numeric(dfPi[,r])
      self$np = apply(self$Pi, 2, function(a) return(length(unique(a))))
      self$nr = ncol(dfPi)
      self$rL = rL
   }


   if(is.null(Tr)){
      self$Tr = matrix(1,self$ns,1)
   } else{
      if(nrow(Tr) != self$ns)
         stop("Hmsc.setData: the number of rows in Tr should be equal to number of species")
      self$Tr = Tr
   }
   self$nt = ncol(self$Tr)

   if(!is.null(C)){
      if(any(dim(C) != self$ns)){
         stop("Hmsc.setData: the size of square matrix C must be equal to number of species")
      }
      self$C = C
   }

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

   if(is.null(trNames)){
      self$trNames = colnames(Tr)
   } else{
      if(length(trNames))
         self$trNames = trNames
      colnames(self$Tr) = trNames
   }

   if(is.null(levelNames)){
      self$levelNames = colnames(dfPi)
   } else{
      self$levelNames = levelNames
      colnames(self$dfPi) = levelNames
      colnames(self$Pi) = levelNames
   }

   switch (distr,
      "normal" = {
         distr = matrix(0,self$ns,4)
         distr[,1] = 1
         distr[,2] = 1
      },
      "probit" = {
         distr = matrix(0,self$ns,4)
         distr[,1] = 2
         distr[,2] = 0
      }
   )
   colnames(distr) = c("family","variance","link","something")
   self$distr = distr

   self$setPriors(setDefault=TRUE)
}

Hmsc$set("public", "setData", setData, overwrite=TRUE)

