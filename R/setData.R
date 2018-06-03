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

setData = function(Y=NULL, XFormula=~., XData=NULL, X=NULL, XScale=TRUE,
   dfPi=NULL, rL=NULL, Xs=NULL, Xv=NULL,
   TrFormula=NULL, TrData=NULL, Tr=NULL, TrScale=NULL,
   phyloTree=NULL, C=NULL,
   distr="normal", spNames=NULL,
   trNames=NULL, covNames=NULL, levelNames=NULL){

   if(!is.matrix(Y)){
      stop("Hmsc.setData: Y argument must be a matrix of sampling units times species")
   }
   self$Y = as.matrix(Y)
   self$ny = nrow(Y)
   self$ns = ncol(Y)
   if(is.null(spNames)){
      self$spNames = colnames(Y)
   } else{
      self$spNames = spNames
      colnames(self$Y) = spNames
   }

   # linear regression covariates
   if(!xor(is.null(XData),is.null(X))){
      stop("Hmsc.setData: only single of XData and X arguments must be specified")
   }
   if(!is.null(XData)){
      if(nrow(XData) != self$ny){
         stop("Hmsc.setData: the number of rows in XData should be equal to number of rows in Y")
      }
      self$XData = XData
      self$XFormula = XFormula
      self$X = model.matrix(XFormula, XData)
   }
   if(!is.null(X)){
      if(!is.matrix(X)){
         stop("Hmsc.setData: X must be a matrix")
      }
      if(nrow(X) != self$ny){
         stop("Hmsc.setData: the number of rows in X should be equal to number of rows in Y")
      }
      self$XData = NULL
      self$X = as.matrix(X)
   }
   self$nc = ncol(self$X)

   if(identical(XScale,FALSE)){
      self$XScalePar = NULL
      self$XScaled = self$X
      self$XScaleFlag = FALSE
      self$XInterceptInd = NULL
   } else{
      XInterceptInd = which(colnames(self$X) %in% c("Intercept","(Intercept)"))
      if(length(XInterceptInd)>1){
         stop("Hmsc.setData: only one column of X matrix could be named Intercept or (Intercept)")
      }
      if(!all(self$X[,XInterceptInd] == 1)){
         stop("Hmsc.setData: intercept column in X matrix must be a column of ones")
      }
      if(length(XInterceptInd)==1){
         self$XInterceptInd = XInterceptInd
      } else
         self$XInterceptInd = NULL
      XScalePar = matrix(0,2,self$nc)
      # XScalePar[1,XInterceptInd] = -1
      XScaled = self$X
      if(identical(XScale,TRUE)){
         scaleInd = apply(self$X, 2, function(a) !all(a %in% c(0,1)))
         self$XScaleFlag = TRUE
      } else{
         scaleInd = XScale
         self$XScaleFlag = NULL
      }
      scaleInd[XInterceptInd] = FALSE
      sc = scale(self$X)
      XScalePar[,scaleInd] = rbind(attr(sc,"scaled:center"), attr(sc,"scaled:scale"))[,scaleInd]
      if(length(XInterceptInd)>0){
         XScaled[,scaleInd] = sc[,scaleInd]
      } else{
         XScalePar[1,scaleInd] = 0
         XScaled[,scaleInd] = XScaled[,scaleInd] / matrix(XScalePar[2,scaleInd],ny,sum(scaleInd),byrow=T)
      }
      self$XScalePar = XScalePar
      self$XScaled = XScaled
   }

   # latent factors bindings
   if(is.null(dfPi)){
      self$dfPi = NULL
      self$Pi = NULL
      self$np = NULL
      self$nr = 0
   } else{
      if(nrow(dfPi) != self$ny){
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

   # traits
   if(!is.null(TrData)){
      if(!is.null(Tr)){
         stop("Hmsc.setData: at maximum one of TrData and Tr arguments can be specified")
      }
      if(is.null(TrFormula)){
         stop("Hmsc.setData: TrFormula argument must be specified if TrData is provided")
      }
   }
   if(!is.null(TrData)){
      if(nrow(TrData) != self$ns){
         stop("Hmsc.setData: the number of rows in TrData should be equal to number of columns in Y")
      }
      self$TrData = TrData
      self$TrFormula = TrFormula
      self$Tr = model.matrix(TrFormula, TrData)
   }
   if(!is.null(Tr)){
      if(!is.matrix(Tr)){
         stop("Hmsc.setData: Tr must be a matrix")
      }
      if(nrow(Tr) != self$ns){
         stop("Hmsc.setData: the number of rows in Tr should be equal to number of columns in Y")
      }
      self$TrData = NULL
      self$Tr = Tr
   }
   if(is.null(self$Tr)){
      self$Tr = matrix(1,self$ns,1)
   }
   self$nt = ncol(self$Tr)

   if(identical(TrScale,FALSE)){
      self$TrScalePar = NULL
      self$TrScaled = self$Tr
      self$TrScaleFlag = FALSE
      self$TrInterceptInd = NULL
   } else{
      TrInterceptInd = which(colnames(self$Tr) %in% c("Intercept","(Intercept)"))
      if(length(TrInterceptInd)>1){
         stop("Hmsc.setData: only one column of Tr matrix could be named Intercept or (Intercept)")
      }
      if(!all(self$Tr[,TrInterceptInd]==1)){
         stop("Hmsc.setData: intercept column in Tr matrix must be a column of ones")
      }
      if(length(TrInterceptInd)==1){
         self$TrInterceptInd = TrInterceptInd
      } else
         self$TrInterceptInd = NULL
      TrScalePar = matrix(0,2,self$nt)
      TrScaled = self$Tr
      if(identical(TrScale,TRUE)){
         scaleInd = apply(self$Tr, 2, function(a) !all(a %in% c(0,1)))
         self$TrScaleFlag = TRUE
      } else{
         scaleInd = TrScale
         self$TrScaleFlag = NULL
      }
      scaleInd[TrInterceptInd] = FALSE
      sc = scale(self$Tr)
      TrScalePar[,scaleInd] = rbind(attr(sc,"scaled:center"), attr(sc,"scaled:scale"))[,scaleInd]
      if(length(TrInterceptInd)>0){
         TrScaled[,scaleInd] = sc[,scaleInd]
      } else{
         TrScalePar[1,scaleInd] = 0
         TrScaled[,scaleInd] = TrScaled[,scaleInd] / matrix(TrScalePar[2,scaleInd],ns,sum(scaleInd),byrow=T)
      }
      self$TrScalePar = TrScalePar
      self$TrScaled = TrScaled
   }


   # phylogeny
   if(!is.null(C) && !is.null(phyloTree)){
      stop("Hmsc.setData: at maximum one of phyloTree and C arguments can be specified")
   }
   if(!is.null(phyloTree)){
      corM = vcv.phylo(phyloTree, model="Brownian", cor=T)
      corM = corM[self$spNames,self$spNames]
      self$phyloTree = phyloTree
      self$C = corM
   }
   if(!is.null(C)){
      if(any(dim(C) != self$ns)){
         stop("Hmsc.setData: the size of square matrix C must be equal to number of species")
      }
      self$C = C
   }

   if(is.null(covNames)){
      self$covNames = colnames(self$X)
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
      },
      "poisson" = {
         distr = matrix(0,self$ns,4)
         distr[,1] = 3
         distr[,2] = 0
      },
      "lognormal poisson" = {
         distr = matrix(0,self$ns,4)
         distr[,1] = 3
         distr[,2] = 1
      }
   )
   colnames(distr) = c("family","variance","link","something")
   self$distr = distr

   self$setPriors(setDefault=TRUE)
}

Hmsc$set("public", "setData", setData, overwrite=TRUE)

