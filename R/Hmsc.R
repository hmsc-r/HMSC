#' @title Hmsc
#'
#' @description Creates a Hmsc object
#'
#' @param Y matrix of species occurences or abundances
#' @param XFormula formula for linear regression component of HMSC
#' @param XData dataframe of measured covariates
#' @param X matrix of measured covariates for direct specification
#' @param XScale (boolean; default is TRUE) scale values in X matrix?
#' @param dfPi matrix of correspondence between sampling units and units on
#'  different levels of latent factors
#' @param rL list of HmscRandomLevel objects, specifying the structure and data for random levels
#' @param Xs
#' @param Xv
#' @param TrFormula
#' @param TrData
#' @param Tr
#' @param TrScale (boolean; default is TRUE) scale values in the trait (TR) matrix?
#' @param phyloTree
#' @param C
#' @param distr
#' @param spNames
#' @param trNames
#' @param covNames
#' @param levelNames
#'
#'
#' @return
#' Object of Hmsc class
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export


Hmsc = function(Y=NULL, XFormula=~., XData=NULL, X=NULL, XScale=TRUE,
   dfPi=NULL, rL=NULL, rNames=colnames(dfPi), Xs=NULL, Xv=NULL,
   TrFormula=NULL, TrData=NULL, Tr=NULL, TrScale=TRUE,
   phyloTree=NULL, C=NULL,
   distr="normal",
   priors=NULL){

   hM = structure(list(
            Y = NULL,
            XData=NULL, XFormula=NULL, X=NULL, XScaled=NULL, XInterceptInd=NULL,
            dfPi = NULL, rL=NULL, rNames=NULL, Pi=NULL,
            Xs = NULL,
            Xv = NULL,
            TrData=NULL, TrFormula=NULL, Tr=NULL, TrScaled=NULL, TrInterceptInd=NULL,
            C=NULL, phyloTree=NULL,
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
            postList=NULL), class="Hmsc")


   if(!is.matrix(Y)){
      stop("Hmsc.setData: Y argument must be a matrix of sampling units times species")
   }
   hM$Y = as.matrix(Y)
   hM$ny = nrow(Y)
   hM$ns = ncol(Y)
   if(is.null(colnames(hM$Y))){
      colnames(hM$Y) = sprintf(sprintf("sp%%.%dd",ceiling(log10(hM$ns))), 1:hM$ns)
   }
   hM$spNames = colnames(hM$Y)

   # linear regression covariates
   if(!xor(is.null(XData),is.null(X))){
      stop("Hmsc.setData: only single of XData and X arguments must be specified")
   }
   if(!is.null(XData)){
      if(nrow(XData) != hM$ny){
         stop("Hmsc.setData: the number of rows in XData should be equal to number of rows in Y")
      }
      if(any(is.na(XData))){
         stop("Hmsc.setData: XData parameter must not contain any NA values")
      }
      hM$XData = XData
      hM$XFormula = XFormula
      hM$X = model.matrix(XFormula, XData)
   }
   if(!is.null(X)){
      if(!is.matrix(X)){
         stop("Hmsc.setData: X must be a matrix")
      }
      if(nrow(X) != hM$ny){
         stop("Hmsc.setData: the number of rows in X should be equal to number of rows in Y")
      }
      if(any(is.na(X))){
         stop("Hmsc.setData: X parameter must not contain any NA values")
      }
      hM$XData = NULL
      hM$X = as.matrix(X)
   }
   hM$nc = ncol(hM$X)
   if(is.null(colnames(hM$X))){
      colnames(hM$X) = sprintf(sprintf("cov%%.%dd",ceiling(log10(hM$nc))), 1:hM$nc)
   }
   hM$covNames = colnames(hM$X)

   if(identical(XScale,FALSE)){
      hM$XScalePar = rbind(rep(0,hM$nc), rep(1,hM$nc))
      hM$XScaled = hM$X
      hM$XInterceptInd = NULL
   } else{
      XInterceptInd = which(colnames(hM$X) %in% c("Intercept","(Intercept)"))
      if(length(XInterceptInd)>1){
         stop("Hmsc.setData: only one column of X matrix could be named Intercept or (Intercept)")
      }
      if(!all(hM$X[,XInterceptInd] == 1)){
         stop("Hmsc.setData: intercept column in X matrix must be a column of ones")
      }
      if(length(XInterceptInd)==1){
         hM$XInterceptInd = XInterceptInd
      } else
         hM$XInterceptInd = NULL
      XScalePar = rbind(rep(0,hM$nc), rep(1,hM$nc))
      XScaled = hM$X
      if(identical(XScale,TRUE)){
         scaleInd = apply(hM$X, 2, function(a) !all(a %in% c(0,1)))
      } else{
         scaleInd = XScale
      }
      scaleInd[XInterceptInd] = FALSE
      if(length(XInterceptInd)>0){
         sc = scale(hM$X)
         XScalePar[,scaleInd] = rbind(attr(sc,"scaled:center"), attr(sc,"scaled:scale"))[,scaleInd]
      } else{
         sc = scale(hM$X, center=FALSE)
         XScalePar[2,scaleInd] = attr(sc,"scaled:scale")[scaleInd]
      }
      XScaled[,scaleInd] = sc[,scaleInd]
      hM$XScalePar = XScalePar
      hM$XScaled = XScaled
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
      if(nrow(TrData) != hM$ns){
         stop("Hmsc.setData: the number of rows in TrData should be equal to number of columns in Y")
      }
      if(any(is.na(TrData))){
         stop("Hmsc.setData: TrData parameter must not contain any NA values")
      }
      hM$TrData = TrData
      hM$TrFormula = TrFormula
      hM$Tr = model.matrix(TrFormula, TrData)
   }
   if(!is.null(Tr)){
      if(!is.matrix(Tr)){
         stop("Hmsc.setData: Tr must be a matrix")
      }
      if(nrow(Tr) != hM$ns){
         stop("Hmsc.setData: the number of rows in Tr should be equal to number of columns in Y")
      }
      if(any(is.na(Tr))){
         stop("Hmsc.setData: Tr parameter must not contain any NA values")
      }
      hM$TrData = NULL
      hM$Tr = Tr
   }
   if(is.null(hM$Tr)){
      hM$Tr = matrix(1,hM$ns,1)
   }
   hM$nt = ncol(hM$Tr)
   if(is.null(colnames(hM$Tr))){
      colnames(hM$Tr) = sprintf(sprintf("tr%%.%dd",ceiling(log10(hM$nt))), 1:hM$nt)
   }
   hM$traitNames = colnames(hM$Tr)

   if(identical(TrScale,FALSE)){
      hM$TrScalePar = rbind(rep(0,hM$nt), rep(1,hM$nt))
      hM$TrScaled = hM$Tr
      hM$TrInterceptInd = NULL
   } else{
      TrInterceptInd = which(colnames(hM$Tr) %in% c("Intercept","(Intercept)"))
      if(length(TrInterceptInd)>1){
         stop("Hmsc.setData: only one column of Tr matrix could be named Intercept or (Intercept)")
      }
      if(!all(hM$Tr[,TrInterceptInd]==1)){
         stop("Hmsc.setData: intercept column in Tr matrix must be a column of ones")
      }
      if(length(TrInterceptInd)==1){
         hM$TrInterceptInd = TrInterceptInd
      } else
         hM$TrInterceptInd = NULL
      TrScalePar = rbind(rep(0,hM$nt), rep(1,hM$nt))
      TrScaled = hM$Tr
      if(identical(TrScale,TRUE)){
         scaleInd = apply(hM$Tr, 2, function(a) !all(a %in% c(0,1)))
      } else{
         scaleInd = TrScale
      }
      scaleInd[TrInterceptInd] = FALSE
      if(length(TrInterceptInd)>0){
         sc = scale(hM$Tr)
         TrScalePar[,scaleInd] = rbind(attr(sc,"scaled:center"), attr(sc,"scaled:scale"))[,scaleInd]
      } else{
         sc = scale(hM$Tr, center=FALSE)
         TrScalePar[2,scaleInd] = attr(sc,"scaled:scale")[scaleInd]
      }
      TrScaled[,scaleInd] = sc[,scaleInd]
      hM$TrScalePar = TrScalePar
      hM$TrScaled = TrScaled
   }

   # phylogeny
   if(!is.null(C) && !is.null(phyloTree)){
      stop("Hmsc.setData: at maximum one of phyloTree and C arguments can be specified")
   }
   if(!is.null(phyloTree)){
      corM = vcv.phylo(phyloTree, model="Brownian", cor=T)
      corM = corM[hM$spNames,hM$spNames]
      hM$phyloTree = phyloTree
      hM$C = corM
   }
   if(!is.null(C)){
      if(any(dim(C) != hM$ns)){
         stop("Hmsc.setData: the size of square matrix C must be equal to number of species")
      }
      hM$C = C
   }

   # latent factors
   if(is.null(dfPi)){
      hM$dfPi = NULL
      hM$Pi = matrix(NA,hM$ny,0)
      hM$np = integer(0)
      hM$nr = 0
      hM$levelNames = character(0)
      if(!is.null(rL)){
         if(length(rL) > 0){
            stop("Hmsc.setData: dfPi is empty, but rL is not")
         }
      }
   } else {
      if(nrow(dfPi) != hM$ny){
         stop("Hmsc.setData: the number of rows in dfPi must be equal to number of rows in Y")
      }
      if(!all(rNames %in% names(rL)) || !all(rNames %in% colnames(dfPi))){
         stop("Hmsc.setData: colnames of dfPi must match names of rL")
      }
      hM$dfPi = dfPi
      hM$rL = rL[rNames]

      hM$Pi = matrix(NA,hM$ny,length(rNames),dimnames=list(NULL,rNames))
      for(r in 1:length(rNames))
         hM$Pi[,r] = as.numeric(dfPi[,rNames[r]])
      hM$np = apply(hM$Pi, 2, function(a) return(length(unique(a))))
      hM$nr = ncol(hM$Pi)
   }
   hM$levelNames = rNames

   switch (distr,
      "normal" = {
         distr = matrix(0,hM$ns,4)
         distr[,1] = 1
         distr[,2] = 1
      },
      "probit" = {
         distr = matrix(0,hM$ns,4)
         distr[,1] = 2
         distr[,2] = 0
      },
      "poisson" = {
         distr = matrix(0,hM$ns,4)
         distr[,1] = 3
         distr[,2] = 0
      },
      "lognormal poisson" = {
         distr = matrix(0,hM$ns,4)
         distr[,1] = 3
         distr[,2] = 1
      }
   )
   colnames(distr) = c("family","variance","link","something")
   hM$distr = distr

   hM = setPriors(hM, setDefault=TRUE)

   return(hM)
}

