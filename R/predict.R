#' @title predict
#'
#' @description Calculates predicted values from a fitted \code{Hmsc} model.
#'
#' @param object a fitted \code{Hmsc} model object
#' @param post a list of posterior samples of the HMSC model. By default uses all samples from the pooled
#' posterior of the hM object.
#' @param XData a dataframe specifying the unpreprocessed covariates for the predictions to be made.
#' Works only if the \code{XFormula} argument was specified in the \code{Hmsc()} model constructor call.
#' Requirements are similar to those in the \code{Hmsc} model constructor.
#' @param X a matrix specifying the covariates for the predictions to be made.
#' Only one of XData and X arguments may be provided.
#' @param XRRRData a dataframe of covariates for reduced-rank regression
#' @param XRRR a matrix of covariates for reduced-rank regression
#' @param studyDesign a matrix, specifying the structure of the study design for the prediction.
#' Requirements are similar to those of the \code{Hmsc} constructor. By default this argument is
#' assigned the study design of the training data in the fitted Hmsc model.
#' @param ranLevels a list of \code{HmscRandomLevel} objects, futher specifying the structure of
#' random levels. Requirements are similar to those of the \code{Hmsc} constructor.
#' Each level must cover all units, specified in the correspondingly named column of \code{studyDesign}
#' argument. By default this argument is assigned the list of \code{HmscRandomLevel} objects
#' specified for fitting Hmsc model.
#' @param Gradient an object returned by \code{\link{constructGradient}}. Providing \code{Gradient} is
#' an alternative for providing \code{XData}, \code{studyDesign} and \code{ranLevels}
#' @param Yc a matrix of the outcomes that are assumed to be known for conditional predictions.
#' @param mcmcStep the number of extra mcmc steps used for updating the random effects
#' @param expected boolean flag indicating whether to return the location parameter of the observation
#' models or sample the values from those.
#' @param predictEtaMean boolean flag indicating whether to use the estimated mean values of posterior
#' predictive distribution for random effets corresponding for the new units.
#' @param predictEtaMeanField boolean flag indicating whether to use draws from the mean-field of the
#' posterior predictive distribution for random effets corresponding for the new units.
#' @param \dots other arguments passed to functions.
#'
#' @details In \code{mcmcStep,the number of extra mcmc steps used for updating the random effects
#' for the Eta parameters, starting from the samples of the fitted Hmsc model in order to
#' account for the conditional infromation provided in the Yc argument. The higher this number is,
#' the more the obtained updated samples are unaffected by the posterior estimates of latent factors
#' in the model fitted to the training data and more resembles the true conditional posterior. However,
#' the elapsed time for conditional prediction grows approximately linearly as this parameter increases.
#' The exact number for sufficient is problem-dependent and should be assessed by e.g. gradually
#' increasing this parameter till the stationarity of the produced predictions.}
#'
#' @return A list of length \code{length(post)}, each element of which contains a sample from the posterior
#' predictive distribution (given the sample of the Hmsc model parameters in the corresponding element of
#' the \code{post} argument)
#'
#' @seealso \code{\link{predictLatentFactor}}
#'
#'
#' @importFrom stats model.matrix rnorm pnorm rpois
#'
#' @export

predict.Hmsc = function(object, post=poolMcmcChains(object$postList), XData=NULL,
                        X=NULL, XRRRData=NULL, XRRR=NULL, # this has to be updated to cov-dependent associations
                        studyDesign=object$studyDesign, ranLevels=object$ranLevels,
                        Gradient=NULL, Yc=NULL, mcmcStep=1, expected=FALSE,
                        predictEtaMean=FALSE, predictEtaMeanField=FALSE, ...)
{

   if(!is.null(Gradient)){
      XData=Gradient$XDataNew
      studyDesign=Gradient$studyDesignNew
      ranLevels=Gradient$rLNew
   }

   if(!is.null(XData) && !is.null(X)){
      stop("Hmsc.predict: only one of XData and X arguments can be specified")
   }
   if(!is.null(XRRRData) && !is.null(XRRR)){
      stop("Hmsc.predict: only one of XRRRData and XRRR arguments can be specified")
   }
   if(predictEtaMean==TRUE && predictEtaMeanField==TRUE)
      stop("Hmsc.predict: predictEtaMean and predictEtaMeanField arguments cannot be TRUE simultanuisly")
   if(!is.null(XData)){
      switch(class(XData)[1L],
             list={
                xlev = lapply(Reduce(rbind,object$XData), levels)[unlist(lapply(Reduce(rbind,object$XData), is.factor))]
                X = lapply(XData, function(a) model.matrix(object$XFormula, a, xlev=xlev))
             },
             data.frame={
                xlev = lapply(object$XData, levels)[unlist(lapply(object$XData, is.factor))]
                X = model.matrix(object$XFormula, XData, xlev=xlev)
             }
      )
   } else{
      if(is.null(X))
         X = object$X
   }
   if(!is.null(XRRRData)){
      xlev = lapply(object$XRRRData, levels)[unlist(lapply(object$XRRRData, is.factor))]
      XRRR = model.matrix(object$XRRRFormula, XRRRData, xlev=xlev)
   } else{
      if(is.null(object$ncRRR)) object$ncRRR=0
      if(is.null(XRRR) && object$ncRRR>0)
         XRRR=object$XRRR
   }
   switch(class(X)[1L],
          list={
             nyNew = nrow(X[[1]])
          },
          matrix={
             nyNew = nrow(X)
          }
   )

   if(!is.null(Yc)){
      if(ncol(Yc) != object$ns){
         stop("hMsc.predict: number of columns in Yc must be equal to ns")
      }
      if(nrow(Yc) != nyNew){
         stop("hMsc.predict: number of rows in Yc and X must be equal")
      }
   }
   if(!all(object$rLNames %in% colnames(studyDesign))){
      stop("hMsc.predict: dfPiNew does not contain all the necessary named columns")
   }
   if(!all(object$rLNames %in% names(ranLevels))){
      stop("hMsc.predict: rL does not contain all the necessary named levels")
   }

   if(!is.null(studyDesign)){
      dfPiNew = studyDesign[,object$rLNames,drop=FALSE]
   } else
      dfPiNew = matrix(NA,nyNew,0)
   rL = ranLevels[object$rLNames]

   predN = length(post)
   predPostEta = vector("list", object$nr)
   PiNew = matrix(NA,nrow(dfPiNew),object$nr)
   for(r in seq_len(object$nr)){
      postEta = lapply(post, function(c) c$Eta[[r]])
      postAlpha = lapply(post, function(c) c$Alpha[[r]])
      predPostEta[[r]] = predictLatentFactor(unitsPred=levels(dfPiNew[,r]),units=levels(object$dfPi[,r]),
                                             postEta=postEta,postAlpha=postAlpha,rL=rL[[r]],predictMean=predictEtaMean,predictMeanField=predictEtaMeanField)
      rowNames = rownames(predPostEta[[r]][[1]])
      PiNew[,r] = sapply(dfPiNew[,r], function(s) which(rowNames==s))
   }
   pred = vector("list",predN)
   for(pN in 1:predN){
      sam = post[[pN]]

      if(object$ncRRR>0){
         XB=XRRR%*%t(sam$wRRR)
      }
      switch(class(X)[1L],
             matrix = {
                X1=X
                if(object$ncRRR>0){
                   X1=cbind(X1,XB)
                }
                LFix = X1 %*% sam$Beta
             },
             list = {
                LFix = matrix(NA,nyNew,object$ns)
                for(j in 1:object$ns){
                   X1=X[[j]]
                   if(object$ncRRR>0){
                      X1=cbind(X1,XB)
                   }
                   LFix[,j] = X1%*%sam$Beta[,j]
                }

             }
      )
      LRan = vector("list",object$nr)
      Eta = vector("list",object$nr)
      for(r in seq_len(object$nr)){
         Eta[[r]] = predPostEta[[r]][[pN]]
         if(rL[[r]]$xDim == 0){
            LRan[[r]] = Eta[[r]][as.character(dfPiNew[,r]),] %*% sam$Lambda[[r]]
         } else{
            LRan[[r]] = matrix(0,object$ny,object$ns)
            for(k in 1:rL[[r]]$xDim)
               LRan[[r]] = LRan[[r]] + (Eta[[r]][as.character(dfPiNew[,r]),]*rL[[r]]$x[as.character(dfPiNew[,r]),k]) %*% sam$Lambda[[r]][,,k]
         }
      }
      if(object$nr > 0){L = LFix + Reduce("+", LRan)} else L = LFix

      if(!is.null(Yc) && any(!is.na(Yc))){
         Z = L
         Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=1/sam$sigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,dfPi=dfPiNew,distr=object$distr,rL=rL)
         for(sN in seq_len(mcmcStep)){
            Eta = updateEta(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=1/sam$sigma,Eta=Eta,Lambda=sam$Lambda,Alpha=sam$Alpha, rLPar=object$rLPar, X=X,Pi=PiNew,dfPi=dfPiNew,rL=rL)
            Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=1/sam$sigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,dfPi=dfPiNew,distr=object$distr,rL=rL)
         }
         for(r in seq_len(object$nr)){
            if(rL[[r]]$xDim == 0){
               LRan[[r]] = Eta[[r]][as.character(dfPiNew[,r]),] %*% sam$Lambda[[r]]
            } else{
               LRan[[r]] = matrix(0,object$ny,object$ns)
               for(k in 1:rL[[r]]$xDim)
                  LRan[[r]] = LRan[[r]] + (Eta[[r]][as.character(dfPiNew[,r]),]*rL[[r]]$x[as.character(dfPiNew[,r]),k]) %*% sam$Lambda[[r]][,,k]
            }
         }
         if(object$nr > 0){L = LFix + Reduce("+", LRan)} else L = LFix
      }
      if(!expected){
         Z = L + matrix(sqrt(sam$sigma),nrow(L),object$ns,byrow=TRUE) * matrix(rnorm(nrow(L)*object$ns),nrow(L),object$ns)
      } else{
         Z = L
      }
      for(j in 1:object$ns){
         if(object$distr[j,"family"] == 2){ # probit
            if(expected){
               Z[,j] = pnorm(Z[,j])
            } else{
               Z[,j] = as.numeric(Z[,j]>0)
            }
         }
         if(object$distr[j,"family"] == 3){ # poisson
            if(expected){
               Z[,j] = exp(Z[,j] + sam$sigma[j]/2)
            } else{
               Z[,j] = rpois(nrow(Z),exp(Z[,j]))
            }
         }
      }
      colnames(Z) = object$spNames

      for(i in 1:object$ns){
         m = object$YScalePar[1,i]
         s = object$YScalePar[2,i]
         if(m!=0 || s!=1){
            Z[,i] = Z[,i]*s + m
         }
      }
      pred[[pN]] = Z
   }
   return(pred)
}

