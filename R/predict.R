#' @title predict.Hmsc
#'
#' @description Calculates predicted values for a given hMSC model
#'
#' @param post
#' @param XData
#' @param X
#' @param dfPiNew
#' @param rL
#' @param expected (boolean; default is FALSE)
#' @param predictEtaMean (boolean; default is FALSE)
#'
#'
#' @return
#'
#'
#' @seealso
#'
#'
#' @examples
#'
#' @export

predict.Hmsc = function(hM, post=poolMcmcChains(hM$postList), XData=NULL, X=NULL,
   studyDesign=hM$studyDesign, ranLevels=hM$ranLevels,
   Yc=NULL, mcmcStep=1, expected=FALSE, predictEtaMean=FALSE){

   if(!is.null(XData) && !is.null(X)){
      stop("Hmsc.predict: only single of XData and X arguments can be specified")
   }
   if(!is.null(XData)){
      switch(class(XData),
         list={
            xlev = lapply(Reduce(rbind,hM$XData), levels)[unlist(lapply(Reduce(rbind,hM$XData), is.factor))]
            X = lapply(XData, function(a) model.matrix(hM$XFormula, a, xlev=xlev))
         },
         data.frame={
            xlev = lapply(hM$XData, levels)[unlist(lapply(hM$XData, is.factor))]
            X = model.matrix(hM$XFormula, XData, xlev=xlev)
         }
      )
   } else{
      if(is.null(X))
         X = hM$X
   }
   switch(class(X),
      list={
         nyNew = nrow(X[[1]])
      },
      matrix={
         nyNew = nrow(X)
      }
   )

   if(!is.null(Yc)){
      if(ncol(Yc) != hM$ns){
         stop("hMsc.predict: number of columns in Yc must be equal to ns")
      }
      if(nrow(Yc) != nyNew){
         stop("hMsc.predict: number of rows in Yc and X must be equal")
      }
   }
   if(!all(hM$rLNames %in% colnames(studyDesign))){
      stop("hMsc.predict: dfPiNew does not contain all the necessary named columns")
   }
   if(!all(hM$rLNames %in% names(ranLevels))){
      stop("hMsc.predict: rL does not contain all the necessary named levels")
   }

   dfPiNew = studyDesign[,hM$rLNames,drop=FALSE]
   rL = ranLevels[hM$rLNames]

   predN = length(post)
   predPostEta = vector("list", hM$nr)
   PiNew = matrix(NA,nrow(dfPiNew),hM$nr)
   for(r in seq_len(hM$nr)){
      postEta = lapply(post, function(c) c$Eta[[r]])
      postAlpha = lapply(post, function(c) c$Alpha[[r]])
      predPostEta[[r]] = predictLatentFactor(unitsPred=levels(dfPiNew[,r]),units=levels(hM$dfPi[,r]),
         postEta=postEta,postAlpha=postAlpha,rL=rL[[r]],predictMean=predictEtaMean)
      rowNames = rownames(predPostEta[[r]][[1]])
      PiNew[,r] = sapply(dfPiNew[,r], function(s) which(rowNames==s))
   }
   pred = vector("list",predN)
   for(pN in 1:predN){
      sam = post[[pN]]
      switch(class(X),
         matrix = {
            LFix = X %*% sam$Beta
         },
         list = {
            LFix = matrix(NA,nyNew,hM$ns)
            for(j in 1:hM$ns)
               LFix[,j] = X[[j]]%*%sam$Beta[,j]
         }
      )
      LRan = vector("list",hM$nr)
      Eta = vector("list",hM$nr)
      for(r in seq_len(hM$nr)){
         LRan[[r]] = predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
         Eta[[r]] = predPostEta[[r]][[pN]]
      }
      if(hM$nr > 0){L = LFix + Reduce("+", LRan)} else L = LFix

      if(!is.null(Yc) && any(!is.na(Yc))){
         Z = L
         Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,distr=hM$distr)
         for(sN in seq_len(mcmcStep)){
            Eta = updateEta(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda,Alpha=sam$Alpha, rLPar=hM$rLPar, X=X,Pi=PiNew,rL=rL)
            Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,distr=hM$distr)
         }
         for(r in seq_len(hM$nr)){
            LRan[[r]] = predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
         }
         if(hM$nr > 0){L = LFix + Reduce("+", LRan)} else L = LFix
      }

      if(!expected){
         Z = L + matrix(sqrt(sam$sigma),nrow(L),hM$ns,byrow=TRUE) * matrix(rnorm(nrow(L)*hM$ns),nrow(L),hM$ns)
      } else{
         Z = L
      }

      for(j in 1:hM$ns){
         if(hM$distr[j,"family"] == 2){ # probit
            if(expected){
               Z[,j] = pnorm(Z[,j])
            } else{
               Z[,j] = as.numeric(Z[,j]>0)
            }
         }
         if(hM$distr[j,"family"] == 3){ # poisson
            if(expected){
               Z[,j] = exp(Z[,j] + sam$sigma[j]/2)
            } else{
               Z[,j] = rpois(nrow(Z),exp(Z[,j]))
            }
         }
      }
      colnames(Z) = hM$spNames
      pred[[pN]] = Z
   }
   return(pred)
}
