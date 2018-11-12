#' @title predictConditional
#'
#' @description Make conditions predictions given HMSC model object 
#'
#' @param Yc
#' @param post
#' @param mcmcStep
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


predictConditional = function(Yc, post=poolMcmcChains(self$postList), 
  mcmcStep=1, XData=NULL, X=NULL, dfPiNew=self$dfPi, rL=self$rL, 
  expected=FALSE, predictEtaMean=FALSE){
   if(!is.null(XData) && !is.null(X)){
      stop("Hmsc.predictConditional: only single of XData and X arguments can be specified")
   }
   if(!is.null(XData)){
      if(nrow(Yc) != nrow(XData)){
         stop("Hmsc.predictConditional: number of rows in Yc and XData must be equal")
      }
      xlev = lapply(self$XData, levels)[unlist(lapply(self$XData, is.factor))]
      X = model.matrix(self$XFormula, XData, xlev=xlev)
   } else{
      if(is.null(X))
         X = self$X
      if(nrow(Yc) != nrow(X)){
         stop("Hmsc.predictConditional: number of rows in Yc and X must be equal")
      }
   }


   predN = length(post)
   predPostEta = vector("list", self$nr)
   PiNew = matrix(NA,nrow(dfPiNew),self$nr)
   for(r in seq_len(self$nr)){
      postEta = lapply(post, function(c) c$Eta[[r]])
      postAlpha = lapply(post, function(c) c$Alpha[[r]])
      predPostEta[[r]] = predictLatentFactor(unitsPred=levels(dfPiNew[,r]),units=levels(self$dfPi[,r]),
         postEta=postEta,postAlpha=postAlpha,rL=rL[[r]],predictMean=predictEtaMean)
      rNames = rownames(predPostEta[[r]][[1]])
      PiNew[,r] = sapply(dfPiNew[,r], function(s) which(rNames==s))
   }
   pred = vector("list",predN)
   for(pN in 1:predN){
      sam = post[[pN]]
      LFix = X %*% sam$Beta
      LRan = matrix(0,nrow(LFix),ncol(LFix))
      Eta = vector("list",self$nr)
      for(r in seq_len(self$nr)){
         LRan = LRan + predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
         Eta[[r]] = predPostEta[[r]][[pN]]
      }
      rownames(LRan) = c()
      Z = LFix + LRan

      iSigma = sam$sigma^-1

      Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,distr=self$distr)
      for(sN in seq_len(mcmcStep)){
         Eta = updateEta(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda,Alpha=sam$Alpha, rLPar=self$rLPar, X=X,Pi=PiNew,rL=rL)
         Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,distr=self$distr)
      }

      for(j in 1:self$ns){
         if(self$distr[j,"family"] == 2){ # probit
            if(expected){
               Z[,j] = pnorm(Z[,j])
            } else{
               Z[,j] = as.numeric(Z[,j]>0)
            }
         }
         if(self$distr[j,"family"] == 3){ # poisson
            if(expected){
               Z[,j] = exp(Z[,j] + sam$sigma[j]/2)
            } else{
               Z[,j] = rpois(nrow(Z),exp(Z[,j]))
            }
         }
      }
      colnames(Z) = self$spNames
      pred[[pN]] = Z
   }
   return(pred)
}

Hmsc$set("public", "predictConditional", predictConditional, overwrite=TRUE)

