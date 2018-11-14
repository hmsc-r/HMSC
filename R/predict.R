#' @title predict.Hmsc
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'

predict.Hmsc = function(post=poolMcmcChains(self$postList), XData=NULL, X=NULL, dfPiNew=self$dfPi, rL=self$rL, Yc=NULL, mcmcStep=1, expected=FALSE, predictEtaMean=FALSE){
   if(!is.null(XData) && !is.null(X)){
      stop("Hmsc.predict: nly single of XData and X arguments can be specified")
   }
   if(!is.null(XData)){
      xlev = lapply(self$XData, levels)[unlist(lapply(self$XData, is.factor))]
      X = model.matrix(self$XFormula, XData, xlev=xlev)
   } else{
      if(is.null(X))
         X = self$X
   }
   if(!is.null(X)){
      if(nrow(Yc) != nrow(X)){
         stop("Hmsc.predict: number of rows in Yc and X must be equal")
      }
      if(ncol(Yc) != self$ns){
         stop("Hmsc.predict: number of columns in Yc must be equal to ns")
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
      LRan = vector("list", nr)
      Eta = vector("list",self$nr)
      for(r in seq_len(self$nr)){
         LRan[[r]] = predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
         Eta[[r]] = predPostEta[[r]][[pN]]
      }
      L = LFix + Reduce("+", LRan)

      if(!is.null(Yc) && any(!is.na(Yc))){
         Z = L
         Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,distr=self$distr)
         for(sN in seq_len(mcmcStep)){
            Eta = updateEta(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda,Alpha=sam$Alpha, rLPar=self$rLPar, X=X,Pi=PiNew,rL=rL)
            Z = updateZ(Y=Yc,Z=Z,Beta=sam$Beta,iSigma=iSigma,Eta=Eta,Lambda=sam$Lambda, X=X,Pi=PiNew,distr=self$distr)
         }
         for(r in seq_len(self$nr)){
            LRan[[r]] = predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
         }
         L = LFix + Reduce("+", LRan)
      }

      if(!expected){
         Z = L + matrix(sqrt(sam$sigma),nrow(L),self$ns,byrow=TRUE) * matrix(rnorm(nrow(L)*self$ns),nrow(L),self$ns)
      } else{
         Z = L
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

Hmsc$set("public", "predict", predict.Hmsc, overwrite=TRUE)

