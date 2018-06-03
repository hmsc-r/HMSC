#' @title predict.Hmsc
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'

predict.Hmsc = function(post=poolMcmcChains(self$postList), XData=NULL, X=NULL, dfPiNew=self$dfPi, rL=self$rL, expected=FALSE, predictEtaMean=FALSE){
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

   predN = length(post)
   predPostEta = vector("list", self$nr)
   for(r in 1:self$nr){
      postEta = lapply(post, function(c) c$Eta[[r]])
      postAlpha = lapply(post, function(c) c$Alpha[[r]])
      predPostEta[[r]] = predictLatentFactor(unitsPred=levels(dfPiNew[,r]),units=levels(self$dfPi[,r]),
         postEta=postEta,postAlpha=postAlpha,rL=rL[[r]],predictMean=predictEtaMean)
   }
   pred = vector("list",predN)
   for(pN in 1:predN){
      sam = post[[pN]]
      LFix = X %*% sam$Beta
      LRan = matrix(0,nrow(LFix),ncol(LFix))
      for(r in 1:self$nr){
         LRan = LRan + predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
      }
      rownames(LRan) = c()
      L = LFix + LRan
      if(!expected){
         Z = L + matrix(rep(sam$sigma,nrow(L)), nrow(L), self$ns, byrow=TRUE)*matrix(rnorm(nrow(L)*self$ns), nrow(L), self$ns)
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
               Z[,j] = exp(Z[,j])
            } else{
               Z[,j] = rpois(nrow(Z)exp(,Z[,j]))
            }
         }
      }
      colnames(Z) = self$spNames
      pred[[pN]] = Z
   }
   return(pred)
}

Hmsc$set("public", "predict", predict.Hmsc, overwrite=TRUE)

