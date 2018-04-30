#' @title predict
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'
#' @export


predict.Hmsc = function(post, X=self$X, dfPiNew=self$dfPi, rL=self$rL, expected=FALSE){
   predN = length(post)
   nr = ncol(dfPiNew)
   predPostEta = vector("list", nr)
   for(r in 1:nr){
      postEta = lapply(post, function(c) c$Eta[[r]])
      postAlpha = lapply(post, function(c) c$Alpha[[r]])
      predPostEta[[r]] = predictLatentFactor(levels(dfPiNew[,r]),levels(self$dfPi[,r]),postEta,postAlpha,rL[[r]],predictMean=FALSE)
   }
   pred = vector("list",predN)
   for(pN in 1:predN){
      sam = post[[pN]]
      LFix = X %*% sam$Beta
      LRan = matrix(0,nrow(LFix),ncol(LFix))
      for(r in 1:nr){
         LRan = LRan + predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
      }
      L = LFix + LRan
      Z =
      if(!expected){
         Z = L + matrix(rep(sam$sigma,nrow(L)), nrow(L), self$ns, byrow=TRUE)
      }

      for(j in 1:self$ns){
         if(m$distr[j,"family"] == 2){
            if(expected){
               Z[,j] = pnorm(Z[,j])
            } else{
               Z[,j] = as.numeric(Z[,j]>0)
            }
         }
      }
      pred[[pN]] = Z
   }
   return(pred)
}

Hmsc$set("public", "predict", predict.Hmsc, overwrite=TRUE)

