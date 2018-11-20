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

predict.Hmsc = function(hM, post=poolMcmcChains(hM$postList), XData=NULL, X=NULL, dfPiNew=hM$dfPi, rL=hM$rL, Yc=NULL, mcmcStep=1, expected=FALSE, predictEtaMean=FALSE){
   if(!is.null(XData) && !is.null(X)){
      stop("hMsc.predict: nly single of XData and X arguments can be specified")
   }
   if(!is.null(XData)){
      xlev = lapply(hM$XData, levels)[unlist(lapply(hM$XData, is.factor))]
      X = model.matrix(hM$XFormula, XData, xlev=xlev)
   } else{
      if(is.null(X))
         X = hM$X
   }
   if(!is.null(Yc)){
      if(ncol(Yc) != hM$ns){
         stop("hMsc.predict: number of columns in Yc must be equal to ns")
      }
      if(nrow(Yc) != nrow(X)){
         stop("hMsc.predict: number of rows in Yc and X must be equal")
      }
   }

   predN = length(post)
   predPostEta = vector("list", hM$nr)
   PiNew = matrix(NA,nrow(dfPiNew),hM$nr)
   for(r in seq_len(hM$nr)){
      postEta = lapply(post, function(c) c$Eta[[r]])
      postAlpha = lapply(post, function(c) c$Alpha[[r]])
      predPostEta[[r]] = predictLatentFactor(unitsPred=levels(dfPiNew[,r]),units=levels(hM$dfPi[,r]),
         postEta=postEta,postAlpha=postAlpha,rL=rL[[r]],predictMean=predictEtaMean)
      rNames = rownames(predPostEta[[r]][[1]])
      PiNew[,r] = sapply(dfPiNew[,r], function(s) which(rNames==s))
   }
   pred = vector("list",predN)
   for(pN in 1:predN){
      sam = post[[pN]]
      LFix = X %*% sam$Beta
      LRan = vector("list",hM$nr)
      Eta = vector("list",hM$nr)
      for(r in seq_len(hM$nr)){
         LRan[[r]] = predPostEta[[r]][[pN]][dfPiNew[,r],] %*% sam$Lambda[[r]]
         Eta[[r]] = predPostEta[[r]][[pN]]
      }
      L = LFix + Reduce("+", LRan)

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
         L = LFix + Reduce("+", LRan)
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
