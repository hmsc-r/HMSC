#' @title predictLatentFactor
#'
#' @description Computes the initial values before the sampling starts
#' @param initPar initial parameters value
#'
#' @examples
#'
#' @export

predictLatentFactor = function(unitsPred, units, postEta, postAlpha, rL, predictMean=FALSE){
   predN = length(postEta)
   indOld = (unitsPred %in% units)
   indNew = ~(indOld)
   N = length(unitsPred)
   postEtaPred = vector("list", predN)
   # np = length(unitsPred)
   for(pN in 1:predN){
      eta = postEta[[pN]]
      nf = ncol(eta)
      etaPred = matrix(NA, N, nf)
      rownames(etaPred) = unitsPred
      etaPred[indOld,] = eta[match(unitsPred[indOld],units),]
      if(rL$sDim == 0){
         if(predictMean){
            etaPred[indNew,] = 0
         } else{
            etaPred[indNew,] = matrix(rnorm(sum(indNew)*nf), sum(indNew), nf)
         }
      } else{
         alphapw = rL$alphapw
      }
      postEtaPred[[pN]] = etaPred
   }
   return(postEtaPred)
}

# Hmsc$set("public", "predictLatentFactor", predictLatentFactor, overwrite=TRUE)

