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
   indNew = !(indOld)
   n = length(unitsPred)
   np = length(units)
   nn = sum(indNew)

   postEtaPred = vector("list", predN)
   for(pN in 1:predN){
      eta = postEta[[pN]]
      nf = ncol(eta)
      etaPred = matrix(NA, n, nf)
      rownames(etaPred) = unitsPred
      etaPred[indOld,] = eta[match(unitsPred[indOld],units),]
      if(nn < 0){
         if(rL$sDim == 0){
            if(predictMean){
               etaPred[indNew,] = 0
            } else{
               etaPred[indNew,] = matrix(rnorm(sum(indNew)*nf), sum(indNew), nf)
            }
         } else{
            alpha = postAlpha[[pN]]
            alphapw = rL$alphapw
            unitsAll = c(units,unitsPred[indNew])
            s = rL$s[unitsAll,]
            distance = as.matrix(dist(s))
            for(h in 1:nf){
               K = exp(-distance/alphapw[alpha[h],1])
               K11 = K[1:np,1:np]
               K12 = K[1:np,np+(1:nn)]
               K22 = K[np+(1:nn),np+(1:nn)]
               iK11 = solve(K11)
               m = crossprod(K12, solve(K11, eta[,h]))
               W = K22 - crossprod(K12, solve(K11, K12))
               L = t(chol(W))
               etaPred[indNew,h] = m + L%*%rnorm(nn)
            }
         }
      }
      postEtaPred[[pN]] = etaPred
   }
   return(postEtaPred)
}

# Hmsc$set("public", "predictLatentFactor", predictLatentFactor, overwrite=TRUE)

