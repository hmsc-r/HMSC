#' @title predictLatentFactor
#'
#' @description Predicts
#'
#' @param unitsPred
#' @param units
#' @param postEta
#' @param postAlpha
#' @param rL
#' @param predictMean (boolean; default is FALSE)
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
      if(nn > 0){
         if(rL$sDim == 0){
            if(predictMean){
               etaPred[indNew,] = 0
            } else{
               etaPred[indNew,] = matrix(rnorm(sum(indNew)*nf), sum(indNew), nf)
            }
         } else{
            alpha = postAlpha[[pN]]
            alphapw = rL$alphapw
            if(predictMean){
               s1 = rL$s[units,]
               s2 = rL$s[unitsPred[indNew],]
               D11 = as.matrix(dist(s1))
               # replace following 3 lines with direct compuation of distances between s1 and s2
               unitsAll = c(units,unitsPred[indNew])
               s = rL$s[unitsAll,]
               D = as.matrix(dist(s))
               D12 = D[1:np,np+(1:nn)]
               for(h in 1:nf){
                  if(alphapw[alpha[h],1] > 0){
                     K11 = exp(-D11/alphapw[alpha[h],1])
                     K12 = exp(-D12/alphapw[alpha[h],1])
                     m = crossprod(K12, solve(K11, eta[,h]))
                     etaPred[indNew,h] = m
                  } else{
                     etaPred[indNew,h] = 0
                  }
               }
            } else{
               unitsAll = c(units,unitsPred[indNew])
               s = rL$s[unitsAll,]
               D = as.matrix(dist(s))
               for(h in 1:nf){
                  if(alphapw[alpha[h],1] > 0){
                     K = exp(-D/alphapw[alpha[h],1])
                     K11 = K[1:np,1:np]
                     K12 = K[1:np,np+(1:nn)]
                     K22 = K[np+(1:nn),np+(1:nn)]
                     m = crossprod(K12, solve(K11, eta[,h]))
                     W = K22 - crossprod(K12, solve(K11, K12))
                     L = t(chol(W))
                     etaPred[indNew,h] = m + L%*%rnorm(nn)
                  } else{
                     etaPred[indNew,h] = rnorm(nn)
                  }
               }
            }
         }
      }
      postEtaPred[[pN]] = etaPred
   }
   return(postEtaPred)
}

