#' @title predictLatentFactor
#'
#' @description Draws samples from the conditional predictive distribution of latent factors
#'
#' @param unitsPred a factor vector with random level units for which predictions are to be made
#' @param units a factor vector with random level units that are conditioned on
#' @param postEta a list containing samples of random factors at conditioned units
#' @param postAlpha a list containing samples of range (lengthscale) parameters for latent factors
#' @param rL a \code{HmscRandomLevel}-class object that describes the random level structure
#' @param predictMean a boolean flag indicating whether to return the mean of the predictive Gaussian process
#'   distribution
#' @param predictMeanField a boolean flag indicating whether to return the samples from the mean-field distribution of
#'   the predictive Gaussian process distribution
#'
#' @return a list of length \code{length(postEta)} containing samples of random factors at \code{unitsPred} from their
#'   predictive distribution conditional on the values at \code{units}
#'
#' @details Length of \code{units} vector and number of rows in \code{postEta} matrix shall be equal. The method assumes
#'   that the i-th row of \code{postEta} correspond to i-th element of \code{units}.
#'
#'   This method uses only the coordinates \code{rL$s} field of the \code{rL$s} argument. This field shall be a matrix
#'   with rownames covering the union of \code{unitsPred} and \code{units} factors.
#'
#'   In case of spatial random level, the computational complexity of the generic method scales cubically as the number
#'   of unobserved units to be predicted. Both \code{predictMean=TRUE} and \code{predictMeanField=TRUE} options decrease
#'   the asymptotic complexity to linear. The \code{predictMeanField=TRUE} option also preserves the uncertainty in
#'   marginal distribution of predicted latent factors, but neglects the inter-dependece between them.
#'
#' @importFrom stats rnorm dist
#' @importFrom pdist pdist
#'
#' @export

predictLatentFactor = function(unitsPred, units, postEta, postAlpha, rL, predictMean=FALSE, predictMeanField=FALSE){
   if(predictMean==TRUE && predictMeanField==TRUE)
      stop("Hmsc.predictLatentFactor: predictMean and predictMeanField arguments cannot be simultaneously TRUE")
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
            if(predictMean || predictMeanField){
               if(!is.null(rL$s)){
                  s1 = rL$s[units,]
                  s2 = rL$s[unitsPred[indNew],]
                  D11 = as.matrix(dist(s1))
                  D12 = as.matrix(pdist(s1,s2))
               } else
               {
                  D11 = rL$distMat[s1,s1]
                  D12 = rL$distMat[s1,s2]
               }
               for(h in 1:nf){
                  if(alphapw[alpha[h],1] > 0){
                     K11 = exp(-D11/alphapw[alpha[h],1])
                     K12 = exp(-D12/alphapw[alpha[h],1])
                     m = crossprod(K12, solve(K11, eta[,h]))
                     if(predictMean)
                        etaPred[indNew,h] = m
                     if(predictMeanField){
                        LK11 = t(chol(K11))
                        iLK11K12 = solve(LK11, K12)
                        v = 1 - colSums(iLK11K12^2)
                        etaPred[indNew,h] = m + rnorm(nn, sd=sqrt(v))
                     }
                  } else{
                     if(predictMean)
                        etaPred[indNew,h] = 0
                     if(predictMeanField)
                        etaPred[indNew,h] = rnorm(nn)
                  }
               }
            } else{
               unitsAll = c(units,unitsPred[indNew])
               if(!is.null(rL$s)){
                  s = rL$s[unitsAll,]
                  D = as.matrix(dist(s))
               } else {
                  D = rL$distMat[unitsAll,unitsAll]
               }
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

