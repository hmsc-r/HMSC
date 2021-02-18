#' @title predictLatentFactor
#'
#' @description Draws samples from the conditional predictive
#'     distribution of latent factors
#'
#' @param unitsPred a factor vector with random level units for which
#'     predictions are to be made
#' @param units a factor vector with random level units that are
#'     conditioned on
#' @param postEta a list containing samples of random factors at
#'     conditioned units
#' @param postAlpha a list containing samples of range (lengthscale)
#'     parameters for latent factors
#' @param rL a \code{HmscRandomLevel}-class object that describes the
#'     random level structure
#' @param predictMean a boolean flag indicating whether to return the
#'     mean of the predictive Gaussian process distribution
#' @param predictMeanField a boolean flag indicating whether to return
#'     the samples from the mean-field distribution of the predictive
#'     Gaussian process distribution
#'
#' @return a list of length \code{length(postEta)} containing samples
#'     of random factors at \code{unitsPred} from their predictive
#'     distribution conditional on the values at \code{units}
#'
#' @details Length of \code{units} vector and number of rows in
#'     \code{postEta} matrix shall be equal. The method assumes that
#'     the i-th row of \code{postEta} correspond to i-th element of
#'     \code{units}.
#'
#'   This method uses only the coordinates \code{rL$s} field of the
#'   \code{rL$s} argument. This field shall be a matrix with rownames
#'   covering the union of \code{unitsPred} and \code{units} factors.
#'
#'   In case of spatial random level, the computational complexity of
#'   the generic method scales cubically as the number of unobserved
#'   units to be predicted. Both \code{predictMean=TRUE} and
#'   \code{predictMeanField=TRUE} options decrease the asymptotic
#'   complexity to linear. The \code{predictMeanField=TRUE} option
#'   also preserves the uncertainty in marginal distribution of
#'   predicted latent factors, but neglects the inter-dependece
#'   between them.
#'
#' @importFrom stats rnorm dist
#' @importFrom methods is
#' @importFrom FNN knnx.index
#' @importFrom sp spDists
#'
#' @export

predictLatentFactor =
   function(unitsPred, units, postEta, postAlpha, rL, predictMean=FALSE,
            predictMeanField=FALSE)
{
   if(predictMean && predictMeanField)
      stop("predictMean and predictMeanField arguments cannot be simultaneously TRUE")
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
                  if (is(s1, "Spatial")) {
                     D11 <- spDists(s1)
                     D12 <- spDists(s1, s2)
                  } else {
                     dim = NCOL(s1)
                     D11 = as.matrix(dist(s1))
                     D12 = sqrt(Reduce("+",
                                       Map(function(i)
                                           outer(s1[,i], s2[,i], "-")^2,
                                           seq_len(dim))))
                  }
               } else {
                  ## s1, s2 are UNDEFINED: this will FAIL
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
               switch(rL$spatialMethod,
                      'Full' = {
                  unitsAll = c(units,unitsPred[indNew])
                  if(!is.null(rL$s)){
                     s = rL$s[unitsAll,]
                     if (is(s, "Spatial"))
                        D <- spDists(s)
                     else
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
               },
               'NNGP' = {
                  unitsAll = c(units,unitsPred[indNew])
                  s = rL$s[unitsAll,]
                  sOld = s[1:np,, drop=FALSE]
                  sNew = as.matrix(s[np+(1:nn),],, drop=FALSE)
                  indNN = knnx.index(sOld,sNew,k=rL$nNeighbours)
                  indices = list()
                  dist12 = matrix(NA,nrow=rL$nNeighbours,ncol=nn)
                  dist11 = array(NA, c(rL$nNeighbours,rL$nNeighbours,nn))
                  for(i in 1:nn){
                     ind = indNN[i,]
                     indices[[i]] = rbind(i*rep(1,length(ind)),ind)
                     if (is(sOld, "Spatial")) {
                        dist12[,i] <- spDists(sOld[ind,,drop=FALSE], sNew[i,])
                        dist11[,,i] = spDists(sOld[ind,, drop=FALSE])
                     } else {
                        das <- 0
                        for (dim in seq_len(rL$sDim))
                           das <- das + (sOld[ind, dim] - sNew[i, dim])^2
                        dist12[,i] <- sqrt(das)
                        dist11[,,i] <- as.matrix(dist(sOld[ind,]))
                     }
                  }
                  BgA = list()
                  FgA = list()
                  for(h in 1:nf){
                     if(alphapw[alpha[h],1] > 0){
                        K12 = exp(-dist12/alphapw[alpha[h],1])
                        ind1 = t(matrix(unlist(indices),nrow=2))[,2]
                        ind2 = t(matrix(unlist(indices),nrow=2))[,1]
                        K11 = exp(-dist11/alphapw[alpha[h],1])
                        K21iK11 = matrix(NA,ncol=nn,nrow=rL$nNeighbours)
                        for(i in 1:nn){
                           iK11 = solve(K11[,,i])
                           K21iK11[,i] = K12[,i]%*%iK11
                        }
                        B = Matrix(0,nrow=nn, ncol=np,sparse=TRUE)
                        B[cbind(ind2,ind1)] = as.vector(K21iK11)
                        Fmat = 1 - colSums(K21iK11*K12)
                        m = B%*%eta[,h]
                        etaPred[indNew,h] = as.numeric(m + sqrt(Fmat)*rnorm(nn))
                     } else{
                        etaPred[indNew,h] = rnorm(nn)
                     }
                  }
               },
               "GPP" = {
                  sKnot = rL$sKnot
                  unitsAll = c(units,unitsPred[indNew])
                  s = rL$s[unitsAll,]
                  if (is(s, "Spatial")) {
                     das <- spDists(s, sKnot)
                     dss <- spDists(sKnot)
                  } else {
                     dim = NCOL(s)
                     dss = as.matrix(dist(sKnot))
                     das = sqrt(Reduce("+",
                                          Map(function(i)
                                              outer(s[,i], sKnot[,i], "-")^2,
                                              seq_len(dim))))
                  }
                  dns = das[np+(1:nn),]
                  dnsOld = das[1:np,]
                  for(h in 1:nf){
                     ag = alpha[h]
                     if(alphapw[ag,1]>0){
                        Wns = exp(-dns/alphapw[ag,1])
                        W12 = exp(-dnsOld/alphapw[ag,1])
                        Wss = exp(-dss/alphapw[ag,1])
                        iWss= solve(Wss)
                        WnsiWss = Wns%*%iWss
                        dDn = 1 - rowSums(WnsiWss*Wns)
                        ## dDn can be numerically 0, but negative, say -2.2e-16
                        if (any(dDn < 0))
                           dDn[dDn < 0] <- 0
                        D =  W12%*%iWss%*%t(W12)
                        dD = 1-diag(D)
                        idD = 1/dD
                        tmp0 = matrix(rep(idD,NROW(sKnot)),ncol=NROW(sKnot))
                        idDW12 = tmp0*W12
                        FMat = Wss + t(W12)%*%idDW12
                        iF= solve(FMat)
                        LiF = chol(iF)
                        muS1 = iF%*%t(idDW12)%*%eta[,h]
                        epsS1 = LiF%*%rnorm(nrow(Wss))
                        m = Wns%*%(muS1+epsS1)
                        etaPred[indNew,h] = as.numeric(m + sqrt(dDn)*rnorm(nn))
                     } else{
                        etaPred[indNew,h] = rnorm(nn)
                     }
                  }
               })
            }
         }
      }
      postEtaPred[[pN]] = etaPred
   }
   return(postEtaPred)
}
