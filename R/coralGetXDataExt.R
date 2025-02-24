#' @title coralGetXDataExt
#'
#' @description Splits Hmsc data into common and rare parts for CORAL analysis
#'
#' @param m fitted \code{Hmsc}-class object
#' @param nf arg2
#' @param varProp arg3
#'
#' @return
#' A named list of extended XData data.frame and extended XFormula
#'
#' @importFrom stats as.formula
#'
#' @export

coralGetXDataExt = function(m, nf=NULL, varProp=NULL){
   # finds a point estimate of latent factors from fitted m.backbone.0 Hmsc object
   # takes either a vector of number of factors per random level to extract (can be a single value, broadcasted to all levels),
   # or a vector of values within [0,1], representing how big proportion of variance to ensure including per random level
   # (can be a single value, broadcasted to all levels)
   # by default takes all estimated factors?
   # returns a lsit of two objexts
   # XDataExt - a new XData by stacking the original XData and extracted LF estimates with appropriate
   # row projection of random level units to sampling units and column naming
   # XFormula - a new XFormula, e.g. in string format, which is cocnatenating the original XFormula with extra part of
   # all added latent factors
   # For example if Hmsc model has XFormula="~X1+X2" and 2 random levels, out of which 3 and 4 latent factors were extracted,
   # then the extended XFormula will be XFormula="~X1+X2 + LF1.1+LF1.2+LF1.3 + LF2.1+LF2.2+LF2.3+LF2.4"
   if(!is.null(nf) && !is.null(varProp)) stop("only one of ny and varProp arguments can be specified")
   if(is.null(nf) && is.null(varProp)) nf = Inf
   post = poolMcmcChains(m$postList)
   postN = length(post)
   R2 = matrix(NA, postN, m$ns)
   for(i in 1:postN){
      p = post[[i]]
      LFix = m$X %*% p$Beta
      LRan = vector("list", m$nr)
      for(r in seq_len(m$nr)){
         LRan[[r]] = p$Eta[[r]][m$Pi[,r],]%*%p$Lambda[[r]]
      }
      L = Reduce("+", c(list(LFix), LRan))
      pred = pnorm(L)
      R2[i,] = colSums(pred*(m$Y==1)) / colSums(m$Y==1) - colSums(pred*(m$Y==0)) / colSums(m$Y==0)
   }

   if(!is.null(varProp)){
      nf = rep(NA, m$nr)
      if(length(varProp) == 1) varProp = rep(varProp, m$nr)
   }
   if(length(nf) == 1) nf = rep(nf, m$nr)

   mR2 = rowMeans(R2)
   # plot(mR2, ylim=c(0,1))
   i = which.max(mR2) # isolate latent factor with highest explanatory power
   p = post[[i]]
   Eta = p$Eta
   Lambda = p$Lambda
   EtaUsed = vector("list", m$nr)
   EtaFormula = rep(NA, m$nr)
   for(r in seq_len(m$nr)){
      va = rep(NA, ncol(Eta[[r]]))
      for(f in seq_len(ncol(Eta[[r]]))){
         EtaLambdaSlice = Eta[[r]][m$Pi[,r],f,drop=FALSE] %*% Lambda[[r]][f,,drop=FALSE]
         va[f] = mean(apply(EtaLambdaSlice, 2, var))
      }
      #TODO can be rewritten without cycle, but for further generalization purpose probably better to leave as such
      # print(va)
      # print(apply(Eta[[r]][m$Pi[,r],,drop=FALSE], 2, var) * rowMeans(Lambda[[r]]^2))
      vp = va / sum(va)
      if(!is.null(varProp)){
         nf[r] = which.min(cumsum(va) >= varProp[r])
      } else{
         nf[r] = min(nf[r], ncol(Eta[[r]]))
      }
      # plot(cumsum(vp), main=sprintf("%s, nf=%d", colnames(m$rLNames[r]), nf[r]))
      # if(!is.null(varProp)) hline(varProp[r], lty=2)
      # plot(nf[r], cumsum(vp)[nf[r]], pch=19, col="red")
      EtaUsed[[r]] = Eta[[r]][m$Pi[,r], 1:nf[r], drop=FALSE]
      colnames(EtaUsed[[r]]) = paste(m$rLNames[r], ".f" , 1:nf[r], sep="")
      EtaFormula[r] = paste(colnames(EtaUsed[[r]]), collapse="+")
   }

   XDataExt = Reduce(cbind, c(list(m$XData), EtaUsed))
   XFormulaExt = as.formula(paste(c(deparse1(m$XFormula), EtaFormula), collapse=" + "))
   return(list(XDataExt=XDataExt, XFormulaExt=XFormulaExt))
}
