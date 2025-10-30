#' @title coralGetXDataExt
#'
#' @description Constructs extended \code{XData} and \code{XFormula} from fitted \code{Hmsc} for consequent CORAL analysis
#'
#' @param m fitted \code{Hmsc}-class object
#' @param nf upper number of leading latent factors to be used in CORAL analysis per \code{HmscRandomLevel}
#' @param varProp proportion of explanatory variance for selection of the number of leading latent factors
#'
#' @return
#' A named list of extended \code{XDataExt} dataframe and corresponding extended \code{XFormulaExt}
#'
#' @details
#' This functions finds a point estimate of latent factors from fitted \code{m} \code{Hmsc}-class object
#' and stacks \code{m$XData} and \code{m$XFormula} with this point estimate.
#' Such point estimate of latent factors is selected among the posterior MCMC samples of the \code{m}, specifically the one with top explanatory power.
#' The selected posterior sample is potentially truncated to a smaller number of leading latent factors, based on \code{nf} or \code{varProp} arguments.
#'
#' Only one of \code{nf} and \code{varProp} can be specified.
#' Each of these arguments can be provided as vector of length \code{m$nr} or scalar value that is broadcasted to such vector.
#' By default the \code{nf} is set to infinity, resulting in selection of all estimated latent factors.
#'
#' For example, if \code{Hmsc} model has \code{XFormula="~X1+X2"} and 2 random levels \code{LF1, LF2}, with 3 and 4 latent factors extracted correspondingly,
#' then the extended \code{XDataExt="~X1+X2 + LF1.1+LF1.2+LF1.3 + LF2.1+LF2.2+LF2.3+LF2.4"}.
#'
#' @importFrom stats as.formula
#'
#' @export

coralGetXDataExt = function(m, nf=NULL, varProp=NULL){
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

   mR2 = rowMeans(R2, na.rm=TRUE)
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
