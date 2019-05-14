#' @title Hmsc$computeWAIC
#'
#' @description Computes WAIC (Widely Applicable Information Criterion) of the HMSC model
#'
#' @return the scalar WAIC
#'
#' @details implemented yet only for probit and normal models
#'
#' @examples
#'
#' WAIC = computeWAIC(hM=m)
#'
#' @export

computeWAIC = function(hM){
   post=poolMcmcChains(hM$postList)
   bind0 = function(...){
      abind(...,along=0)
   }
   valList = lapply(post, function(a) a$L)
   val = do.call(bind0, valList)
   Bl = -log(colMeans(exp(val)))
   V = rep(0,hM$ny)
   for (i in 1:hM$ny){
      V[i] = var(val[,i])
   }
   WAIC = mean(Bl + V)
   return(WAIC)
}
