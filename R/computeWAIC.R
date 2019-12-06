#' @title computeWAIC
#'
#' @description Computes the value of WAIC (Widely Applicable Information Criterion) for the \code{Hmsc} model
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param ghN order of Gauss-Hermite quadrature for approximate numerical integration
#'
#' @details The result is exact for normal and probit observational models. For Poisson-type
#' observational model the result is obtained through numerical integration using Gauss-Hermite quadrature.
#'
#' @return the scalar WAIC
#'
#' @examples
#' # Compute WAIC of previously sampled Hmsc object
#' WAIC = computeWAIC(TD$m)
#'
#'
#' @importFrom stats dpois
#' @importFrom stats var
#' @importFrom abind abind
#' @importFrom statmod gauss.quad
#'
#' @export

computeWAIC = function(hM, ghN=11){
   post=poolMcmcChains(hM$postList)
   bind0 = function(...){
      abind(...,along=0)
   }
   Y = hM$Y
   Pi = hM$Pi
   dfPi = hM$dfPi
   distr = hM$distr
   X = hM$X
   rL = hM$rL
   ny = hM$ny
   ns = hM$ns
   nr = hM$nr
   np = hM$np

   indColNormal = (distr[,1]==1)
   indColProbit = (distr[,1]==2)
   indColPoisson = (distr[,1]==3)
   gq = gauss.quad(ghN, kind="hermite")
   gw = gq$weights
   gx = gq$nodes

   valList = vector("list", length(post))
   for(sN in 1:length(post)){
      Beta = post[[sN]]$Beta
      Eta = post[[sN]]$Eta
      Lambda = post[[sN]]$Lambda
      sigma = post[[sN]]$sigma
      switch(class(X)[1L],
             matrix = {
                LFix = X%*%Beta
             },
             list = {
                LFix = matrix(NA,ny,ns)
                for(j in 1:ns)
                   LFix[,j] = X[[j]]%*%Beta[,j]
             }
      )
      LRan = vector("list", nr)
      for(r in seq_len(nr)){
         if(rL[[r]]$xDim == 0){
            LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
         } else{
            LRan[[r]] = matrix(0,ny,ns)
            for(k in 1:rL[[r]]$xDim)
               LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),k]) %*% Lambda[[r]][,,k]
         }
      }
      if(nr > 0){
         E = LFix + Reduce("+", LRan)
      } else
         E = LFix

      indNA = is.na(Y)
      std = matrix(sigma^-0.5,ny,ns,byrow=TRUE)

      L = rep(0,ny)

      cN = sum(indColNormal)
      if(cN > 0){
         tmp = dnorm(x=Y[,indColNormal], mean=E[,indColNormal], sd=std[indColNormal], log=TRUE)
         tmp[is.na(Y[,indColNormal])] = 0
         if(cN > 1){
            L = L + rowSums(tmp)
         } else{
            L = L + tmp
         }
      }

      cN = sum(indColProbit)
      if(cN > 0){
         pz0 = pnorm(-E[,indColProbit], log.p=TRUE)
         pz1 = pnorm(E[,indColProbit], log.p=TRUE)
         tmp = pz1*Y[,indColProbit] + pz0*(1-Y[,indColProbit]) # this formula stands for unit std only, better to replace it with std-dependent one
         tmp[is.na(Y[,indColProbit])] = 0
         if(cN > 1){
            L = L + rowSums(tmp)
         } else{
            L = L + tmp
         }
      }

      cN = sum(indColPoisson)
      if(cN > 0){
         gX = array(E[,indColPoisson],c(ny,cN,ghN)) + sqrt(2)*array(gx,c(ny,cN,ghN))*array(std[indColPoisson],c(ny,cN,ghN))
         likeArray = dpois(Y, exp(gX))
         likeIntegral = log(sqrt(pi)^-1 * rowSums(likeArray*array(rep(gw,each=ny*cN),c(ny,cN,ghN)), dims=2))
         likeIntegral[is.na(Y[,indColPoisson])] = 0
         if(cN > 1){
            L = L + rowSums(likeIntegral)
         } else{
            L = L + likeIntegral
         }
      }
      valList[[sN]] = L
   }

   val = do.call(bind0, valList)
   Bl = -log(colMeans(exp(val)))
   V = rep(0,hM$ny)
   for (i in 1:hM$ny){
      V[i] = var(val[,i])
   }
   WAIC = mean(Bl + V)
   return(WAIC)
}
