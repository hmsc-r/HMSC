#' @title computeSAIR
#'
#' @description Computes the shared and idiosyncratic responses to measured and latent predictors
#'
#' @param hM a fitted \code{Hmsc} model object
#' @param X a design matrix to be used in the computations (as defaul hM$X)
#' 
#' @return
#' returns the posterior distribution of the parameters mu.X2,V.XX,mu.Omega2,V.OmegaOmega,mu.tot,V.tot,s
#' described Ovaskainen and Abrego (manuscript):
#' Measuring niche overlap with joint species distribution models: shared and idiosyncratic responses of the species to measured and latent predictors
#' 
#' @details
#' The shared and idiosyncratic responses are computed only for models without traits.
#'
#' @examples
#' # Simulate a small dataset, fit Hmsc model to it, compute SAIR, and show posterior means 
#'
#' nc = 2
#' ns = 5
#' ny = 10
#' mu = rnorm(n = nc)
#' X = matrix(rnorm(n=nc*ny),nrow=ny)
#' X[,1] = 1
#' eps = matrix(rnorm(nc*ns),nrow=nc)
#' L = matrix(rep(X%*%mu,ns),nrow=ny)  + X%*%eps
#' Y = pnorm(L)
#' m = Hmsc(Y = Y, XData = data.frame(env = X[,2]), distr = "probit")
#' m = sampleMcmc(m,samples=100,transient=50,verbose = 0)
#' SI = computeSAIR(m)
#' colMeans(SI)
#' 
#' @importFrom stats cov
#' @export


computeSAIR = function(hM, X = NULL){
  
  if(is.null(X)) X=hM$X

  if(hM$nt>1) return(NULL)
  
  if(hM$ncRRR==0){
    switch(class(X)[1L],
           matrix = {
             CX = cov(X)
           }, list = {
             CX = lapply(X, cov)
           }
    )
  }
  
  postList = poolMcmcChains(hM$postList)
  
  n = length(postList)
  SI = matrix(NA,ncol=7,nrow=n)
  colnames(SI) = c("mu.X2","V.XX","mu.Omega2","V.OmegaOmega","mu.tot","V.tot","s")
  for(i in 1:n){
    post = postList[[i]]
    
    if(hM$ncRRR>0){
      XB=hM$XRRR%*%t(post$wRRR)
      switch(class(X)[1L],
             matrix = {
               CX = cov(cbind(X,XB))
             }, list = {
               XX = X
               for(j in 1:length(XX)){
                 XX[[j]] = cbind(XX[[j]],XB)
               }
               CX = lapply(XX, cov)
             }
      )
    }    
    
    mu = post$Gamma
    V = post$V
    switch(class(X)[1L], matrix = {
      mu.X2 = t(mu)%*%CX%*%mu
      V.XX = sum(V*CX)
    }, list = {
      mu.X2 = 0
      V.XX = 0
      for(j in 1:hM$ns){
        mu.X2 =  mu.X2 + t(mu)%*%CX[[j]]%*%mu
        V.XX = V.XX + sum(V*CX[[j]])
      }
      mu.X2 =   mu.X2/hM$ns
      V.XX =   V.XX/hM$ns
    }
    )
    mu.Omega2 =  0
    V.OmegaOmega = 0
    if(hM$nr>0){
      for(h in 1:hM$nr){
        la = post$Lambda[[h]]
        Omega = t(la)%*%la
        Omega.diag = mean(diag(Omega))
        diag(Omega) = NA
        Omega.offdiag = mean(Omega, na.rm = TRUE)
        mu.Omega2 =  mu.Omega2 + Omega.offdiag
        V.OmegaOmega = V.OmegaOmega + (Omega.diag - Omega.offdiag)
      }
    }
    mu.tot = mu.X2 + mu.Omega2
    V.tot = V.XX + V.OmegaOmega
    s = mu.tot/(mu.tot+V.tot)
    SI[i,] = c(mu.X2,V.XX,mu.Omega2,V.OmegaOmega,mu.tot,V.tot,s)
  }
  return(SI)
}
