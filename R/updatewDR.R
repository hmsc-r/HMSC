#' @title updatewDR
#'
#' @description updates wDR
#'
updatewDR = function(Z,Beta,iSigma,Eta,Lambda,X1A,XDR,Pi,dfPi,rL,PsiDR,DeltaDR){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   if(class(X1A)=="matrix"){
       ncNDR = dim(X1A)[2]
   } else {
       ncNDR = dim(X1A[[1]])[2]
   }
   ncDR = length(DeltaDR)
   ncODR = dim(XDR)[2]

   BetaNDR = Beta[1:ncNDR,]
   BetaDR = Beta[-c(1:ncNDR),]

   switch(class(X1A),
          matrix = {
             LFix = X1A%*%BetaNDR
          },
          list = {
             LFix = matrix(NA,ny,ns)
             for(j in 1:ns)
                LFix[,j] = X1A[[j]]%*%BetaNDR[,j]
          }
   )
   LRan = vector("list", nr)
   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         LRan[[r]] = Eta[[r]][Pi[,r],]%*%Lambda[[r]]
      } else{
         LRan[[r]] = matrix(0,ny,ns)
         for(k in 1:rL[[r]]$xDim)
            LRan[[r]] = LRan[[r]] + (Eta[[r]][Pi[,r],]*rL[[r]]$x[as.character(dfPi[,r]),r]) %*% Lambda[[r]][,,r]
      }
   }
   if(nr > 1){
      S = Z - (LFix + Reduce("+", LRan))
   } else{
      S = Z - LFix
   }

   A1 = BetaDR%*%diag(iSigma)%*%t(BetaDR)
   A2 = t(XDR)%*%XDR
   QtiSigmaQ = kronecker(A2,A1)
   tauDR = matrix(apply(DeltaDR, 2, cumprod), ncDR, 1)
   tauMatDR = matrix(tauDR,ncDR,ncODR)
   iU=diag(as.vector(PsiDR*tauMatDR))+QtiSigmaQ
   RiU = chol(iU)
   U = chol2inv(RiU)
   mu1 = as.vector(BetaDR%*%diag(iSigma)%*%t(S)%*%XDR)
   mu = U %*% (mu1)
   we = mu + backsolve(RiU, rnorm(ncDR*ncODR))
   wDR = matrix(we,nrow = ncDR)

   X = X1A

   if(ncDR>0){
      XB=XDR%*%t(wDR)
      if(class(X)=="matrix"){
         X=cbind(X,XB)
      } else {
         for (j in 1:ns){
            X[[j]] = cbind(X[[j]],XB)
         }
      }
   }

   wDRXList=list()
   wDRXList$wDR = wDR
   wDRXList$X = X

   return(wDRXList)
}
