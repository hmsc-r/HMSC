# @title updatewRRR
#
# @description updates wRRR
#
#' @importFrom stats rnorm
#
updatewRRR = function(Z,Beta,iSigma,Eta,Lambda,X1A,XRRR,Pi,dfPi,rL,PsiRRR,DeltaRRR){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   np = apply(Pi, 2, function(a) length(unique(a)))
   if(is.matrix(X1A)){
       ncNRRR = dim(X1A)[2]
   } else {
       ncNRRR = dim(X1A[[1]])[2]
   }
   ncRRR = length(DeltaRRR)
   ncORRR = dim(XRRR)[2]

   BetaNRRR = Beta[1:ncNRRR,]
   BetaRRR = matrix(Beta[-c(1:ncNRRR),],nrow=ncRRR)

   switch(class(X1A)[1L],
          matrix = {
             LFix = X1A%*%BetaNRRR
          },
          list = {
             LFix = matrix(NA,ny,ns)
             for(j in 1:ns)
                LFix[,j] = X1A[[j]]%*%BetaNRRR[,j]
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

   A1 = BetaRRR%*%diag(iSigma,nrow = length(iSigma))%*%t(BetaRRR)
   A2 = t(XRRR)%*%XRRR
   QtiSigmaQ = kronecker(A2,A1)
   tauRRR = matrix(apply(DeltaRRR, 2, cumprod), ncRRR, 1)
   tauMatRRR = matrix(tauRRR,ncRRR,ncORRR)
   iU=diag(as.vector(PsiRRR*tauMatRRR))+QtiSigmaQ
   RiU = chol(iU)
   U = chol2inv(RiU)
   mu1 = as.vector(BetaRRR%*%diag(iSigma,nrow = length(iSigma))%*%t(S)%*%XRRR)
   mu = U %*% (mu1)
   we = mu + backsolve(RiU, rnorm(ncRRR*ncORRR))
   wRRR = matrix(we,nrow = ncRRR)

   X = X1A

   if(ncRRR>0){
      XB=XRRR%*%t(wRRR)
      if(is.matrix(X)){
         X=cbind(X,XB)
      } else {
         for (j in 1:ns){
            X[[j]] = cbind(X[[j]],XB)
         }
      }
   }

   wRRRXList=list()
   wRRRXList$wRRR = wRRR
   wRRRXList$X = X

   return(wRRRXList)
}
