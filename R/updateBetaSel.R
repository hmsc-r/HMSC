#' @importFrom stats dnorm runif
#'
updateBetaSel = function(Z,XSelect, BetaSel, Beta, iSigma,
                         Lambda, Eta, Loff,X1,Pi,dfPi,rL){

   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   nc = ncol(X1[[1]])
   np = apply(Pi, 2, function(a) length(unique(a)))
   ncsel = length(XSelect)

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

   std = iSigma^-0.5

   X0 = list()
   for(j in 1:ns){
      X0[[j]] = matrix(0,nrow = ny, ncol = nc)
   }
   X = X1
   for (i in 1: ncsel){
      XSel = XSelect[[i]]
      for (spg in 1:length(XSel$q)){
         if(!BetaSel[[i]][spg]){
            fsp = which(XSel$spGroup==spg)
            for (j in fsp){
               X[[j]][,XSel$covGroup]=0
            }
         }
      }
   }
   LFix = matrix(NA,ny,ns)
   for(j in 1:ns)
      LFix[,j] = X[[j]]%*%Beta[,j]
   E = Reduce("+", c(list(LFix), LRan))
   if(!is.null(Loff)) E = E + Loff

   ll = matrix(NA,ny,ns)
   for (j in 1:ns){
      ll[,j]= dnorm(Z[,j], mean=E[,j], sd=std[j], log=TRUE)
   }

   BetaSelNew = BetaSel
   for (i in 1:ncsel){
      XSel = XSelect[[i]]
      for (spg in 1:length(XSel$q)){
         BetaSelNew[[i]][spg] = !(BetaSel[[i]][spg])
         fsp = which(XSel$spGroup==spg)

         X2 = X0
         for (j in fsp){
            X2[[j]][,XSel$covGroup]=X1[[j]][,XSel$covGroup]
         }
         LFix1 = matrix(0,ny,ns)
         for(j in fsp)
            LFix1[,j] = X2[[j]]%*%Beta[,j]
         if(BetaSelNew[[i]][spg]){
            ENew = E + LFix1
         } else {
            ENew = E - LFix1
         }

         llNew = ll
         for (j in fsp){
            llNew[,j]= dnorm(Z[,j], mean=ENew[,j], sd=std[j], log=TRUE)
         }
         lldif = sum(llNew[,fsp])-sum(ll[,fsp])
         q = XSel$q[spg]
         if(BetaSelNew[[i]][spg]){
            pridif = log(q)-log(1-q)
         } else {
            pridif = log(1-q)-log(q)
         }

         if(exp(lldif + pridif)>runif(1)){
            BetaSel[[i]][spg]=BetaSelNew[[i]][spg]
            E = ENew
            ll = llNew
         }
      }
   }

   X = X1
   for (i in 1:ncsel){
      XSel = XSelect[[i]]
      for (spg in 1:length(XSel$q)){
         if(!BetaSel[[i]][spg]){
            fsp = which(XSel$spGroup==spg)
            for (j in fsp){
               X[[j]][,XSel$covGroup]=0
            }
         }
      }
   }

   BetaSelXList = list()
   BetaSelXList$BetaSel = BetaSel
   BetaSelXList$X = X

   return(BetaSelXList)
}
