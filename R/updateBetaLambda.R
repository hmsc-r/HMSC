# @title updateBetaLambda
#
# @description updates beta lambda
#
#' @importFrom stats rnorm
#' @importFrom Matrix bdiag Diagonal sparseMatrix tcrossprod forceSymmetric
#'
# NEW version
# Psi is the local continuous scale
# Delta is the column continuous scale
# Varphi is the local bernoulli scale (that could be fixed to one in case of non-SIS)
# Vartheta is the remaining part of the column scale that can be
# either bernoulli (CUSP) or product of delta (MGP, in this case Vartheta = NULL)

updateBetaLambda = function(Y,Z,Gamma,iV,iSigma,Eta,Psi,Delta,rho,Vartheta,Varphi, phyloPar, X,Tr,Pi,dfPi,C,rL, rhopw){
   ny = nrow(Z)
   ns = ncol(Z)
   nc = nrow(Gamma)
   nt = ncol(Tr)
   nr = ncol(Pi)
   S = Z
   Lambda = vector("list", nr)

   if (nr > 0) {
      EtaFull = vector("list", nr)
      nf = rep(NA, nr)
      ncr = rep(NA, nr)
      for (r in seq_len(nr)) {
         if (rL[[r]]$xDim == 0) {
            EtaFull[[r]] = Eta[[r]][Pi[, r], ]
         }
         else {
            EtaFull[[r]] = vector("list", rL[[r]]$xDim)
            for (k in 1:rL[[r]]$xDim) EtaFull[[r]][[k]] = Eta[[r]][Pi[,r], ] * rL[[r]]$x[as.character(dfPi[, r]),k]
         }
         nf[r] = ncol(Eta[[r]])
         ncr[r] = max(rL[[r]]$xDim, 1)
      }
      nfSum = sum(nf * ncr)
      EtaSt = matrix(unlist(EtaFull), ny, nfSum)

      LambdaPriorDiscT = vector("list", nr)
      for(r in seq_len(nr)){
         if(rL[[r]]$xDim == 0){
            LambdaPriorDiscT[[r]] = Varphi[[r]]
         } else { #GT: lets omit this more complicated case for now
            stop("Capacity for covariate-dependent associations is currently disabled")
            if(rL[[r]]$progShrinkType=="MGP"){
               LambdaPriorDiscT[[r]] = aperm(Varphi[[r]], c(2, 1, 3))
            } else if(rL[[r]]$progShrinkType=="CUSP"){
               LambdaPriorDiscT[[r]]
               LambdaPriorDiscT[[r]] = array(NA,dim=c(ns,nf[r],rL[[r]]$xDim))
               for(h in 1:nf[r]){
                  LambdaPriorDiscT[[r]][,h,] = Varphi[[r]][h,,]*Vartheta[[r]]
               }
            }
         }
      }

      # making a list of ns vectors of length nfSum
      # LambdaPriorDiscSt = lapply(seq_len(ns), function(i) matrix(unlist(LambdaPriorDiscT), ns, nfSum)[i,] )
      LambdaPriorDiscSt = split(do.call(rbind, LambdaPriorDiscT), rep(1:ns,each=nfSum))

      # Build XEtaList with Eta correctly adjusted
      switch(class(X)[1L], matrix = {
         XEtaList = lapply(LambdaPriorDiscSt, function(a) cbind(X, EtaSt*matrix(a,ny,nfSum,byrow=TRUE)))
      }, list = {
         XEtaList = mapply(function(X,Y){ cbind(X, EtaSt*matrix(Y,ny,nfSum,byrow=TRUE)) }, X=X, Y=LambdaPriorDiscSt)
      })
      PsiT = TauCt = vector("list", nr)
      for (r in seq_len(nr)) {
         if (rL[[r]]$xDim == 0) {
            PsiT[[r]] = t(Psi[[r]])
         }
         else {
            PsiT[[r]] = aperm(Psi[[r]], c(2,1,3))
         }
         # check if column specification is MGP or CUSP and compute the continuous part of Tau
         if(rL[[r]]$progShrinkType=="MGP"){
            TauCt[[r]] = apply(Delta[[r]], 2, cumprod)
         } else if(rL[[r]]$progShrinkType=="CUSP"){
            TauCt[[r]] = Delta[[r]] / Vartheta[[r]]
         }
      }
      psiSt = matrix(unlist(PsiT), nfSum, ns, byrow=TRUE)

      tauSt = matrix(unlist(TauCt), nfSum, 1)
      priorLambda = psiSt * matrix(tauSt, nfSum, ns)    # only the continuous part of lambda variance
   }
   else {
      nf = c()
      ncr = c()
      nfSum = 0
      switch(class(X)[1L], matrix = {
         XEtaList = rep(list(X),ns)
      }, list = {
         XEtaList = X
      })
      priorLambda = matrix(numeric(0), 0, ns)
   }
   Mu = rbind(tcrossprod(Gamma, Tr), matrix(0, nfSum, ns))
   Yx = !is.na(Y)
   S = Yx*S
   XEtaTXEtaList = vector("list", ns)
   isXTS = matrix(NA, nc+nfSum, ns)
   for(j in 1:ns){
      tmp = Yx[,j]*XEtaList[[j]]
      XEtaTXEtaList[[j]] = crossprod(tmp)
      isXTS[,j] = crossprod(XEtaList[[j]], S[,j]) * iSigma[j]
   }
   if (is.null(C)) {
      diagiV = diag(iV)
      P0 = matrix(0, nc + nfSum, nc + nfSum)
      P0[1:nc, 1:nc] = iV
      BetaLambda = matrix(NA, nc + nfSum, ns)
      for (j in 1:ns) {
         P = P0
         diag(P) = c(diagiV, priorLambda[, j])
         iU = P + XEtaTXEtaList[[j]] * iSigma[j]
         RiU = chol(iU)
         U = chol2inv(RiU)
         m = U %*% (P %*% Mu[, j] + isXTS[,j])
         BetaLambda[, j] = m + backsolve(RiU, rnorm(nc + nfSum))
      }
   }
   else { # available phylogeny information
      if(phyloPar$phyloFast == FALSE){
         if(length(rho)==1){
            # eiQ05 = (rhopw[rho,1]*eC + (1-rhopw[rho,1]))^-0.5
            # iQ = tcrossprod(VC*matrix(eiQ05,ns,ns,byrow=TRUE))
            # iQBar = kronecker(iV,iQ)
            iQBar = kronecker(iV, phyloPar$iQg[,,rho])
         } else{
            eiQ05List = lapply(as.list(rho), function(r) (rhopw[r,1]*eC + (1-rhopw[r,1]))^-0.5)
            VCeiQ05List = lapply(eiQ05List, function(eiQ05) VC*matrix(eiQ05,ns,ns,byrow=TRUE))
            tmp1 = bdiag(VCeiQ05List)
            LiV = t(chol(iV))
            tmp2 = kronecker(LiV,Diagonal(ns))
            tmp3 = tmp1 %*% tmp2
            iQBar = as.matrix(Matrix::tcrossprod(tmp3))
         }
         P = bdiag(iQBar, Diagonal(x=t(priorLambda)))
         tmp = vector("list", ns)
         for (j in 1:ns)
            tmp[[j]] = XEtaTXEtaList[[j]] * iSigma[j]
         tmpMat = do.call(rbind, tmp)
         ind1 = rep(rep(1:ns, each = nc + nfSum) + ns * rep(0:(nc+nfSum-1), ns), nc + nfSum)
         ind2 = rep(1:((nc + nfSum) * ns), each=nc+nfSum)
         mat = sparseMatrix(ind1, ind2, x = as.vector(tmpMat))
         RiU = Matrix::chol(mat + P)
         m1 = backsolve(RiU, P %*% as.vector(t(Mu)) + as.vector(t(isXTS)), transpose=TRUE)
         BetaLambda = matrix(backsolve(RiU, m1 + 0*rnorm(ns*(nc+nfSum))), nc+nfSum, ns, byrow=TRUE)
         # BetaLambda1 = BetaLambda
      } else{
         if(length(rho)==1){
            iQHat = phyloPar$iQHatList[[rho]]
            tipIndVec = phyloPar$iQHatIndList[[rho]]
            tmp1 = kronecker(iQHat, bdiag(iV,0*Diagonal(nfSum)))
            priorLambdaHat = matrix(0,nfSum,length(tipIndVec))
            priorLambdaHat[,which(tipIndVec>0)] = priorLambda[,tipIndVec[which(tipIndVec>0)]]
            priorLambdaHatExt = rbind(matrix(0,nc,length(tipIndVec)), priorLambdaHat)
            tmp2 = Diagonal(x=as.vector(priorLambdaHatExt))
            tmp3List = rep(list(matrix(0,nc+nfSum,nc+nfSum)),length(tipIndVec))
            P = tmp1 + tmp2
            for(i in which(tipIndVec>0)){
               tmp3List[[i]] = XEtaTXEtaList[[tipIndVec[i]]] * iSigma[tipIndVec[i]]
            }
            iUHat = P + bdiag(tmp3List)
            indVec = rep(tipIndVec>0,each=nc+nfSum)
            indAllZero = Matrix::diag(iUHat)==0
            P = P[!indAllZero,!indAllZero]
            iUHat = iUHat[!indAllZero,!indAllZero]
            indVec = indVec[!indAllZero]
            RiUHat = Matrix::chol(forceSymmetric(iUHat))
            MuHat = isXTSHat = matrix(0,nc+nfSum,length(tipIndVec))
            MuHat[,which(tipIndVec>0)] = Mu[,tipIndVec[which(tipIndVec>0)]]
            isXTSHat[,which(tipIndVec>0)] = isXTS[,tipIndVec[which(tipIndVec>0)]]
            m0Hat = P %*% as.vector(MuHat)[!indAllZero] + as.vector(isXTSHat)[!indAllZero]
            m1Hat = solve(t(RiUHat), m0Hat)
            m2Hat = solve(RiUHat, m1Hat + 0*rnorm(length(m1Hat)))
            BetaLambdaHat = matrix(0, nc+nfSum, length(tipIndVec))
            BetaLambdaHat[!indAllZero] = m2Hat
            BetaLambda = matrix(NA, nc+nfSum, ns)
            BetaLambda[,tipIndVec[which(tipIndVec>0)]] = m2Hat[indVec]#BetaLambdaHat[,which(tipIndVec>0)]
            # BetaLambda2 = BetaLambda
         } else{
            stop("Vector rho is not implemented with phyloFast option")
         }
      }
         # plot(BetaLambda1, BetaLambda2, main=rho)
   }
   Beta = BetaLambda[1:nc,,drop=FALSE]
   nfCumSum = c(0, cumsum(nf * ncr)) + nc
   LambdaTilde = vector("list", nr)
   for(r in seq_len(nr)) {
      if(rL[[r]]$xDim == 0) {
         LambdaTilde[[r]] = BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]), , drop = FALSE]
         # compute the real lambda by multiplying for the discrete scale
         if(rL[[r]]$progShrinkType=="MGP"){
            Lambda[[r]] = LambdaTilde[[r]]*Varphi[[r]]
         } else if(rL[[r]]$progShrinkType=="CUSP"){
            Lambda[[r]] = LambdaTilde[[r]]*Varphi[[r]]*Vartheta[[r]]
         }
      }
      else{
         LambdaTilde[[r]] = aperm(array(BetaLambda[(nfCumSum[r]+1):(nfCumSum[r+1]),,drop=FALSE], c(nf[r],ncr[r],ns)), c(1,3,2))
         if(rL[[r]]$progShrinkType=="MGP"){
            Lambda[[r]] = LambdaTilde[[r]]*Varphi[[r]]
         } else if(rL[[r]]$progShrinkType=="CUSP"){
            Varphitheta = Varphi[[r]]
            for(h in 1:nf[r]){
               Varphitheta[h,,] = Varphi[[r]][h,,]*Vartheta[[r]]
            }
            Lambda[[r]] = LambdaTilde[[r]]*Varphitheta
         }
      }
   }
   return(list(Beta=Beta, Lambda=Lambda, LambdaTilde=LambdaTilde))
}

