# id = diagonal of inverse residual variations

updateGammaEta = function(Z,V,id,Eta,Lambda,Alpha, X,Tr,Pi,rL, rLPar,Q,U){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   nc = ncol(X)
   nt = ncol(Tr)
   np = apply(Pi, 2, function(a) length(unique(a)))

   LRan = vector("list", nr)
   for(r in seq_len(nr)){
      LRan[[r]] = Eta[[r]][Pi[,r],,drop=FALSE]%*%Lambda[[r]]
   }

   Eta = vector("list", nr)
   iD = Diagonal(x=id)
   iDT = matrix(id,ns,nt)*Tr
   iD05T = matrix(sqrt(id),ns,nt)*Tr
   XtX = crossprod(X)
   iDT_XtX = kronecker(iDT,XtX)
   A = as.matrix(Matrix::tcrossprod(Matrix::tcrossprod(kronecker(Tr,Diagonal(nc)),chol(U)))) + kronecker(Q,V)
   iA = chol2inv(chol(A))

   for(r in seq_len(nr)){
      if(nr > 1){
         S = Z - Reduce("+", LRan[setdiff(1:nr,r)])
      } else{
         S = Z
      }
      lambda = Lambda[[r]]
      nf = nrow(lambda)
      lPi = Pi[,r]

      if(rL[[r]]$sDim == 0){
         cat("HMSC.updateGammaEta to be implemented for non-structured latent factors\n")
      } else{
         K = bdiag(lapply(seq_len(nf), function(x) rLPar[[r]]$Wg[,,Alpha[[r]][x]]))
         iK = bdiag(lapply(seq_len(nf), function(x) rLPar[[r]]$iWg[,,Alpha[[r]][x]]))
         P = sparseMatrix(i=1:ny,j=lPi)

         iD05Lamt = matrix(sqrt(id),ns,nf)*t(lambda)
         LamiD = lambda*matrix(id,nf,ns,byrow=TRUE)
         LamiD05 = t(iD05Lamt)
         LamiDLam = tcrossprod(lambda*matrix(sqrt(id),nf,ns,byrow=TRUE))
         PtX = Matrix::crossprod(P, X)
         PtP = Diagonal(x=Matrix::colSums(P))
         LamiDLam_PtP = kronecker(LamiDLam, PtP)
         LamiD_PtX = kronecker(LamiD, PtX)
         LamiDT_PtX = kronecker(LamiD%*%Tr,PtX)
         XtS = crossprod(X,S)

         W = iK + LamiDLam_PtP
         RW = chol(W)

         iLW.LamiDLam_PtP = backsolve(RW, LamiDLam_PtP, transpose=TRUE)
         iLW.LamiD_PtX = backsolve(RW, LamiD_PtX, transpose=TRUE)
         iDLamt_XtP.iW.LamiD_PtX = crossprod(iLW.LamiD_PtX)
         M = iA + kronecker(iD,XtX) - iDLamt_XtP.iW.LamiD_PtX
         RM = chol(M)

         mg10 = as.vector(XtS %*% iDT)
         mg21 = as.vector(Matrix::tcrossprod(Matrix::crossprod(P,S), LamiD))
         mg22 = backsolve(RW, backsolve(RW, mg21, transpose=TRUE))
         mg20 = Matrix::crossprod(LamiDT_PtX, mg22)
         mg31 = as.vector(XtS %*% iD) - Matrix::crossprod(LamiD_PtX, mg22)
         mg32 = backsolve(RM, backsolve(RM, mg31, transpose=TRUE))
         tmp1 = iDT_XtX - iDLamt_XtP.iW.LamiD_PtX %*% kronecker(Tr,Diagonal(nc))
         mg30 = Matrix::crossprod(tmp1, mg32)
         mg = U %*% (mg10 - mg20 - mg30)

         me10 = mg21
         me20 = LamiDLam_PtP %*% mg22
         tmp2 = LamiD_PtX - crossprod(iLW.LamiDLam_PtP, iLW.LamiD_PtX)
         me30 = tmp2 %*% mg32
         me = K %*% (me10 - me20 - me30)

         # F1g = kronecker(chol(crossprod(iD05T)), chol(crossprod(X))) %*% U
         # F1e = kronecker(chol(LamiDLam),sqrt(PtP)) %*% K
         F2g = backsolve(RW, LamiDT_PtX, transpose=TRUE) %*% U
         F2e = iLW.LamiDLam_PtP %*% K
         F3g = backsolve(RM, tmp1, transpose=TRUE) %*% U
         F3e = backsolve(RM, Matrix::t(tmp2), transpose=TRUE) %*% K

         # H1gg = crossprod(F1g)
         # H1ee = Matrix::crossprod(F1e)
         H1gg = U %*% kronecker(crossprod(iD05T), crossprod(X)) %*% U
         H1gg = (H1gg + t(H1gg)) / 2
         H1ee = K %*% kronecker(LamiDLam, PtP) %*% K
         H1ee = (H1ee + Matrix::t(H1ee)) / 2
         H1eg = K %*% kronecker(LamiD%*%Tr,PtX) %*% U
         H1 = rbind(cbind(H1gg,Matrix::t(H1eg)), cbind(H1eg,H1ee))
         H2 = Matrix::crossprod(cbind(F2g, F2e))
         H3 = Matrix::crossprod(cbind(F3g, F3e))

         m = as.vector(rbind(mg,me))
         G = bdiag(U,K) - (H1 - H2 - H3)
         gammaEta = m + crossprod(chol(G), rnorm(nc*nt+np*nf))
      }
      Gamma = matrix(gammaEta[1:(nc*nt)],nc,nt)
      Eta[[r]] = matrix(gammaEta[(nc*nt+1):(nc*nt+np[r]*nf)],np[r],nf)
      LRan[[r]] = Eta[[r]][Pi[,r],,drop=FALSE]%*%Lambda[[r]]
   }
   return(list(Gamma=Gamma, Eta=Eta))
}
