# d = diagonal of inverse residual variations

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
   iDTr = matrix(id,ns,nt) * Tr
   A = Matrix::tcrossprod(Matrix::tcrossprod(kronecker(Tr,Diagonal(nc)),chol(U))) + kronecker(Q,V)
   iA = chol2inv(chol(A))
   XtX = crossprod(X)

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
         LamiDLam = tcrossprod(lambda*matrix(sqrt(id),nf,ns,byrow=TRUE))
         LamiD = lambda*matrix(id,nf,ns,byrow=TRUE)
         P = sparseMatrix(i=1:ny,j=lPi)
         PtX = Matrix::crossprod(P, X)

         tmp1 = kronecker(LamiDLam, Diagonal(x=Matrix::colSums(P)))
         W = iK + tmp1
         RW = chol(W)
         LamiD_PtX = kronecker(LamiD,PtX)
         tmp2 = backsolve(RW, LamiD_PtX, transpose=TRUE)
         tmp3 = crossprod(tmp2)
         M = iA + kronecker(iD,XtX) - tmp3
         RM = chol(M)

         XtS = crossprod(X,S)
         mg10 = as.vector(XtS %*% iDTr)
         mg21 = as.vector(Matrix::tcrossprod(Matrix::crossprod(P,S), LamiD))
         mg22 = backsolve(RW, backsolve(RW, mg21, transpose=TRUE))
         mg20 = Matrix::crossprod(kronecker(LamiD%*%Tr,PtX), mg22)
         mg31 = as.vector(XtS %*% iD) - Matrix::crossprod(LamiD_PtX, mg22)
         mg32 = backsolve(RM, backsolve(RM, mg31, transpose=TRUE))
         mg30 = (kronecker(t(iDTr),XtX) - Matrix::crossprod(kronecker(Tr,Diagonal(nc)),tmp3)) %*% mg32
         mg = U %*% (mg10 - mg20 - mg30)

         me10 = mg21
         me20 = tmp1 %*% mg22
         # me31 = (LamiD_PtX - tmp1 %*% (backsolve(RW, backsolve(RW, LamiD_PtX, transpose=TRUE))))
         # me30 = me31 %*% mg32
         me30 = LamiD_PtX %*% mg32 - tmp1 %*% backsolve(RW, tmp2%*%mg32)
         me = K %*% (me10 - me20 - me30)
      }
      # Sigma = Matrix::crossprod(chol(A)%*%kronecker(Diagonal(ns),t(X)))
      Gamma = matrix(mg,nc,nt)
      Eta[[r]] = matrix(me,np[r],nf)
      LRan[[r]] = Eta[[r]][Pi[,r],,drop=FALSE]%*%Lambda[[r]]
   }
   return(list(Gamma=Gamma, Eta=Eta))
}
