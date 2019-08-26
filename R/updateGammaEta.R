### id = diagonal of inverse residual variations

#' @importFrom methods as
#' @importFrom stats rnorm
#' @importFrom Matrix Diagonal sparseMatrix bdiag
#'
updateGammaEta = function(Z,Gamma,V,iV,id,Eta,Lambda,Alpha, X,Tr,Pi,dfPi,rL, rLPar,Q,iQ,RQ,U,iU){
   ny = nrow(Z)
   ns = ncol(Z)
   nr = ncol(Pi)
   nc = ncol(X)
   nt = ncol(Tr)
   np = apply(Pi, 2, function(a) length(unique(a)))

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

   EtaNew = vector("list", nr)
   iD = Diagonal(x=id)
   iDT = matrix(id,ns,nt)*Tr
   iD05T = matrix(sqrt(id),ns,nt)*Tr
   XtX = crossprod(X)
   iDT_XtX = kronecker(iDT,XtX)
   A = as.matrix(Matrix::tcrossprod(Matrix::tcrossprod(kronecker(Tr,Diagonal(nc)),chol(U)))) + kronecker(Q,V)
   iA = chol2inv(chol(A))

   for(r in seq_len(nr)){
      if(rL[[r]]$xDim == 0){
         if(nr > 1){
            S = Z - Reduce("+", LRan[setdiff(1:nr,r)])
         } else{
            S = Z
         }
         lambda = Lambda[[r]]
         nf = nrow(lambda)
         lPi = Pi[,r]
         LamiD = lambda*matrix(id,nf,ns,byrow=TRUE)
         LamiDLam = tcrossprod(lambda*matrix(sqrt(id),nf,ns,byrow=TRUE))
         XtS = crossprod(X,S)

         if(rL[[r]]$sDim == 0){ # non-spatial level
            if(np[r] == ny){ # observation-corresponding LF
               W0 = LamiDLam + diag(nf)
               RW0 = chol(W0)
               iW0 = chol2inv(RW0)
               iLW0LamiD = backsolve(RW0,LamiD,transpose=TRUE)
               tmp1 = diag(id,ns) - crossprod(iLW0LamiD)
               M = iA + kronecker(tmp1, XtX)
               RM = chol(M)

               mb10 = as.vector(XtS * matrix(id,nc,ns,byrow=TRUE))
               mb20 = as.vector((tcrossprod(XtS, LamiD) %*% iW0) %*% LamiD)
               mb31 = backsolve(RM, backsolve(RM, mb10-mb20, transpose=TRUE))
               mb30 = kronecker(tmp1, XtX) %*% mb31
               mb = A %*% (mb10-mb20-mb30)
               Beta = matrix(mb + backsolve(RM,rnorm(nc*ns)),nc,ns)

               # update Gamma conditional on Beta
               R = chol(iU + kronecker(crossprod(backsolve(RQ,Tr,transpose=TRUE)), iV))
               mg = chol2inv(R) %*% as.vector((iV%*%Beta)%*%(iQ%*%Tr))
               GammaNew = matrix(mg + backsolve(R,rnorm(nc*nt)),nc,nt)

               # update Eta conditional on Beta, S
               S1 = S - X%*%Beta
               me = tcrossprod(S1,LamiD) %*% iW0
               EtaNew[[r]] = matrix(NA,ny,nf)
               EtaNew[[r]][lPi,] = me + t(backsolve(RW0,matrix(rnorm(ny*nf),nf,ny)))
            } else{ # non-observation-corresponding LF
               P = sparseMatrix(i=1:ny,j=lPi)
               PtX = Matrix::crossprod(P, X)
               colSumP = Matrix::colSums(P)
               PtP = Diagonal(x=colSumP)
               LamiDLam_PtP = kronecker(LamiDLam, PtP)
               LamiD_PtX = kronecker(LamiD, PtX)

               # W = Diagonal(nf*np[r]) + LamiDLam_PtP
               # RW = chol(W)
               # iW = Matrix::chol2inv(RW)

               WList = vector("list",np[r])
               RWList = vector("list",np[r])
               iWList = vector("list",np[r])
               LiWList = vector("list",np[r])
               for(p in 1:np[r]){
                  WList[[p]] = diag(nf) + colSumP[p]*LamiDLam
                  RWList[[p]] = chol(WList[[p]])
                  iWList[[p]] = chol2inv(RWList[[p]])
                  LiWList[[p]] = solve(RWList[[p]])
               }
               indR = rep((0:(nf-1))*np[r], nf*np[r]) + rep(rep(1:np[r],each=nf),nf)
               indC = rep(1:(nf*np[r]), each=nf)
               W = sparseMatrix(indR, indC, x=as.vector(Reduce(rbind, WList)))
               # RW = sparseMatrix(indR, indC, x=as.vector(Reduce(rbind, RWList)))
               iW = sparseMatrix(indR, indC, x=as.vector(Reduce(rbind, iWList)))
               LiW = sparseMatrix(indR, indC, x=as.vector(Reduce(rbind, LiWList)))

               # iLW.LamiD_PtX = Matrix::solve(t(RW), LamiD_PtX)
               # iLW.LamiD_PtX = backsolve(RW,LamiD_PtX,transpose=TRUE)
               iLW.LamiD_PtX = as.matrix(Matrix::crossprod(LiW, LamiD_PtX))
               iDLamt_XtP.iW.LamiD_PtX = crossprod(iLW.LamiD_PtX)
               tmp1 = kronecker(diag(id,ns),XtX) - iDLamt_XtP.iW.LamiD_PtX
               M = iA + tmp1
               RM = chol(M)

               mb10 = as.vector(XtS * matrix(id,nc,ns,byrow=TRUE))
               mb21 = as.vector(Matrix::tcrossprod(Matrix::crossprod(P,S), LamiD))
               # mb22 = as.vector(Matrix::solve(RW, Matrix::solve(t(RW),mb21)))
               mb22 = as.vector(iW %*% mb21)
               mb20 = as.vector(Matrix::crossprod(PtX,matrix(mb22,np[r],nf)) %*% LamiD)
               mb31 = backsolve(RM, backsolve(RM, mb10-mb20, transpose=TRUE))
               mb30 = tmp1 %*% mb31
               mb = A %*% (mb10-mb20-mb30)
               Beta = matrix(mb + backsolve(RM,rnorm(nc*ns)),nc,ns)

               # update Gamma conditional on Beta
               R = chol(iU + kronecker(crossprod(backsolve(RQ,Tr,transpose=TRUE)), iV))
               mg = chol2inv(R) %*% as.vector((iV%*%Beta)%*%(iQ%*%Tr))
               GammaNew = matrix(mg + backsolve(R,rnorm(nc*nt)),nc,nt)

               # update Eta conditional on Beta, S
               S1 = S - X%*%Beta
               PtS1 = Matrix::crossprod(P,S1)
               me10 = as.vector(Matrix::tcrossprod(PtS1,LamiD))
               # me21 = as.vector(Matrix::solve(RW, Matrix::solve(t(RW),me10)))
               me21 = as.vector(iW %*% me10)
               me20 = as.vector(PtP %*% matrix(me21,np[r],nf) %*% LamiDLam)
               me = me10 - me20
               # Eta[[r]] = matrix(me + backsolve(RW,rnorm(np[r]*nf)), np[r],nf)
               EtaNew[[r]] = matrix(me + LiW %*% rnorm(np[r]*nf), np[r],nf)
            }
         } else{ # spatial level
            P = sparseMatrix(i=1:ny,j=lPi)
            iD05Lamt = matrix(sqrt(id),ns,nf)*t(lambda)
            PtX = Matrix::crossprod(P, X)
            colSumP = Matrix::colSums(P)
            PtP = Diagonal(x=colSumP)
            LamiDLam_PtP = kronecker(LamiDLam, PtP)
            LamiD_PtX = kronecker(LamiD, PtX)
            LamiDT_PtX = kronecker(LamiD%*%Tr,PtX)
            switch(rL[[r]]$spatialMethod,
                   "Full" = {
                      K = bdiag(lapply(seq_len(nf), function(x) rLPar[[r]]$Wg[,,Alpha[[r]][x]]))
                      iK = bdiag(lapply(seq_len(nf), function(x) rLPar[[r]]$iWg[,,Alpha[[r]][x]]))
                   },
                   "NNGP" = {
                      stop("updataGammaEta: no method implemented yet for nearest neighbour gaussian process with GammaEta updater")
                   },
                   "GPP" = {
                      stop("updataGammaEta: no method implemented yet for gaussian predictive process with GammaEta updater")
                   }
            )
            W = iK + LamiDLam_PtP
            RW = chol(W)

            iLW.LamiD_PtX = as.matrix(backsolve(RW, LamiD_PtX, transpose=TRUE))
            iDLamt_XtP.iW.LamiD_PtX = crossprod(iLW.LamiD_PtX)
            M = iA + kronecker(diag(id,ns),XtX) - iDLamt_XtP.iW.LamiD_PtX
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
            me30 = LamiD_PtX %*% mg32 - LamiDLam_PtP %*% backsolve(RW, iLW.LamiD_PtX %*% mg32) # CHECK TYPE OF iLW.LamiD_PtX
            me = K %*% (me10 - me20 - me30)

            H = kronecker(iQ, iV) + as(kronecker(iD,XtX), "symmetricMatrix")
            iG1 = bdiag(iU, iK)
            iG2 = Matrix::crossprod(cbind(kronecker(iD05T,X), kronecker(iD05Lamt,P)))
            tmp = backsolve(Matrix::chol(H), cbind(iDT_XtX, Matrix::t(LamiD_PtX)), transpose=TRUE)
            iG3 = crossprod(tmp)
            if(dim(iG3)[1]==1){
               iG3 = crossprod(t(as.matrix(tmp)))
            }
            iG = iG1 + iG2 - iG3

            m = as.vector(rbind(mg,me))
            gammaEta = m + backsolve(chol(iG), rnorm(nc*nt+np[r]*nf))

            GammaNew = matrix(gammaEta[1:(nc*nt)],nc,nt)
            EtaNew[[r]] = matrix(gammaEta[(nc*nt+1):(nc*nt+np[r]*nf)],np[r],nf)
         }
         LRan[[r]] = EtaNew[[r]][Pi[,r],,drop=FALSE]%*%Lambda[[r]]
      } else{
         EtaNew[[r]] = Eta[[r]]
         GammaNew = Gamma
      }
   }
   return(list(Gamma=GammaNew, Eta=EtaNew))
}
